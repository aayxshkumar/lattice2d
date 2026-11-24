#include "lattice/scenes.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace lattice {

namespace {

constexpr double kBox = 0.5;  // half extent of the standard box (1 m cube)

void add_ground(World& w, double half_width = 40) {
  BodyDef g;
  g.shape = Shape::box(half_width, 0.5);
  g.position = {0, -0.5};  // top surface at y = 0
  g.density = 0;
  g.friction = 0.6;
  w.add(g);
}

BodyDef box_at(Vec2 p) {
  BodyDef d;
  d.shape = Shape::box(kBox, kBox);
  d.position = p;
  d.friction = 0.6;
  return d;
}

}  // namespace

std::vector<std::string> scene_names() {
  return {"drop", "stack", "pyramid", "ramp", "billiards", "bullet", "rain"};
}

Scene make_scene(const std::string& name, int n, std::uint64_t seed) {
  Scene s;
  s.name = name;
  World& w = s.world;
  Rng rng(seed);

  if (name == "drop") {
    add_ground(w);
    w.add(box_at({0, 4}));
    s.steps = 240;
  } else if (name == "stack") {
    const int h = n > 0 ? n : 10;
    add_ground(w);
    for (int i = 0; i < h; ++i) w.add(box_at({0, kBox + i * 2 * kBox}));
    s.steps = 600;
  } else if (name == "pyramid") {
    const int rows = n > 0 ? n : 10;
    add_ground(w);
    for (int row = 0; row < rows; ++row) {
      const int count = rows - row;
      const double x0 = -(count - 1) * kBox;
      for (int i = 0; i < count; ++i) w.add(box_at({x0 + i * 2 * kBox, kBox + row * 2 * kBox}));
    }
    s.steps = 600;
  } else if (name == "ramp") {
    const double deg = n > 0 ? n : 20;
    const double a = deg * kPi / 180.0;
    BodyDef slope;
    slope.shape = Shape::box(20, 0.5);
    slope.angle = a;
    slope.position = {0, 0};
    slope.density = 0;
    slope.friction = 0.5;
    w.add(slope);
    // Rest the box on the slope surface: up the slope normal by 0.5 + 0.5.
    BodyDef b = box_at({0, 0});
    b.angle = a;
    b.position = Vec2{-std::sin(a), std::cos(a)} * (0.5 + kBox);
    b.friction = 0.5;
    w.add(b);
    s.steps = 180;
  } else if (name == "billiards") {
    const int count = n > 0 ? n : 12;
    w = World(WorldConfig{{0, 0}, 10, 0.2, 0.005, 0.0});
    for (int i = 0; i < count; ++i) {
      BodyDef d;
      d.shape = Shape::circle(0.3);
      d.position = {(i % 4) * 1.5 + rng.uniform(-0.2, 0.2), (i / 4) * 1.5 + rng.uniform(-0.2, 0.2)};
      d.velocity = {rng.uniform(-3, 3), rng.uniform(-3, 3)};
      d.restitution = 1.0;
      d.friction = 0.0;
      w.add(d);
    }
    s.steps = 600;
  } else if (name == "bullet") {
    w = World(WorldConfig{{0, 0}, 10, 0.2, 0.005, 1.0});
    BodyDef wall;
    wall.shape = Shape::box(0.05, 3);
    wall.position = {5, 0};
    wall.density = 0;
    w.add(wall);
    BodyDef b;
    b.shape = Shape::circle(0.1);
    b.position = {0, 0};
    b.velocity = {300, 0};
    w.add(b);
    s.steps = 60;
  } else if (name == "rain") {
    const int count = n > 0 ? n : 400;
    add_ground(w, 15);
    for (int side : {-1, 1}) {
      BodyDef wall;
      wall.shape = Shape::box(0.5, 30);
      wall.position = {side * 15.5, 30};
      wall.density = 0;
      w.add(wall);
    }
    for (int i = 0; i < count; ++i) {
      BodyDef d;
      d.position = {rng.uniform(-13, 13), 2 + i * 0.12 + rng.uniform(0, 0.1)};
      d.angle = rng.uniform(0, 2 * kPi);
      const double r = rng.uniform(0.2, 0.45);
      switch (i % 3) {
        case 0: d.shape = Shape::circle(r); break;
        case 1: d.shape = Shape::box(r, r * rng.uniform(0.5, 1.0)); break;
        default: {
          std::vector<Vec2> pts;
          const int k = 3 + static_cast<int>(rng.next() % 5);
          for (int j = 0; j < k; ++j) {
            const double t = 2 * kPi * j / k + rng.uniform(0, 0.3);
            pts.push_back({r * std::cos(t), r * std::sin(t)});
          }
          d.shape = Shape::polygon(pts);
        }
      }
      w.add(d);
    }
    s.steps = 600;
  } else {
    throw std::invalid_argument("unknown scene '" + name + "'");
  }
  return s;
}

namespace {

CheckResult check(std::string name, std::string what, double value, double threshold, bool pass, bool gap = false) {
  return CheckResult{std::move(name), std::move(what), value, threshold, pass, gap};
}

void run(Scene& s, int steps) {
  for (int i = 0; i < steps; ++i) s.world.step(s.dt);
}

}  // namespace

std::vector<CheckResult> run_checks() {
  std::vector<CheckResult> out;

  {  // Free fall must match the closed form of semi-implicit Euler exactly.
    World w;
    BodyDef d = box_at({0, 100});
    w.add(d);
    const double dt = 1.0 / 60.0;
    const int n = 120;
    for (int i = 0; i < n; ++i) w.step(dt);
    const double g = w.config().gravity.y;
    const double expect = 100 + g * dt * dt * n * (n + 1) / 2.0;
    const double err = std::fabs(w.bodies()[0].position.y - expect);
    out.push_back(check("freefall", "|y - closed form| after 2 s (m)", err, 1e-9, err < 1e-9));
  }

  {  // A dropped box comes to rest on the ground, barely overlapping it.
    Scene s = make_scene("drop");
    run(s, 180);
    const Body& b = s.world.bodies()[1];
    const double speed = b.velocity.length() + std::fabs(b.angular_velocity);
    out.push_back(check("drop-rests", "speed + |spin| after 3 s", speed, 0.01, speed < 0.01));
    const double pen = s.world.last_step().max_penetration;
    out.push_back(check("drop-penetration", "resting penetration (m)", pen, 0.02, pen < 0.02));
    // A box dropped flat should stay level. Solving its two corner contacts
    // one after the other, without warm starting or a block solver, tips it.
    const double tilt = std::fabs(b.angle);
    out.push_back(check("drop-level", "|angle| of a box dropped flat, after 3 s (rad)", tilt, 1e-6, tilt < 1e-6, true));
  }

  {  // Momentum is conserved exactly by pairwise impulses.
    Scene s = make_scene("billiards");
    const Vec2 p0 = s.world.linear_momentum();
    const double e0 = s.world.kinetic_energy();
    run(s, s.steps);
    const Vec2 p1 = s.world.linear_momentum();
    const double dp = (p1 - p0).length() / std::max(1e-12, p0.length());
    out.push_back(check("billiards-momentum", "relative momentum change over 10 s", dp, 1e-9, dp < 1e-9));
    const double de = std::fabs(s.world.kinetic_energy() - e0) / e0;
    out.push_back(check("billiards-energy", "relative kinetic energy change over 10 s (elastic)", de, 0.02, de < 0.02));
  }

  {  // Static friction holds a box on a 20 degree slope (tan 20 = 0.36 < 0.5).
    Scene s = make_scene("ramp", 20);
    const Vec2 p0 = s.world.bodies()[1].position;
    run(s, s.steps);
    const double moved = (s.world.bodies()[1].position - p0).length();
    out.push_back(check("ramp-holds", "box displacement on a 20 deg slope, mu 0.5, 3 s (m)", moved, 0.05, moved < 0.05));
  }

  {  // On 35 degrees (tan 35 = 0.70 > 0.5) it slides with a = g (sin - mu cos).
    Scene s = make_scene("ramp", 35);
    const double a = 35 * kPi / 180.0;
    const Vec2 p0 = s.world.bodies()[1].position;
    run(s, 120);
    const double t = 120 * s.dt;
    const double moved = (s.world.bodies()[1].position - p0).length();
    const double expect = 0.5 * 9.81 * (std::sin(a) - 0.5 * std::cos(a)) * t * t;
    const double rel = std::fabs(moved - expect) / expect;
    out.push_back(check("ramp-slides", "relative error vs g(sin - mu cos) t^2/2 on 35 deg, 2 s", rel, 0.05, rel < 0.05));
  }

  // Towers and pyramids: measure how far the top box moves sideways. The
  // vertical sink isn't counted, since every resting contact is allowed
  // linear_slop of overlap, so a perfect 10-box tower still sinks ~5 cm.
  auto top_sideways = [](const std::string& scene, int n) {
    Scene s = make_scene(scene, n);
    const double x0 = s.world.bodies().back().position.x;
    run(s, 600);
    return std::fabs(s.world.bodies().back().position.x - x0);
  };
  {
    const double d = top_sideways("stack", 3);
    out.push_back(check("stack-3", "top box sideways drift after 10 s, 3 boxes (m)", d, 0.1, d < 0.1));
  }
  {  // These need warm starting: without it, 10 iterations don't converge
     // (a 10-box tower needs ~100, a 20-box tower ~300).
    const double d10 = top_sideways("stack", 10);
    out.push_back(check("stack-10", "top box sideways drift after 10 s, 10 boxes (m)", d10, 0.1, d10 < 0.1, true));
    const double d20 = top_sideways("stack", 20);
    out.push_back(check("stack-20", "top box sideways drift after 10 s, 20 boxes (m)", d20, 0.1, d20 < 0.1, true));
    const double p20 = top_sideways("pyramid", 20);
    out.push_back(check("pyramid-20", "top box sideways drift after 10 s, 20 rows (m)", p20, 0.1, p20 < 0.1, true));
  }

  {  // A fast, small body must hit a thin wall, not pass through it.
    Scene s = make_scene("bullet");
    run(s, s.steps);
    const double x = s.world.bodies()[1].position.x;
    out.push_back(check("bullet", "bullet x after 1 s (wall at x = 5; must stay below 5)", x, 5.0, x < 5.0, true));
  }

  {  // Same scene twice gives bit-identical state.
    Scene a = make_scene("rain", 200, 7), b = make_scene("rain", 200, 7);
    run(a, 300);
    run(b, 300);
    const bool same = a.world.state_hash() == b.world.state_hash();
    out.push_back(check("determinism", "rain(200) run twice: hashes differ (0 = identical)", same ? 0 : 1, 0.5, same));
  }

  return out;
}

}  // namespace lattice
