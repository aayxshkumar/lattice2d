#include "lattice/world.hpp"

#include <algorithm>
#include <cstring>
#include <utility>

namespace lattice {

std::uint32_t World::add(const BodyDef& def) {
  Body b;
  b.id = static_cast<std::uint32_t>(bodies_.size());
  b.shape = def.shape;
  b.position = def.position;
  b.angle = def.angle;
  b.velocity = def.velocity;
  b.angular_velocity = def.angular_velocity;
  b.friction = def.friction;
  b.restitution = def.restitution;
  if (def.density > 0) {
    const MassProperties mp = b.shape.mass(def.density);
    b.mass = mp.mass;
    b.inv_mass = 1.0 / mp.mass;
    b.inertia = mp.inertia;
    b.inv_inertia = mp.inertia > 0 ? 1.0 / mp.inertia : 0;
  } else {
    // Static: never moves.
    b.velocity = {};
    b.angular_velocity = 0;
  }
  bodies_.push_back(std::move(b));
  return bodies_.back().id;
}

void World::find_contacts() {
  contacts_.clear();
  stats_ = {};
  // Sort and sweep along x. Ties and order are by body index, so the pair
  // order (and therefore the solver order) is deterministic.
  struct Entry {
    AABB box;
    std::uint32_t index;
  };
  std::vector<Entry> entries;
  entries.reserve(bodies_.size());
  std::vector<Transform> xf(bodies_.size());
  for (std::uint32_t i = 0; i < bodies_.size(); ++i) {
    xf[i] = bodies_[i].transform();
    entries.push_back({bodies_[i].shape.bounds(xf[i]), i});
  }
  std::sort(entries.begin(), entries.end(), [](const Entry& p, const Entry& q) {
    return p.box.lo.x < q.box.lo.x || (p.box.lo.x == q.box.lo.x && p.index < q.index);
  });
  std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    for (std::size_t j = i + 1; j < entries.size() && entries[j].box.lo.x <= entries[i].box.hi.x; ++j) {
      const std::uint32_t a = entries[i].index, b = entries[j].index;
      if (bodies_[a].is_static() && bodies_[b].is_static()) continue;
      if (!entries[i].box.overlaps(entries[j].box)) continue;
      pairs.emplace_back(std::min(a, b), std::max(a, b));
    }
  }
  std::sort(pairs.begin(), pairs.end());
  stats_.broadphase_pairs = pairs.size();

  for (const auto& [a, b] : pairs) {
    Manifold m = collide(bodies_[a].shape, xf[a], bodies_[b].shape, xf[b]);
    if (m.count == 0) continue;
    Contact c;
    c.a = a;
    c.b = b;
    c.manifold = m;
    c.friction = std::sqrt(bodies_[a].friction * bodies_[b].friction);
    c.restitution = std::max(bodies_[a].restitution, bodies_[b].restitution);
    contacts_.push_back(c);
    stats_.contact_points += static_cast<std::size_t>(m.count);
    for (int k = 0; k < m.count; ++k)
      stats_.max_penetration = std::max(stats_.max_penetration, m.points[k].penetration);
  }
  stats_.contacts = contacts_.size();
}

void World::prepare_contacts(double dt) {
  for (Contact& c : contacts_) {
    const Body& A = bodies_[c.a];
    const Body& B = bodies_[c.b];
    const Vec2 n = c.manifold.normal;
    const Vec2 t = cross(n, 1.0);
    for (int k = 0; k < c.manifold.count; ++k) {
      Contact::Point& p = c.points[k];
      const ContactPoint& cp = c.manifold.points[k];
      p.ra = cp.point - A.position;
      p.rb = cp.point - B.position;
      const double rna = cross(p.ra, n), rnb = cross(p.rb, n);
      const double kn = A.inv_mass + B.inv_mass + A.inv_inertia * rna * rna + B.inv_inertia * rnb * rnb;
      p.normal_mass = kn > 0 ? 1.0 / kn : 0;
      const double rta = cross(p.ra, t), rtb = cross(p.rb, t);
      const double kt = A.inv_mass + B.inv_mass + A.inv_inertia * rta * rta + B.inv_inertia * rtb * rtb;
      p.tangent_mass = kt > 0 ? 1.0 / kt : 0;
      p.normal_impulse = 0;
      p.tangent_impulse = 0;

      // Bounce only for real impacts, measured before solving.
      const Vec2 dv = B.velocity + cross(B.angular_velocity, p.rb) - A.velocity - cross(A.angular_velocity, p.ra);
      const double vn = dot(dv, n);
      p.velocity_bias = vn < -cfg_.restitution_threshold ? -c.restitution * vn : 0;
      p.position_bias = cfg_.baumgarte / dt * std::max(0.0, cp.penetration - cfg_.linear_slop);
    }
  }
}

void World::solve_velocities() {
  for (Contact& c : contacts_) {
    Body& A = bodies_[c.a];
    Body& B = bodies_[c.b];
    const Vec2 n = c.manifold.normal;
    const Vec2 t = cross(n, 1.0);
    for (int k = 0; k < c.manifold.count; ++k) {
      Contact::Point& p = c.points[k];

      // Normal: push apart until the closing speed reaches the target.
      Vec2 dv = B.velocity + cross(B.angular_velocity, p.rb) - A.velocity - cross(A.angular_velocity, p.ra);
      const double vn = dot(dv, n);
      const double target = std::max(p.velocity_bias, p.position_bias);
      double dpn = p.normal_mass * (target - vn);
      // Clamp the accumulated impulse, not the increment, so later
      // iterations can take back impulse applied earlier.
      const double old_n = p.normal_impulse;
      p.normal_impulse = std::max(old_n + dpn, 0.0);
      dpn = p.normal_impulse - old_n;
      const Vec2 pn = dpn * n;
      A.velocity -= A.inv_mass * pn;
      A.angular_velocity -= A.inv_inertia * cross(p.ra, pn);
      B.velocity += B.inv_mass * pn;
      B.angular_velocity += B.inv_inertia * cross(p.rb, pn);

      // Friction, bounded by the Coulomb cone of the current normal impulse.
      dv = B.velocity + cross(B.angular_velocity, p.rb) - A.velocity - cross(A.angular_velocity, p.ra);
      const double vt = dot(dv, t);
      double dpt = -p.tangent_mass * vt;
      const double max_t = c.friction * p.normal_impulse;
      const double old_t = p.tangent_impulse;
      p.tangent_impulse = std::clamp(old_t + dpt, -max_t, max_t);
      dpt = p.tangent_impulse - old_t;
      const Vec2 pt = dpt * t;
      A.velocity -= A.inv_mass * pt;
      A.angular_velocity -= A.inv_inertia * cross(p.ra, pt);
      B.velocity += B.inv_mass * pt;
      B.angular_velocity += B.inv_inertia * cross(p.rb, pt);
    }
  }
}

void World::step(double dt) {
  find_contacts();

  // Forces to velocities.
  for (Body& b : bodies_) {
    if (b.is_static()) continue;
    b.velocity += dt * (cfg_.gravity + b.inv_mass * b.force);
    b.angular_velocity += dt * b.inv_inertia * b.torque;
  }

  prepare_contacts(dt);
  for (int i = 0; i < cfg_.velocity_iterations; ++i) solve_velocities();

  // Velocities to positions.
  for (Body& b : bodies_) {
    if (b.is_static()) continue;
    b.position += dt * b.velocity;
    b.angle += dt * b.angular_velocity;
    b.force = {};
    b.torque = 0;
  }
  ++steps_;
}

double World::kinetic_energy() const {
  double e = 0;
  for (const Body& b : bodies_)
    if (!b.is_static()) e += b.kinetic_energy();
  return e;
}

Vec2 World::linear_momentum() const {
  Vec2 p;
  for (const Body& b : bodies_)
    if (!b.is_static()) p += b.mass * b.velocity;
  return p;
}

std::uint64_t World::state_hash() const {
  std::uint64_t h = 0xcbf29ce484222325ULL;
  auto mix = [&h](double d) {
    unsigned char bytes[sizeof d];
    std::memcpy(bytes, &d, sizeof d);
    for (unsigned char c : bytes) {
      h ^= c;
      h *= 0x100000001b3ULL;
    }
  };
  for (const Body& b : bodies_) {
    mix(b.position.x);
    mix(b.position.y);
    mix(b.angle);
    mix(b.velocity.x);
    mix(b.velocity.y);
    mix(b.angular_velocity);
  }
  return h;
}

}  // namespace lattice
