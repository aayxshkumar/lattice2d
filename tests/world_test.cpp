#include <gtest/gtest.h>

#include <cmath>

#include "lattice/scenes.hpp"
#include "lattice/world.hpp"

using namespace lattice;

TEST(World, StaticBodiesNeverMove) {
  Scene s = make_scene("drop");
  const Vec2 g0 = s.world.bodies()[0].position;
  for (int i = 0; i < 120; ++i) s.world.step(s.dt);
  EXPECT_EQ(s.world.bodies()[0].position, g0);
  EXPECT_EQ(s.world.bodies()[0].angle, 0.0);
}

TEST(World, HeadOnElasticCollisionSwapsVelocities) {
  World w(WorldConfig{{0, 0}, 10, 0.2, 0.005, 0.0});
  BodyDef a, b;
  a.shape = b.shape = Shape::circle(0.5);
  a.position = {0, 0};
  b.position = {1.05, 0};
  a.velocity = {3, 0};
  a.restitution = b.restitution = 1;
  a.friction = b.friction = 0;
  w.add(a);
  w.add(b);
  for (int i = 0; i < 10; ++i) w.step(1.0 / 60);
  EXPECT_NEAR(w.bodies()[0].velocity.x, 0, 1e-9);
  EXPECT_NEAR(w.bodies()[1].velocity.x, 3, 1e-9);
}

TEST(World, InelasticCollisionKeepsMomentumLosesEnergy) {
  World w(WorldConfig{{0, 0}, 10, 0.2, 0.005, 0.0});
  BodyDef a, b;
  a.shape = Shape::box(0.5, 0.5);
  b.shape = Shape::box(0.5, 0.5);
  a.position = {0, 0};
  b.position = {1.2, 0.1};
  a.velocity = {4, 0};
  b.density = 3;
  w.add(a);
  w.add(b);
  const Vec2 p0 = w.linear_momentum();
  const double e0 = w.kinetic_energy();
  for (int i = 0; i < 60; ++i) w.step(1.0 / 60);
  EXPECT_NEAR((w.linear_momentum() - p0).length(), 0.0, 1e-9);
  EXPECT_LT(w.kinetic_energy(), e0);
}

TEST(World, FrictionStopsASlidingBox) {
  Scene s = make_scene("drop");
  Body& box = s.world.bodies()[1];
  box.position = {0, 0.5};
  box.velocity = {3, 0};
  for (int i = 0; i < 180; ++i) s.world.step(s.dt);
  const Body& b = s.world.bodies()[1];
  EXPECT_LT(b.velocity.length(), 0.01);
  // v^2 / (2 mu g) with mu = 0.6: about 0.76 m.
  EXPECT_NEAR(b.position.x, 9.0 / (2 * 0.6 * 9.81), 0.08);
}

TEST(World, SameSceneSameHash) {
  Scene a = make_scene("rain", 100, 3), b = make_scene("rain", 100, 3), c = make_scene("rain", 100, 4);
  for (int i = 0; i < 200; ++i) {
    a.world.step(a.dt);
    b.world.step(b.dt);
    c.world.step(c.dt);
  }
  EXPECT_EQ(a.world.state_hash(), b.world.state_hash());
  EXPECT_NE(a.world.state_hash(), c.world.state_hash());
}

TEST(World, EveryNonGapCheckPasses) {
  for (const CheckResult& r : run_checks()) {
    if (r.known_gap) continue;
    EXPECT_TRUE(r.pass) << r.name << ": " << r.what << " = " << r.value << " (threshold " << r.threshold << ")";
  }
}
