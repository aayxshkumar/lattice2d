#include <gtest/gtest.h>

#include <cmath>

#include "lattice/collision.hpp"

using namespace lattice;

namespace {
Transform at(double x, double y, double a = 0) { return {{x, y}, Rot::from_angle(a)}; }
}  // namespace

TEST(Collision, CirclesSeparatedAndOverlapping) {
  const Shape c = Shape::circle(1);
  EXPECT_EQ(collide(c, at(0, 0), c, at(2.01, 0)).count, 0);
  const Manifold m = collide(c, at(0, 0), c, at(1.5, 0));
  ASSERT_EQ(m.count, 1);
  EXPECT_NEAR(m.normal.x, 1, 1e-12);
  EXPECT_NEAR(m.points[0].penetration, 0.5, 1e-12);
  EXPECT_NEAR(m.points[0].point.x, 0.75, 1e-12);
}

TEST(Collision, BoxRestingOnBoxGivesTwoPoints) {
  const Shape ground = Shape::box(5, 0.5), box = Shape::box(0.5, 0.5);
  // Box sinks 0.01 into the ground's top face (y = 0).
  const Manifold m = collide(ground, at(0, -0.5), box, at(0, 0.49));
  ASSERT_EQ(m.count, 2);
  EXPECT_NEAR(m.normal.x, 0, 1e-12);
  EXPECT_NEAR(m.normal.y, 1, 1e-12);  // from ground to box
  for (int k = 0; k < 2; ++k) {
    EXPECT_NEAR(m.points[k].penetration, 0.01, 1e-9);
    EXPECT_NEAR(std::fabs(m.points[k].point.x), 0.5, 1e-9);  // the box's two bottom corners
  }
}

TEST(Collision, NormalAlwaysPointsFromAToB) {
  const Shape ground = Shape::box(5, 0.5), box = Shape::box(0.5, 0.5);
  const Manifold ab = collide(ground, at(0, -0.5), box, at(0.3, 0.45, 0.1));
  const Manifold ba = collide(box, at(0.3, 0.45, 0.1), ground, at(0, -0.5));
  ASSERT_GT(ab.count, 0);
  ASSERT_GT(ba.count, 0);
  EXPECT_GT(ab.normal.y, 0.9);
  EXPECT_LT(ba.normal.y, -0.9);
}

TEST(Collision, SeparatedPolygonsDontTouch) {
  const Shape a = Shape::box(1, 1), b = Shape::polygon({{0, 0}, {1, 0}, {0, 1}});
  EXPECT_EQ(collide(a, at(0, 0), b, at(3, 0)).count, 0);
  // Close but rotated away along a diagonal gap.
  EXPECT_EQ(collide(a, at(0, 0, kPi / 4), a, at(2.9, 0, kPi / 4)).count, 0);
}

TEST(Collision, TiltedBoxCornerOnGroundIsOnePoint) {
  const Shape ground = Shape::box(5, 0.5), box = Shape::box(0.5, 0.5);
  const double a = kPi / 6;
  // Lowest corner of the tilted box 0.02 below y = 0.
  const double half_diag_y = 0.5 * (std::cos(a) + std::sin(a));
  const Manifold m = collide(ground, at(0, -0.5), box, at(0, half_diag_y - 0.02, a));
  ASSERT_EQ(m.count, 1);
  EXPECT_NEAR(m.points[0].penetration, 0.02, 1e-9);
}

TEST(Collision, CircleOnPolygonFaceAndCorner) {
  const Shape box = Shape::box(1, 1), ball = Shape::circle(0.5);
  Manifold face = collide(box, at(0, 0), ball, at(0, 1.4));
  ASSERT_EQ(face.count, 1);
  EXPECT_NEAR(face.normal.y, 1, 1e-12);
  EXPECT_NEAR(face.points[0].penetration, 0.1, 1e-12);

  // Near the corner (1, 1), along the diagonal.
  const double d = 1 + 0.45 / std::sqrt(2.0);
  Manifold corner = collide(box, at(0, 0), ball, at(d, d));
  ASSERT_EQ(corner.count, 1);
  EXPECT_NEAR(corner.normal.x, 1 / std::sqrt(2.0), 1e-9);
  EXPECT_NEAR(corner.points[0].penetration, 0.05, 1e-9);
  EXPECT_EQ(collide(box, at(0, 0), ball, at(1.5, 1.5)).count, 0);

  // Same contact with the arguments swapped: normal flips.
  Manifold swapped = collide(ball, at(0, 1.4), box, at(0, 0));
  EXPECT_NEAR(swapped.normal.y, -1, 1e-12);
}

TEST(Collision, CircleCentreInsidePolygon) {
  const Shape box = Shape::box(1, 1), ball = Shape::circle(0.2);
  const Manifold m = collide(box, at(0, 0), ball, at(0.9, 0));
  ASSERT_EQ(m.count, 1);
  EXPECT_NEAR(m.normal.x, 1, 1e-12);
  EXPECT_NEAR(m.points[0].penetration, 0.3, 1e-12);
}
