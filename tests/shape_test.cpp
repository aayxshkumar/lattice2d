#include <gtest/gtest.h>

#include <cmath>

#include "lattice/shape.hpp"

using namespace lattice;

TEST(Shape, BoxMassAndInertia) {
  const Shape s = Shape::box(1.0, 0.5);  // 2 x 1
  const MassProperties m = s.mass(3.0);
  EXPECT_NEAR(m.mass, 6.0, 1e-12);
  // m (w^2 + h^2) / 12
  EXPECT_NEAR(m.inertia, 6.0 * (4.0 + 1.0) / 12.0, 1e-12);
  EXPECT_NEAR(m.centroid.length(), 0.0, 1e-12);
}

TEST(Shape, CircleMassAndInertia) {
  const MassProperties m = Shape::circle(2.0).mass(1.0);
  EXPECT_NEAR(m.mass, kPi * 4.0, 1e-12);
  EXPECT_NEAR(m.inertia, 0.5 * m.mass * 4.0, 1e-12);
}

TEST(Shape, PolygonIsRecentredOnItsCentroid) {
  // A right triangle: centroid at (1/3, 1/3) of the legs.
  const Shape t = Shape::polygon({{0, 0}, {3, 0}, {0, 3}});
  ASSERT_EQ(t.vertices.size(), 3u);
  Vec2 sum;
  for (const Vec2& v : t.vertices) sum += v;
  EXPECT_NEAR(sum.length(), 0.0, 1e-12);
  EXPECT_NEAR(t.mass(1.0).mass, 4.5, 1e-12);
  // Triangle inertia about the centroid: m (a^2 + b^2 + c^2) / 36.
  const double sides = 9 + 9 + 18;
  EXPECT_NEAR(t.mass(1.0).inertia, 4.5 * sides / 36.0, 1e-9);
}

TEST(Shape, HullDropsInteriorAndCollinearPoints) {
  const Shape s = Shape::polygon({{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}, {0, 1}});
  EXPECT_EQ(s.vertices.size(), 4u);
  // Counter-clockwise: positive signed area.
  double area2 = 0;
  for (std::size_t i = 0; i < s.vertices.size(); ++i)
    area2 += cross(s.vertices[i], s.vertices[(i + 1) % s.vertices.size()]);
  EXPECT_NEAR(area2, 8.0, 1e-12);
}

TEST(Shape, NormalsAreOutwardUnitVectors) {
  const Shape s = Shape::polygon({{0, 0}, {4, 0}, {5, 2}, {1, 3}});
  for (std::size_t i = 0; i < s.vertices.size(); ++i) {
    EXPECT_NEAR(s.normals[i].length(), 1.0, 1e-12);
    // Every vertex lies on the inner side of every face.
    for (const Vec2& v : s.vertices) EXPECT_LE(dot(s.normals[i], v - s.vertices[i]), 1e-12);
  }
}

TEST(Shape, RejectsDegenerateInput) {
  EXPECT_THROW(Shape::polygon({{0, 0}, {1, 1}}), std::invalid_argument);
  EXPECT_THROW(Shape::polygon({{0, 0}, {1, 1}, {2, 2}}), std::invalid_argument);
  EXPECT_THROW(Shape::circle(0), std::invalid_argument);
}

TEST(Shape, BoundsOfRotatedBox) {
  const Shape s = Shape::box(1, 1);
  const AABB b = s.bounds({{5, 0}, Rot::from_angle(kPi / 4)});
  EXPECT_NEAR(b.hi.x - b.lo.x, 2 * std::sqrt(2.0), 1e-12);
  EXPECT_NEAR(b.lo.x, 5 - std::sqrt(2.0), 1e-12);
}
