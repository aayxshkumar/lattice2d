#pragma once

#include <stdexcept>
#include <vector>

#include "lattice/math.hpp"

namespace lattice {

struct MassProperties {
  double mass = 0;
  double inertia = 0;  // about the centroid
  Vec2 centroid;
};

// A circle or a convex polygon, in the body's local frame. Polygons are
// recentred on their centroid when built, so a body's position is its centre
// of mass.
struct Shape {
  enum class Kind { Circle, Polygon };
  Kind kind = Kind::Circle;
  double radius = 0;
  std::vector<Vec2> vertices;  // counter-clockwise
  std::vector<Vec2> normals;   // outward; normals[i] belongs to edge i -> i+1

  static Shape circle(double r);
  static Shape box(double half_width, double half_height);
  // Convex polygon from points in any order (the convex hull is taken).
  // Throws std::invalid_argument for fewer than 3 non-collinear points.
  static Shape polygon(std::vector<Vec2> points);

  MassProperties mass(double density) const;
  AABB bounds(const Transform& xf) const;
};

constexpr std::size_t kMaxPolygonVertices = 16;

}  // namespace lattice
