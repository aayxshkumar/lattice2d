#pragma once

#include <cstdint>

#include "lattice/body.hpp"

namespace lattice {

struct ContactPoint {
  Vec2 point;          // world space, midway between the two surfaces
  double penetration;  // > 0 when overlapping
  // Which features touch (reference edge, incident edge or vertex), packed,
  // so contacts can be matched across steps. Unused until warm starting is
  // implemented.
  std::uint32_t feature = 0;
};

// Up to two contact points between two shapes. `normal` points from A to B.
struct Manifold {
  Vec2 normal;
  ContactPoint points[2];
  int count = 0;
};

// Narrow phase. Returns a manifold with count 0 if the shapes don't touch.
Manifold collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb);

Manifold collide_circles(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb);
Manifold collide_polygon_circle(const Shape& poly, const Transform& xp, const Shape& circle, const Transform& xc);
Manifold collide_polygons(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb);

}  // namespace lattice
