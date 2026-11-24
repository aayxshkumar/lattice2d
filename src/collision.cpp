#include "lattice/collision.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace lattice {

namespace {

constexpr double kLinearSlop = 0.005;

Manifold flipped(Manifold m) {
  m.normal = -m.normal;
  return m;
}

}  // namespace

Manifold collide(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb) {
  using K = Shape::Kind;
  if (a.kind == K::Circle && b.kind == K::Circle) return collide_circles(a, xa, b, xb);
  if (a.kind == K::Polygon && b.kind == K::Circle) return collide_polygon_circle(a, xa, b, xb);
  if (a.kind == K::Circle && b.kind == K::Polygon) return flipped(collide_polygon_circle(b, xb, a, xa));
  return collide_polygons(a, xa, b, xb);
}

Manifold collide_circles(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb) {
  Manifold m;
  const Vec2 d = xb.p - xa.p;
  const double dist2 = d.length_squared();
  const double r = a.radius + b.radius;
  if (dist2 > r * r) return m;
  const double dist = std::sqrt(dist2);
  m.normal = dist > 0 ? d * (1.0 / dist) : Vec2{0, 1};
  const Vec2 pa = xa.p + a.radius * m.normal;
  const Vec2 pb = xb.p - b.radius * m.normal;
  m.points[0] = {0.5 * (pa + pb), r - dist, 0};
  m.count = 1;
  return m;
}

Manifold collide_polygon_circle(const Shape& poly, const Transform& xp, const Shape& circle, const Transform& xc) {
  Manifold m;
  // Work in the polygon's frame.
  const Vec2 c = xp.apply_inverse(xc.p);
  const double r = circle.radius;
  const std::size_t n = poly.vertices.size();
  double best = -std::numeric_limits<double>::infinity();
  std::size_t edge = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double s = dot(poly.normals[i], c - poly.vertices[i]);
    if (s > r) return m;  // separating axis
    if (s > best) {
      best = s;
      edge = i;
    }
  }
  const Vec2 v1 = poly.vertices[edge], v2 = poly.vertices[(edge + 1) % n];
  Vec2 normal_local, closest;
  if (best < std::numeric_limits<double>::epsilon()) {
    // Centre inside the polygon: push out along the nearest face.
    normal_local = poly.normals[edge];
    closest = c - best * normal_local;
  } else {
    const double u1 = dot(c - v1, v2 - v1);
    const double u2 = dot(c - v2, v1 - v2);
    if (u1 <= 0) {
      if ((c - v1).length_squared() > r * r) return m;
      normal_local = (c - v1).normalized();
      closest = v1;
    } else if (u2 <= 0) {
      if ((c - v2).length_squared() > r * r) return m;
      normal_local = (c - v2).normalized();
      closest = v2;
    } else {
      normal_local = poly.normals[edge];
      closest = c - dot(c - v1, normal_local) * normal_local;
    }
  }
  const double sep = dot(c - closest, normal_local) - r;
  m.normal = xp.q.apply(normal_local);
  const Vec2 surface_circle = xc.p - r * m.normal;
  const Vec2 surface_poly = xp.apply(closest);
  m.points[0] = {0.5 * (surface_circle + surface_poly), -sep, static_cast<std::uint32_t>(edge)};
  m.count = 1;
  return m;
}

namespace {

// Largest separation of b from a along a's face normals, and that face.
double max_separation(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb, std::size_t& face) {
  double best = -std::numeric_limits<double>::infinity();
  face = 0;
  for (std::size_t i = 0; i < a.vertices.size(); ++i) {
    // Face normal and vertex of a, expressed in b's frame.
    const Vec2 n = xb.q.apply_inverse(xa.q.apply(a.normals[i]));
    const Vec2 v = xb.apply_inverse(xa.apply(a.vertices[i]));
    double deepest = std::numeric_limits<double>::infinity();
    for (const Vec2& w : b.vertices) deepest = std::min(deepest, dot(n, w - v));
    if (deepest > best) {
      best = deepest;
      face = i;
    }
  }
  return best;
}

struct ClipVertex {
  Vec2 v;
  std::uint32_t id;
};

// Keeps the part of segment [in0, in1] on the inner side of the plane
// dot(normal, x) = offset.
int clip_segment(ClipVertex out[2], const ClipVertex in[2], Vec2 normal, double offset, std::uint32_t clip_id) {
  int n = 0;
  const double d0 = dot(normal, in[0].v) - offset;
  const double d1 = dot(normal, in[1].v) - offset;
  if (d0 <= 0) out[n++] = in[0];
  if (d1 <= 0) out[n++] = in[1];
  if (d0 * d1 < 0) {
    const double t = d0 / (d0 - d1);
    out[n].v = in[0].v + t * (in[1].v - in[0].v);
    out[n].id = clip_id;
    ++n;
  }
  return n;
}

}  // namespace

Manifold collide_polygons(const Shape& a, const Transform& xa, const Shape& b, const Transform& xb) {
  Manifold m;
  std::size_t face_a, face_b;
  const double sep_a = max_separation(a, xa, b, xb, face_a);
  if (sep_a > 0) return m;
  const double sep_b = max_separation(b, xb, a, xa, face_b);
  if (sep_b > 0) return m;

  // Reference face: the one with the larger separation, preferring A when
  // they are nearly equal so the choice doesn't flicker between steps.
  const Shape* ref = &a;
  const Shape* inc = &b;
  const Transform* xr = &xa;
  const Transform* xi = &xb;
  std::size_t ref_face = face_a;
  bool flip = false;
  if (sep_b > sep_a + 0.1 * kLinearSlop) {
    std::swap(ref, inc);
    std::swap(xr, xi);
    ref_face = face_b;
    flip = true;
  }

  // Incident edge: the face of `inc` most anti-parallel to the reference
  // normal.
  const Vec2 ref_normal = xr->q.apply(ref->normals[ref_face]);
  std::size_t inc_edge = 0;
  double min_dot = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < inc->normals.size(); ++i) {
    const double d = dot(ref_normal, xi->q.apply(inc->normals[i]));
    if (d < min_dot) {
      min_dot = d;
      inc_edge = i;
    }
  }
  const std::size_t ni = inc->vertices.size();
  ClipVertex incident[2] = {
      {xi->apply(inc->vertices[inc_edge]), static_cast<std::uint32_t>(inc_edge << 8)},
      {xi->apply(inc->vertices[(inc_edge + 1) % ni]), static_cast<std::uint32_t>(((inc_edge + 1) % ni) << 8)},
  };

  const std::size_t nr = ref->vertices.size();
  const Vec2 r1 = xr->apply(ref->vertices[ref_face]);
  const Vec2 r2 = xr->apply(ref->vertices[(ref_face + 1) % nr]);
  const Vec2 tangent = (r2 - r1).normalized();

  // Clip the incident edge to the side planes of the reference face.
  ClipVertex clip1[2], clip2[2];
  if (clip_segment(clip1, incident, -tangent, -dot(tangent, r1), static_cast<std::uint32_t>(ref_face) | 0x10000u) < 2)
    return m;
  if (clip_segment(clip2, clip1, tangent, dot(tangent, r2), static_cast<std::uint32_t>(ref_face) | 0x20000u) < 2)
    return m;

  m.normal = flip ? -ref_normal : ref_normal;
  const double front = dot(ref_normal, r1);
  for (const ClipVertex& cv : clip2) {
    const double s = dot(ref_normal, cv.v) - front;
    if (s > 0) continue;  // above the reference face: not touching
    ContactPoint& p = m.points[m.count++];
    // Midpoint between the incident point and its projection on the face.
    p.point = cv.v - 0.5 * s * ref_normal;
    p.penetration = -s;
    p.feature = cv.id | (flip ? 0x80000000u : 0u) | static_cast<std::uint32_t>(ref_face) << 16;
  }
  return m;
}

}  // namespace lattice
