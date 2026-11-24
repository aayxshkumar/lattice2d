#include "lattice/shape.hpp"

#include <algorithm>

namespace lattice {

Shape Shape::circle(double r) {
  if (!(r > 0)) throw std::invalid_argument("circle radius must be positive");
  Shape s;
  s.kind = Kind::Circle;
  s.radius = r;
  return s;
}

Shape Shape::box(double hw, double hh) {
  if (!(hw > 0 && hh > 0)) throw std::invalid_argument("box half extents must be positive");
  return polygon({{-hw, -hh}, {hw, -hh}, {hw, hh}, {-hw, hh}});
}

Shape Shape::polygon(std::vector<Vec2> pts) {
  // Andrew's monotone chain, giving a counter-clockwise hull without
  // collinear points.
  std::sort(pts.begin(), pts.end(), [](Vec2 a, Vec2 b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
  if (pts.size() < 3) throw std::invalid_argument("polygon needs at least 3 distinct points");
  std::vector<Vec2> hull(2 * pts.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    while (k >= 2 && cross(hull[k - 1] - hull[k - 2], pts[i] - hull[k - 2]) <= 0) --k;
    hull[k++] = pts[i];
  }
  for (std::size_t i = pts.size() - 1, t = k + 1; i-- > 0;) {
    while (k >= t && cross(hull[k - 1] - hull[k - 2], pts[i] - hull[k - 2]) <= 0) --k;
    hull[k++] = pts[i];
  }
  hull.resize(k - 1);
  if (hull.size() < 3) throw std::invalid_argument("polygon points are collinear");
  if (hull.size() > kMaxPolygonVertices) throw std::invalid_argument("polygon has too many vertices");

  Shape s;
  s.kind = Kind::Polygon;
  // Recentre on the centroid.
  const MassProperties mp = [&] {
    Shape tmp;
    tmp.kind = Kind::Polygon;
    tmp.vertices = hull;
    return tmp.mass(1.0);
  }();
  for (Vec2& v : hull) v -= mp.centroid;
  s.vertices = hull;
  for (std::size_t i = 0; i < hull.size(); ++i) {
    const Vec2 e = hull[(i + 1) % hull.size()] - hull[i];
    s.normals.push_back(Vec2{e.y, -e.x}.normalized());
  }
  return s;
}

MassProperties Shape::mass(double density) const {
  MassProperties m;
  if (kind == Kind::Circle) {
    m.mass = density * kPi * radius * radius;
    m.inertia = 0.5 * m.mass * radius * radius;
    return m;
  }
  // Split into triangles from the first vertex; accumulate area, centroid and
  // second moment about the origin, then shift to the centroid.
  double area = 0, i_origin = 0;
  Vec2 c;
  const Vec2 o = vertices[0];
  for (std::size_t i = 1; i + 1 < vertices.size(); ++i) {
    const Vec2 e1 = vertices[i] - o, e2 = vertices[i + 1] - o;
    const double d = cross(e1, e2);
    const double tri = 0.5 * d;
    area += tri;
    c += tri * (1.0 / 3.0) * (e1 + e2);
    const double intx2 = e1.x * e1.x + e2.x * e1.x + e2.x * e2.x;
    const double inty2 = e1.y * e1.y + e2.y * e1.y + e2.y * e2.y;
    i_origin += (0.25 / 3.0 * d) * (intx2 + inty2);
  }
  c = c * (1.0 / area);
  m.mass = density * area;
  m.centroid = o + c;
  // Inertia about the centroid (parallel axis theorem, relative to o).
  m.inertia = density * i_origin - m.mass * dot(c, c);
  return m;
}

AABB Shape::bounds(const Transform& xf) const {
  if (kind == Kind::Circle) return {xf.p - Vec2{radius, radius}, xf.p + Vec2{radius, radius}};
  Vec2 lo = xf.apply(vertices[0]), hi = lo;
  for (const Vec2& v : vertices) {
    const Vec2 w = xf.apply(v);
    lo = {std::min(lo.x, w.x), std::min(lo.y, w.y)};
    hi = {std::max(hi.x, w.x), std::max(hi.y, w.y)};
  }
  return {lo, hi};
}

}  // namespace lattice
