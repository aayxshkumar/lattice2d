#pragma once

#include <cmath>
#include <cstdint>

namespace lattice {

struct Vec2 {
  double x = 0, y = 0;

  constexpr Vec2() = default;
  constexpr Vec2(double x_, double y_) : x(x_), y(y_) {}

  constexpr Vec2 operator+(Vec2 o) const { return {x + o.x, y + o.y}; }
  constexpr Vec2 operator-(Vec2 o) const { return {x - o.x, y - o.y}; }
  constexpr Vec2 operator-() const { return {-x, -y}; }
  constexpr Vec2 operator*(double s) const { return {x * s, y * s}; }
  constexpr Vec2& operator+=(Vec2 o) { x += o.x; y += o.y; return *this; }
  constexpr Vec2& operator-=(Vec2 o) { x -= o.x; y -= o.y; return *this; }
  constexpr bool operator==(const Vec2&) const = default;

  double length() const { return std::sqrt(x * x + y * y); }
  constexpr double length_squared() const { return x * x + y * y; }
  Vec2 normalized() const {
    const double l = length();
    return l > 0 ? Vec2{x / l, y / l} : Vec2{};
  }
};

constexpr Vec2 operator*(double s, Vec2 v) { return {s * v.x, s * v.y}; }
constexpr double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
// 2D cross product: the z component of the 3D cross product.
constexpr double cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }
// cross(v, s) = v x (s * z): rotates v by -90 degrees and scales it.
constexpr Vec2 cross(Vec2 v, double s) { return {s * v.y, -s * v.x}; }
// cross(s, v) = (s * z) x v: angular velocity times a lever arm.
constexpr Vec2 cross(double s, Vec2 v) { return {-s * v.y, s * v.x}; }

// A rotation stored as cosine and sine.
struct Rot {
  double c = 1, s = 0;

  static Rot from_angle(double a) { return {std::cos(a), std::sin(a)}; }
  constexpr Vec2 apply(Vec2 v) const { return {c * v.x - s * v.y, s * v.x + c * v.y}; }
  constexpr Vec2 apply_inverse(Vec2 v) const { return {c * v.x + s * v.y, -s * v.x + c * v.y}; }
};

struct Transform {
  Vec2 p;
  Rot q;
  constexpr Vec2 apply(Vec2 v) const { return q.apply(v) + p; }
  constexpr Vec2 apply_inverse(Vec2 v) const { return q.apply_inverse(v - p); }
};

struct AABB {
  Vec2 lo, hi;
  constexpr bool overlaps(const AABB& o) const {
    return lo.x <= o.hi.x && o.lo.x <= hi.x && lo.y <= o.hi.y && o.lo.y <= hi.y;
  }
};

constexpr double kPi = 3.14159265358979323846;

// xoshiro256** seeded through splitmix64, so a scene seed builds the same
// scene with any standard library.
class Rng {
 public:
  explicit Rng(std::uint64_t seed) {
    std::uint64_t x = seed;
    for (auto& w : s_) {
      std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      w = z ^ (z >> 31);
    }
  }
  std::uint64_t next() {
    const std::uint64_t r = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return r;
  }
  // Uniform in [lo, hi).
  double uniform(double lo, double hi) { return lo + (hi - lo) * (static_cast<double>(next() >> 11) * 0x1.0p-53); }

 private:
  static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  std::uint64_t s_[4];
};

}  // namespace lattice
