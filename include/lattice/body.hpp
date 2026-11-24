#pragma once

#include <cstdint>

#include "lattice/shape.hpp"

namespace lattice {

struct BodyDef {
  Shape shape;
  Vec2 position;
  double angle = 0;
  Vec2 velocity;
  double angular_velocity = 0;
  double density = 1;       // 0 makes the body static
  double friction = 0.5;
  double restitution = 0;
};

struct Body {
  std::uint32_t id = 0;
  Shape shape;
  Vec2 position;  // centre of mass
  double angle = 0;
  Vec2 velocity;
  double angular_velocity = 0;
  Vec2 force;
  double torque = 0;
  double mass = 0, inv_mass = 0;
  double inertia = 0, inv_inertia = 0;
  double friction = 0.5;
  double restitution = 0;

  bool is_static() const { return inv_mass == 0; }
  Transform transform() const { return {position, Rot::from_angle(angle)}; }
  double kinetic_energy() const {
    return 0.5 * mass * velocity.length_squared() + 0.5 * inertia * angular_velocity * angular_velocity;
  }
};

}  // namespace lattice
