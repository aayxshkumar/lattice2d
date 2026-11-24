#pragma once

#include <cstdint>
#include <vector>

#include "lattice/body.hpp"
#include "lattice/collision.hpp"

namespace lattice {

struct WorldConfig {
  Vec2 gravity{0, -9.81};
  int velocity_iterations = 10;
  // Baumgarte stabilisation: the fraction of penetration (beyond the slop)
  // removed per step through the velocity solver.
  double baumgarte = 0.2;
  double linear_slop = 0.005;  // penetration allowed without correction (m)
  // Closing speeds below this don't bounce, so resting contacts don't
  // jitter.
  double restitution_threshold = 1.0;  // m/s
};

struct StepStats {
  std::size_t broadphase_pairs = 0;
  std::size_t contacts = 0;       // manifolds with at least one point
  std::size_t contact_points = 0;
  double max_penetration = 0;
};

// A contact being solved this step.
struct Contact {
  std::uint32_t a, b;  // body indices
  Manifold manifold;
  double friction, restitution;
  struct Point {
    Vec2 ra, rb;  // from each centre of mass to the contact point
    double normal_mass, tangent_mass;
    double normal_impulse = 0, tangent_impulse = 0;  // accumulated this step
    double velocity_bias = 0;                        // restitution target
    double position_bias = 0;                        // Baumgarte target
  } points[2];
};

// Fixed-step 2D rigid-body world.
//
// Each step: find overlapping pairs (sort and sweep on x), build contact
// manifolds, apply gravity, solve contacts with sequential impulses, then
// integrate positions (semi-implicit Euler).
//
// Accumulated impulses start from zero every step: there is no warm
// starting yet, which is why tall stacks need many iterations and still
// creep (see README "Known gaps").
class World {
 public:
  explicit World(WorldConfig cfg = {}) : cfg_(cfg) {}

  std::uint32_t add(const BodyDef& def);
  void step(double dt);

  std::vector<Body>& bodies() { return bodies_; }
  const std::vector<Body>& bodies() const { return bodies_; }
  const std::vector<Contact>& contacts() const { return contacts_; }
  const StepStats& last_step() const { return stats_; }
  const WorldConfig& config() const { return cfg_; }
  std::uint64_t steps() const { return steps_; }

  // Sum over dynamic bodies.
  double kinetic_energy() const;
  Vec2 linear_momentum() const;
  // FNV-1a over every body's position, angle and velocities, bit for bit.
  std::uint64_t state_hash() const;

 private:
  void find_contacts();
  void prepare_contacts(double dt);
  void solve_velocities();

  WorldConfig cfg_;
  std::vector<Body> bodies_;
  std::vector<Contact> contacts_;
  StepStats stats_;
  std::uint64_t steps_ = 0;
};

}  // namespace lattice
