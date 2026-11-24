#pragma once

#include <functional>
#include <string>
#include <vector>

#include "lattice/world.hpp"

namespace lattice {

struct Scene {
  std::string name;
  World world;
  double dt = 1.0 / 60.0;
  int steps = 600;
};

// Named scenes. `n` sizes the scene where that makes sense (stack height,
// pyramid rows, body count); 0 picks the default.
//   drop       one box falling onto the ground
//   stack      a vertical tower of n boxes (default 10)
//   pyramid    a pyramid with n boxes on the bottom row (default 10)
//   ramp       a box on a slope; n is the slope angle in degrees (default 20)
//   billiards  n circles in free space, no gravity, elastic, frictionless (default 12)
//   bullet     a small circle fired at 300 m/s at a thin wall
//   rain       n mixed shapes falling into a container (default 400)
Scene make_scene(const std::string& name, int n = 0, std::uint64_t seed = 1);
std::vector<std::string> scene_names();

// A physical property measured against a threshold.
struct CheckResult {
  std::string name;
  std::string what;      // what is measured, in words
  double value = 0;      // the measurement
  double threshold = 0;  // the bound it must meet
  bool pass = false;
  // A check for something the engine is known not to handle yet. It is
  // reported but doesn't fail the run (see README "Known gaps").
  bool known_gap = false;
};

std::vector<CheckResult> run_checks();

}  // namespace lattice
