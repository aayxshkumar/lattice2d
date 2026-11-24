// lattice_sim: run scenes, physical checks and a benchmark.
//
//   lattice_sim --check                     every physical check; exit 1 on a failure
//   lattice_sim --scene pyramid --n 15      run a scene, print a summary each second
//   lattice_sim --scene rain --trace t.json write a trace for tools/viewer.html
//   lattice_sim --bench [--n 1000]          time the rain scene
//
// Checks marked KNOWN GAP are reported but don't fail --check.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#include "lattice/scenes.hpp"

using namespace lattice;

namespace {

int usage() {
  std::fprintf(stderr,
               "usage: lattice_sim --check\n"
               "       lattice_sim --scene NAME [--n N] [--steps S] [--seed K] [--trace FILE] [--every K]\n"
               "       lattice_sim --bench [--n N] [--steps S]\n"
               "scenes:");
  for (const auto& s : scene_names()) std::fprintf(stderr, " %s", s.c_str());
  std::fprintf(stderr, "\n");
  return 2;
}

int do_check() {
  int failed = 0, gaps = 0, gaps_passing = 0;
  std::printf("%-20s %-8s %14s %12s  %s\n", "check", "result", "value", "threshold", "measures");
  for (const CheckResult& c : run_checks()) {
    const char* status = c.pass ? (c.known_gap ? "PASS*" : "PASS") : (c.known_gap ? "KNOWN GAP" : "FAIL");
    std::printf("%-20s %-8s %14.6g %12.6g  %s\n", c.name.c_str(), status, c.value, c.threshold, c.what.c_str());
    if (!c.pass && !c.known_gap) ++failed;
    if (c.known_gap) ++gaps;
    if (c.known_gap && c.pass) ++gaps_passing;
  }
  std::printf("\n%d failed, %d known gaps", failed, gaps - gaps_passing);
  if (gaps_passing) std::printf(" (%d known-gap check(s) now pass: PASS* — update the README)", gaps_passing);
  std::printf("\n");
  return failed ? 1 : 0;
}

void write_trace(std::FILE* f, const World& w, bool first_frame, bool header) {
  if (header) {
    std::fprintf(f, "{\"bodies\":[");
    for (std::size_t i = 0; i < w.bodies().size(); ++i) {
      const Body& b = w.bodies()[i];
      std::fprintf(f, "%s{\"static\":%s,", i ? "," : "", b.is_static() ? "true" : "false");
      if (b.shape.kind == Shape::Kind::Circle) {
        std::fprintf(f, "\"circle\":%.6g}", b.shape.radius);
      } else {
        std::fprintf(f, "\"poly\":[");
        for (std::size_t k = 0; k < b.shape.vertices.size(); ++k)
          std::fprintf(f, "%s[%.6g,%.6g]", k ? "," : "", b.shape.vertices[k].x, b.shape.vertices[k].y);
        std::fprintf(f, "]}");
      }
    }
    std::fprintf(f, "],\"frames\":[");
    return;
  }
  std::fprintf(f, "%s[", first_frame ? "" : ",");
  for (std::size_t i = 0; i < w.bodies().size(); ++i) {
    const Body& b = w.bodies()[i];
    std::fprintf(f, "%s%.5g,%.5g,%.5g", i ? "," : "", b.position.x, b.position.y, b.angle);
  }
  std::fprintf(f, "]");
}

int do_scene(const std::string& name, int n, int steps, std::uint64_t seed, const char* trace, int every) {
  Scene s = make_scene(name, n, seed);
  if (steps > 0) s.steps = steps;
  std::FILE* tf = nullptr;
  if (trace) {
    tf = std::fopen(trace, "w");
    if (!tf) {
      std::perror(trace);
      return 1;
    }
    write_trace(tf, s.world, true, true);
    write_trace(tf, s.world, true, false);
  }
  std::printf("scene %s: %zu bodies, dt %.5f, %d steps\n", name.c_str(), s.world.bodies().size(), s.dt, s.steps);
  std::printf("%7s %12s %10s %9s %12s\n", "t (s)", "KE (J)", "contacts", "points", "max pen (m)");
  const int per_second = static_cast<int>(std::lround(1.0 / s.dt));
  for (int i = 1; i <= s.steps; ++i) {
    s.world.step(s.dt);
    if (tf && i % every == 0) write_trace(tf, s.world, false, false);
    if (i % per_second == 0 || i == s.steps) {
      const StepStats& st = s.world.last_step();
      std::printf("%7.2f %12.5g %10zu %9zu %12.5f\n", i * s.dt, s.world.kinetic_energy(), st.contacts,
                  st.contact_points, st.max_penetration);
    }
  }
  if (tf) {
    std::fprintf(tf, "],\"dt\":%.8g}\n", s.dt * every);
    std::fclose(tf);
    std::printf("trace written to %s (open tools/viewer.html and load it)\n", trace);
  }
  std::printf("state hash %016llx\n", static_cast<unsigned long long>(s.world.state_hash()));
  return 0;
}

int do_bench(int n, int steps) {
  Scene s = make_scene("rain", n > 0 ? n : 1000, 1);
  if (steps <= 0) steps = 600;
  const auto t0 = std::chrono::steady_clock::now();
  std::size_t pairs = 0, points = 0;
  for (int i = 0; i < steps; ++i) {
    s.world.step(s.dt);
    pairs += s.world.last_step().broadphase_pairs;
    points += s.world.last_step().contact_points;
  }
  const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("rain: %zu bodies, %d steps, %.3f s, %.3f ms/step, %.0f broadphase pairs/step, %.0f contact points/step\n",
              s.world.bodies().size(), steps, sec, 1000.0 * sec / steps, double(pairs) / steps, double(points) / steps);
  std::printf("state hash %016llx\n", static_cast<unsigned long long>(s.world.state_hash()));
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  std::string scene;
  int n = 0, steps = 0, every = 1;
  std::uint64_t seed = 1;
  const char* trace = nullptr;
  bool check = false, bench = false;
  for (int i = 1; i < argc; ++i) {
    auto value = [&]() -> const char* {
      if (i + 1 >= argc) std::exit(usage());
      return argv[++i];
    };
    if (!std::strcmp(argv[i], "--check")) check = true;
    else if (!std::strcmp(argv[i], "--bench")) bench = true;
    else if (!std::strcmp(argv[i], "--scene")) scene = value();
    else if (!std::strcmp(argv[i], "--n")) n = std::atoi(value());
    else if (!std::strcmp(argv[i], "--steps")) steps = std::atoi(value());
    else if (!std::strcmp(argv[i], "--seed")) seed = std::strtoull(value(), nullptr, 10);
    else if (!std::strcmp(argv[i], "--trace")) trace = value();
    else if (!std::strcmp(argv[i], "--every")) every = std::max(1, std::atoi(value()));
    else return usage();
  }
  try {
    if (check) return do_check();
    if (bench) return do_bench(n, steps);
    if (!scene.empty()) return do_scene(scene, n, steps, seed, trace, every);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return usage();
}
