#!/usr/bin/env python3
"""Plant known physics bugs one at a time and see what notices.

Each mutant is applied to a scratch copy of the repo (the working tree is
never touched) and built in its own build dir. Then the unit tests and the
physical checks (lattice_sim --check) run against it.

  python3 tools/mutants.py                 all mutants
  python3 tools/mutants.py friction-sign   one mutant
  python3 tools/mutants.py --list

Exit code 0 if every mutant was caught, 1 otherwise.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MUTANTS = {
    "friction-sign": (
        "src/world.cpp",
        "      double dpt = -p.tangent_mass * vt;",
        "      double dpt = p.tangent_mass * vt;",
        "friction pushes along the sliding direction instead of against it",
    ),
    "clamp-increment": (
        "src/world.cpp",
        "      p.normal_impulse = std::max(old_n + dpn, 0.0);\n      dpn = p.normal_impulse - old_n;",
        "      dpn = std::max(dpn, 0.0);\n      p.normal_impulse = old_n + dpn;",
        "each normal impulse increment is clamped, so later iterations can't take impulse back",
    ),
    "no-baumgarte": (
        "src/world.cpp",
        "      p.position_bias = cfg_.baumgarte / dt * std::max(0.0, cp.penetration - cfg_.linear_slop);",
        "      p.position_bias = 0 * dt;",
        "penetration is never corrected, so resting bodies sink",
    ),
    "gravity-after-solve": (
        "src/world.cpp",
        "  prepare_contacts(dt);\n  for (int i = 0; i < cfg_.velocity_iterations; ++i) solve_velocities();\n",
        "  prepare_contacts(dt);\n  for (int i = 0; i < cfg_.velocity_iterations; ++i) solve_velocities();\n"
        "  for (Body& b : bodies_)\n    if (!b.is_static()) b.velocity += dt * cfg_.gravity;\n",
        "gravity is applied twice, once after the contacts were solved",
    ),
    "inertia-sign": (
        "src/shape.cpp",
        "  m.inertia = density * i_origin - m.mass * dot(c, c);",
        "  m.inertia = density * i_origin + m.mass * dot(c, c);",
        "polygon inertia shifted the wrong way by the parallel axis theorem",
    ),
    "circle-normal-unflipped": (
        "src/collision.cpp",
        "  if (a.kind == K::Circle && b.kind == K::Polygon) return flipped(collide_polygon_circle(b, xb, a, xa));",
        "  if (a.kind == K::Circle && b.kind == K::Polygon) return flipped(flipped(collide_polygon_circle(b, xb, a, xa)));",
        "circle-vs-polygon normal points the wrong way when the circle is body A",
    ),
    "reference-face-flip-lost": (
        "src/collision.cpp",
        "  m.normal = flip ? -ref_normal : ref_normal;",
        "  m.normal = ref_normal;",
        "polygon normal not flipped when B's face is the reference face",
    ),
    "broadphase-neighbours-only": (
        "src/world.cpp",
        "    for (std::size_t j = i + 1; j < entries.size() && entries[j].box.lo.x <= entries[i].box.hi.x; ++j) {",
        "    for (std::size_t j = i + 1; j < entries.size() && j < i + 2 && entries[j].box.lo.x <= entries[i].box.hi.x; ++j) {",
        "the sweep only pairs each body with its next neighbour along x",
    ),
    "restitution-min": (
        "src/world.cpp",
        "    c.restitution = std::max(bodies_[a].restitution, bodies_[b].restitution);",
        "    c.restitution = std::min(bodies_[a].restitution, bodies_[b].restitution);",
        "contact restitution takes the smaller of the two bodies' values",
    ),
}


def run(cmd, cwd, timeout=900):
    try:
        return subprocess.run(cmd, cwd=cwd, text=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        return subprocess.CompletedProcess(cmd, -999, e.stdout or "", "timeout")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    args = ap.parse_args()
    if args.list:
        for n, m in MUTANTS.items():
            print(f"{n:28} {m[3]}")
        return 0
    names = args.names or list(MUTANTS)
    for n in names:
        if n not in MUTANTS:
            print(f"unknown mutant {n}", file=sys.stderr)
            return 2

    ignore = shutil.ignore_patterns("build*", ".git", "__pycache__")
    survivors = []
    print(f"{'mutant':28} {'unit tests':>11} {'checks':>10}")
    with tempfile.TemporaryDirectory(prefix="lattice2d-mutants-") as tmp:
        for name in names:
            path, old, new, what = MUTANTS[name]
            src = os.path.join(tmp, name)
            shutil.copytree(ROOT, src, ignore=ignore)
            target = os.path.join(src, path)
            text = open(target).read()
            if text.count(old) != 1:
                print(f"{name}: anchor not found exactly once in {path}; update tools/mutants.py")
                return 2
            open(target, "w").write(text.replace(old, new))
            build = os.path.join(src, "build")
            c = run(["cmake", "-S", src, "-B", build, "-DCMAKE_BUILD_TYPE=Release"], src)
            b = run(["cmake", "--build", build, "-j", str(args.jobs)], src)
            if c.returncode or b.returncode:
                print(f"{name}: does not compile\n{b.stdout[-2500:]}{b.stderr[-1500:]}")
                return 2
            t = run([os.path.join(build, "lattice_tests"), "--gtest_brief=1"], src)
            failed_tests = sorted(set(re.findall(r"\[  FAILED  \] (\S+) \(", t.stdout)))
            unit = "pass" if t.returncode == 0 else (f"{len(failed_tests)} fail" if failed_tests else "crash")
            k = run([os.path.join(build, "lattice_sim"), "--check"], src)
            failed_checks = [l.split()[0] for l in k.stdout.splitlines() if re.match(r"^\S+\s+FAIL\b", l)]
            checks = "pass" if k.returncode == 0 else (f"{len(failed_checks)} fail" if failed_checks else "crash")
            caught = t.returncode != 0 or k.returncode != 0
            print(f"{name:28} {unit:>11} {checks:>10}   {'' if caught else 'SURVIVED  '}{what}")
            if failed_checks:
                print(f"{'':28} checks: {', '.join(failed_checks)}")
            if failed_tests:
                print(f"{'':28} tests: {', '.join(failed_tests[:4])}{' ...' if len(failed_tests) > 4 else ''}")
            if not caught:
                survivors.append(name)
    if survivors:
        print(f"\n{len(survivors)} mutant(s) survived: {', '.join(survivors)}")
        return 1
    print("\nall mutants caught")
    return 0


if __name__ == "__main__":
    sys.exit(main())
