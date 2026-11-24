# lattice2d

A 2D rigid-body physics engine in C++20: convex polygons and circles, SAT collision with two-point contact manifolds, a sequential-impulse solver with friction and restitution, a sort-and-sweep broadphase, and deterministic fixed-step simulation. Every claim is a measured check against physics or against a known answer, and the checks the engine can't pass yet are listed openly.

```
$ ./build/lattice_sim --check
check                result            value    threshold  measures
freefall             PASS        1.42109e-14        1e-09  |y - closed form| after 2 s (m)
drop-rests           PASS        3.50474e-05         0.01  speed + |spin| after 3 s
drop-penetration     PASS         0.00500095         0.02  resting penetration (m)
drop-level           KNOWN GAP    2.89324e-05        1e-06  |angle| of a box dropped flat, after 3 s (rad)
billiards-momentum   PASS        2.71889e-16        1e-09  relative momentum change over 10 s
billiards-energy     PASS        1.88332e-16         0.02  relative kinetic energy change over 10 s (elastic)
ramp-holds           PASS         0.00279286         0.05  box displacement on a 20 deg slope, mu 0.5, 3 s (m)
ramp-slides          PASS         0.00877947         0.05  relative error vs g(sin - mu cos) t^2/2 on 35 deg, 2 s
stack-3              PASS          0.0185216          0.1  top box sideways drift after 10 s, 3 boxes (m)
stack-10             KNOWN GAP        11.9265          0.1  top box sideways drift after 10 s, 10 boxes (m)
stack-20             KNOWN GAP        8.13493          0.1  top box sideways drift after 10 s, 20 boxes (m)
pyramid-20           KNOWN GAP         3.1207          0.1  top box sideways drift after 10 s, 20 rows (m)
bullet               KNOWN GAP            300            5  bullet x after 1 s (wall at x = 5; must stay below 5)
determinism          PASS                  0          0.5  rain(200) run twice: hashes differ (0 = identical)

0 failed, 5 known gaps
```

No dependencies beyond googletest. Built on Apple clang 17 and GCC 13.

## Build and run

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/lattice_sim --check                       # the physical checks above
./build/lattice_sim --scene pyramid --n 12        # run a scene, one summary line per second
./build/lattice_sim --scene rain --trace t.json   # write a trace, then open tools/viewer.html and load it
./build/lattice_sim --bench --n 1000              # time the rain scene
ctest --test-dir build                            # unit tests + the checks
python3 tools/mutants.py                          # plants physics bugs, reports what catches each
```

Scenes: `drop stack pyramid ramp billiards bullet rain`. `--n` sizes them (stack height, pyramid rows, ramp angle in degrees, body count). googletest 1.15.2 comes from `third_party/` when present (offline machines), else GitHub, pinned to one SHA-256. Sanitizers: `-DLATTICE_SANITIZERS=address,undefined`.

`tools/viewer.html` is a single offline page: load a trace and it replays it on a canvas, with play, pause and scrubbing.

## How a step works

```
find_contacts   AABBs → sort and sweep on x → overlapping pairs (sorted, deterministic)
                → narrow phase → contact manifolds (up to 2 points each)
integrate v     v += dt · (g + F/m),  ω += dt · τ/I
prepare         per contact point: effective normal and tangent mass, restitution target
                (only for impacts above 1 m/s), Baumgarte target for penetration beyond 5 mm
solve × 10      sequential impulses: normal impulse clamped ≥ 0 (accumulated, not per
                increment), then friction clamped to μ · normal impulse
integrate x     x += dt · v,  θ += dt · ω
```

- **Narrow phase** (`src/collision.cpp`): circle–circle; polygon–circle with face and vertex regions; polygon–polygon by the separating axis test, then clipping the incident edge against the reference face's side planes for a two-point manifold. Reference-face selection is biased toward A, so it doesn't flicker between steps. Each contact point gets a feature id, but nothing uses it yet.
- **Shapes** (`src/shape.cpp`): polygons are built from any point set (convex hull, collinear points removed) and recentred on the centroid. Mass and inertia come from a triangle fan with the parallel-axis shift.
- **Determinism:** one thread, sorted pair order, no FMA contraction (`-ffp-contract=off`). The same scene on the same build gives a bit-identical state hash. It is *not* bit-identical across platforms yet, because `sin` and `cos` come from the system libm.

## Known gaps

These are measured by `--check` and reported as KNOWN GAP. They don't fail the build:

1. **No warm starting.** Accumulated impulses restart from zero every step, so the solver never converges on tall stacks. Measured by sweeping iterations: a 10-box tower falls over at 10 iterations, wobbles at 30 and stands at 100; a 20-box tower needs about 300. The contact feature ids are already computed for matching contacts across steps.
2. **Sequential two-point solve.** Even one box dropped flat comes to rest tilted (2.9e-5 rad), because its two corner contacts are solved one after the other. At 50 iterations the tilt is 1e-18, so the cause is the order of solving, not the geometry. A block solver (solving both points together) or warm starting fixes it.
3. **No continuous collision.** A 0.1 m ball at 300 m/s moves 5 m per step and passes straight through a 0.1 m wall.
4. **Baumgarte stabilisation** turns penetration correction into velocity, which adds energy and makes stacks bouncy. Split impulses or position-level correction (NGS) avoid that.

## What catches what

`tools/mutants.py` plants each bug in a scratch copy and runs the unit tests and `--check`:

| Mutant | Bug | Unit tests | Checks |
|---|---|---|---|
| friction-sign | friction pushes along the slide | 3 fail | 4 fail |
| no-baumgarte | penetration never corrected | 1 fail | 2 fail |
| gravity-after-solve | gravity applied twice | 2 fail | 3 fail |
| broadphase-neighbours-only | sweep pairs only adjacent bodies | 1 fail | 1 fail |
| inertia-sign | parallel-axis shift backwards | 2 fail | pass |
| circle-normal-unflipped | circle–polygon normal wrong when the circle is A | 1 fail | pass |
| reference-face-flip-lost | polygon normal wrong when B holds the reference face | 1 fail | pass |
| clamp-increment | each impulse increment clamped instead of the total | pass | pass |
| restitution-min | restitution mixing takes the smaller value | pass | pass |

Two survive. No scene mixes bodies with different restitution, and clamp-increment hurts convergence, which only shows on the stack checks. Those are known gaps, so they can't fail. Three more are caught only by unit tests: the rain scene, which mixes all shape pairs, is checked for determinism but not for penetration or energy.

## Performance (Apple M4, Release)

`lattice_sim --bench --n 1000`: 1,003 bodies piling into a container, 600 steps, **0.61 ms per step** (about 1,840 broadphase pairs and 1,770 contact points per step).

## Roadmap

1. **Warm starting** with contact matching by feature id, so 10 iterations hold the 10- and 20-box towers and the 20-row pyramid. Then a block solver for two-point manifolds.
2. **Continuous collision**: time of impact by conservative advancement for fast bodies, and the bullet check passing.
3. **Joints**: revolute, prismatic, distance and weld, with motors and limits; a ragdoll and a chain bridge scene.
4. **Sleeping and islands**: bodies at rest stop costing time; islands solved independently.
5. **A faster broadphase** (dynamic AABB tree) and **multithreaded islands**, with determinism kept.
6. **Cross-platform determinism**: deterministic sin/cos, so a seed gives the same hash on every machine (what lockstep multiplayer needs).
7. **3D**: GJK/EPA, boxes, spheres and capsules.
8. **Soft bodies and fluids** (position-based dynamics, SPH) in the same world.
9. **WebAssembly build** of the viewer that simulates live in the browser.

## Layout

```
include/lattice/   math, shape, body, collision, world, scenes
src/               shape, collision, world, scenes (scenes + physical checks)
tools/             lattice_sim.cpp (CLI), viewer.html, mutants.py
tests/             shape, collision and world unit tests
```

MIT licence.
