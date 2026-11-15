# Flow-rate partition objective: verification (measured 2026-11-15, FoamWin)

`objective { type flowRatePartition; }` (objective.H, objectiveUpdate.H, `objectiveOutletPressure`,
`objectiveOutletVelocity`). All runs serial, `bin/adjointFoam.exe` built with `build.sh` against the
FoamWin tree, hand-written cases (blockMesh + snappyHexMesh, no Stokes code). Case generator,
runners and analysis scripts were kept outside this repository.

## V0: regression (no objective block)

Laminar channel 60 x 20 x 20 mm, sphere r = 4 mm (snappy, level 1), 0.1 m/s, nu = 1.5e-5,
pressure-loss adjoint with the legacy outlet conditions; converged in 141 iterations.

| Binary | `objective.dat` | `U`, `p`, `phi`, `Ua`, `pa`, `phia`, `shapeSensitivity` at 141 | log |
|---|---|---|---|
| before this change vs after, no `objective` block | byte-identical | byte-identical | identical except the clock lines |
| before vs after, with `etc/adjointDict` (`objective { type ptLosses; }`) | byte-identical | byte-identical | |

## Fail-fast checks

Each refused at start-up with a FatalIOError naming the problem: one outlet; `targetFractions`
not summing to 1; an outlet `pa` of type `adjointOutletPressure`; an inlet with Ua = -U_in;
an unknown objective type.

## The case (V1, V2)

T-junction, 10 x 10 mm square section: inlet leg 40 mm (y from -40 to 0 mm, U = 0.2 m/s, Re about
130), cross channel x in [-45, 45] mm, outlets `out1` (x = -45 mm) and `out2` (x = +45 mm), fixed
p = 0 and `pressureInletOutletVelocity`. A spherical bump on the top wall (y = 10 mm, z = 5 mm),
centred 20 mm from the junction in each branch, radius R1 (branch 1) and R2 (branch 2), meshed by
snappyHexMesh (castellate + snap, surface level 1, no layers). Laminar air. Adjoint: `UaGradU`,
`div(-phi,Ua) bounded Gauss upwind`, pa relaxation 0.3, `viscousTerm true`, `tangential robin`,
residual control 1e-6 on U, p, Ua, pa.

R1 = 3 mm, R2 = 2 mm gave a split of 0.488 / 0.512 (inside 0.48 to 0.52), so R1 was raised to
4 mm as the design says: split 0.466 / 0.534 at 1 mm, 0.465 / 0.535 at 0.5 mm; branch 2 over-fed.

dJ/dR_i (adjoint) = -integral over bump i of `shapeSensitivity` dA (a larger sphere moves its
surface into the fluid, against the normal that points out of the fluid).

## V1: invariants (every run of V2, 16 runs plus variants)

| Check | Target | Measured (worst run) |
|---|---|---|
| J logged by the solver vs J from the written `phi` (via Stokes' `port_flows`) | 1e-6 relative | 2.2e-12 |
| Euler sum, \|sum c_i Q_i\| / (max\|c_i\| Q) | 1e-10 | 1.1e-12 |
| adjoint port flux \|sum phia_out\| / sum \|phia_out\| | 1e-3 | 6.2e-9 (1 mm), 9.1e-10 (0.5 mm) |
| U, p, Ua, pa reach the residual target | all | all: 172 to 179 iterations at 1 mm, 334 to 336 at 0.5 mm |
| max\|Ua\| steady over the last 50 iterations | within 1 % | 1 mm baseline run on to 600 iterations: max\|Ua\| = 1.323371e4 and both dJ/dR fixed to 7 digits from iteration 225 on; 0.0004 % between 125 and 175 |
| backflow at an outlet | none | none |

No convergence trouble: none of the design's cuts (viscous term off, tangential zero, pa
relaxation, startIteration, ATC limiter) was needed. `startIteration` is implemented (default 0)
but was not exercised.

## V2: adjoint vs central finite differences

Noise floor: the baseline solved twice gives the same J to every digit (serial snappy is
deterministic); the floor is taken as the residual of a 0.02 mm step about the central slope.
Conclusive = forward and backward within 30 % of central and |J+ - J-| >= 5 floors. The
0.25 mm step was `curvature`-inconclusive for R1 at both meshes (forward and backward 33 % apart),
so the step was halved to 0.125 mm as ADR-0040's rule says; the table is at 0.125 mm.

| Mesh (cells) | J | | central | fwd / bwd | signal / floor | adjoint (solver `shapeSensitivity`) | ratio | Stokes agent reconstruction | ratio |
|---|---|---|---|---|---|---|---|---|---|
| 1 mm (15,349) | 2.3585e-3 | dJ/dR1 | 4.057 | 4.370 / 3.745 | 202 | 4.973 | **1.226** | 6.003 | 1.480 |
| | | dJ/dR2 | -1.019 | -1.090 / -0.947 | 307 | -1.157 | **1.136** | -1.424 | 1.399 |
| 0.5 mm (110,944) | 2.4741e-3 | dJ/dR1 | 4.318 | 4.675 / 3.961 | 8,203 | 5.010 | **1.160** | 5.773 | 1.337 |
| | | dJ/dR2 | -1.070 | -1.112 / -1.027 | 415 | -1.328 | **1.241** | -1.597 | 1.493 |

(units 1/m; J dimensionless.) Cosine of the two-component gradients: 0.9998 (1 mm), 0.9999 (0.5 mm).
At the 0.25 mm step the central differences were 4.063 / -1.033 (1 mm) and 4.329 / -1.079 (0.5 mm),
within 1 % of the 0.125 mm ones.

Against the design's pass list:

| Criterion | Result |
|---|---|
| signs: dJ/dR1 > 0, dJ/dR2 < 0 | **pass** (both meshes, adjoint and FD) |
| check conclusive | **pass** at 0.125 mm (inconclusive, `curvature`, at 0.25 mm for R1) |
| 0.5 mm: both ratios in [0.8, 1.25], cosine >= 0.95 | **pass**, narrowly for R2 (1.241) |
| \|ratio - 1\| at 0.5 mm not larger than at 1 mm | **R1 pass** (0.226 -> 0.160), **R2 fail** (0.136 -> 0.241) |
| aim [0.9, 1.1] at 0.5 mm | not met |

So V2 does **not** pass as written: the adjoint is high by 14 to 24 % on both meshes with the right
direction (cosine 0.9999), and R2's error grows with refinement.

### Which terms matter (baseline only; the primal and so the FD reference are unchanged)

| Variant | 1 mm R1 / R2 ratio | 0.5 mm R1 / R2 ratio |
|---|---|---|
| as above | 1.226 / 1.136 | 1.160 / 1.241 |
| `viscousTerm false` | 1.226 / 1.136 | 1.159 / 1.240 |
| `tangential zero` | 1.223 / 1.134 | 1.156 / 1.237 |
| U outlets `inletOutlet` (outlet cells no longer counted as inlets, so the ATC is not zeroed there) | 1.209 / 1.119 | 1.149 / 1.229 |
| `div(-phi,Ua) bounded Gauss linearUpwind grad(Ua)` | 1.166 / 1.111 | 1.057 / 1.143 |
| both of the last two | 1.150 / 1.096 | 1.047 / 1.123 |

- The outlet conditions' optional terms are not where the error is (below 0.4 %).
- **ATC at the outlets**: `adjointControls.H` counts every patch whose U `fixesValue()` and is not
  `mixed` as an inlet; `pressureInletOutletVelocity` is `directionMixed`, so the outlets of a
  Stokes case are counted too and the ATC is zeroed in their cells (about 1 % here). Left as is
  because it is the pressure-loss path's behaviour (V0); `objective.H` uses a stricter rule for its
  own inlet check.
- **Upwind adjoint convection** is the largest single term (6 to 10 %). With `linearUpwind` the
  R1 ratio converges (1.166 -> 1.057); R2 still grows slightly (1.111 -> 1.143).
- What remains (about 5 % for R1 and 12 to 14 % for R2 at 0.5 mm with both changes) is most
  likely the first-order wall gradient on snapped cells (the design's first risk; ADR-0038 saw 1.24
  at 1 mm for the pressure loss). Not proved: the R2 trend with refinement is not what a pure
  first-order wall error would do, and the bump-to-wall junction (snapped at a 90 degree corner,
  no layers) was not studied.

### Stokes' agent reconstruction is worse on snapped meshes

`stokes_agent.sensitivity.wall_sensitivity` (the formula `job.sensitivity` uses) is 15 to 20 %
above the solver's own `shapeSensitivity` here (ratios 1.34 to 1.49 against FD), though they agreed
within 0.4 % on the step-1 obstacle case. It takes the owner cell's centre as the mean of its face
centres; OpenFOAM uses the volume centroid, and on snapped polyhedra the two differ. Reading the
solver's written `shapeSensitivity` instead of reconstructing it is the cheaper fix.

## Not done here

- V3 (the manifold through Stokes) and `tools/tjunction_demo.py`: the Stokes half.
- A finite-difference check with boundary layers or a level-2 bump surface, to separate the wall
  gradient from the R2 trend.
- `startIteration`, parallel runs, turbulent cases.
