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

## Phase B (measured 2026-11-15, later the same day)

Scripts: `stokes-windows/.local/adjoint-flowsplit/tools` (gitignored; `v3_mesh.py`, `v3_sweep.ps1`,
`v3_analyse.py` for V3; `bodyfit.py`, `bf_sweep.ps1`, `bf_compare.py`, `bf_variants.py`,
`rerun.ps1` for the experiments). An adjoint variant is the same mesh solved again with one
setting changed: the primal, so J and the finite differences, is unchanged.

### Solver options added (all opt-in; defaults keep the earlier behaviour)

| `system/adjointDict` key | Values | Default | What |
|---|---|---|---|
| `objective { convection ...; }` | an fvSchemes convection scheme, or `fvSchemes` | `fvSchemes` (the case's `div(-phi,Ua)`) | the objective path's adjoint convection scheme |
| `objective { ATCInlets ...; }` | `legacy`, `strict` | `legacy` | `strict`: the ATC limiter is zeroed only next to the inlets objective.H finds, not in the outlet cells that `adjointControls.H` counts as inlets (hypothesis (c)) |
| `sensitivityWallGradient` (top level) | `uncorrected`, `corrected` | `uncorrected` | `corrected`: the wall gradients in `shapeSensitivity` are taken along the face normal, the cell value moved to the foot of the normal with the cell gradient, (X_w - X_c - k.grad X_c)/(n.d), k = d - n(n.d), d = Cf - Cn (Foundation's `patch.delta()` is already n(n.d), so its boundary snGrad ignores k) |

**V0 with these options in the binary:** the obstacle case without an `objective` block and with
`etc/adjointDict` (`ptLosses`), before (`bin-phaseB0`, the binary of the V2 runs) and after:
`U`, `p`, `phi`, `Ua`, `pa`, `phia`, `shapeSensitivity` at 141 and `objective.dat`
byte-identical, the logs identical except clock lines, and also identical to the run of the
binary from before the objective existed (`v0-old`). With none of the new keys, the final binary
also reproduces the old one byte for byte on the flow-split path (all fields, `shapeSensitivity`
and `objective.dat` of the snapped T-junction at 1 mm and the body-fitted case at 1 mm), so the
case Stokes compiles behaves exactly as in V2.

### V3: the manifold (tools/manifold_demo.py), solver shapeSensitivity, Stokes' mapping

Built and meshed through the Stokes API (`manifold_demo.build`, `parameter.assign` and a rebuild
of every implicit operation and a re-tag for each point, `mesh.generate` with snappy, 1 mm, no
layers: 71,903 cells at (16, 16, 16)); the case is what `stokes.compile.of14.adjoint.write_adjoint`
writes for `adjoint_objective: flow_split` (outlets out1 to out4, targets 0.25, residual 1e-5,
`div(-phi,Ua) bounded Gauss upwind`), run by hand with `adjointFoam` (the program is not beside
the FoamWin install here, so `run.start` would refuse; `shape.check`/`shape.optimise` were not
run). dJ/dh_k = sum over the wall faces of s V_n (n_phi . n) A, with s the solver's written
`shapeSensitivity` and V_n from `stokes.geometry.sensitivity.wall_gradient` (the implicit field
differentiated in h_k at the face centres; it covered 99.8 % of the wall area).

Baseline: fractions 0.2645 / 0.2663 / 0.2419 / 0.2274, J = 1.05253e-3, 325 iterations, 141 s.
Steps of +-0.5 mm in h2, h3, h4, a repeat (identical J) and h2 + 0.05 mm: the tiny step moves J
by +3.5e-6 where the slope predicts -5.8e-6, so the floor is 9.3e-6 (snapped-mesh noise, applied
to every variable).

| Variable | central | fwd / bwd | signal / floor | conclusive | adjoint | ratio |
|---|---|---|---|---|---|---|
| h2 | -0.1162 | -0.1256 / -0.1068 | 12.5 | yes | -0.0873 | **0.752** |
| h3 | -0.1761 | -0.1819 / -0.1702 | 19.0 | yes | -0.1419 | **0.806** |
| h4 | -0.0194 | -0.0010 / -0.0378 | 2.1 | no (fwd and bwd 190 % apart, 2 floors) | -0.0327 | (1.69) |

(1/m.) Cosine over the three: 0.9943; over the two conclusive rows: 0.9995. No wrong sign. The
agent's reconstruction (`wall_sensitivity`, `source="reconstructed"`) gives 0.730 / 0.783.
Adjoint variants on the same baseline:

| Variant | h2 | h3 | h4 (inconclusive) | cosine (3) |
|---|---|---|---|---|
| as Stokes compiles it | 0.752 | 0.806 | 1.69 | 0.9943 |
| `convection` linearUpwind | 0.782 | 0.827 | 1.70 | 0.9949 |
| `sensitivityWallGradient corrected` | 0.772 | 0.826 | 1.74 | 0.9942 |
| both | 0.802 | 0.847 | 1.75 | 0.9948 |
| `ATCInlets strict` | 0.740 | 0.793 | 1.66 | 0.9943 |
| linearUpwind + strict | 0.771 | 0.815 | 1.68 | 0.9949 |

V1 on the nine V3 runs: J from the written phi vs the solver's J within 6.1e-8 (Stokes writes 8
digits); Euler sum from the written phi up to 1.4e-9, limited by those 8 digits (the solver's own
log: 1e-18); adjoint port flux <= 3.8e-8; all converged (313 to 326 iterations); no backflow.
max|Ua| over the last 50 iterations was not recorded.

**Against ADR-0041's V3 pass list, as Stokes compiles the case: fail.** No wrong sign and the
cosine is 0.994, but h2's ratio (0.752) is below 0.8; h3 (0.806) is just inside; h4 is
inconclusive at +-0.5 mm (its J is nearly flat there). Here the adjoint is low, where on the
T-junction it was high. linearUpwind with the corrected wall gradient brings both conclusive rows
inside [0.8, 1.25] (0.802, 0.847), but not by a margin anyone should rely on, and that pair makes
the T-junction worse (below).

### R2 on the T-junction: the three hypotheses

Snapped T-junction of V2 (central differences at 0.125 mm from the V2 sweeps):

| Variant | 1 mm R1 / R2 | 0.5 mm R1 / R2 |
|---|---|---|
| as V2 (upwind, uncorrected, legacy) | 1.226 / 1.136 | 1.160 / 1.241 |
| (a) `convection` linearUpwind | 1.166 / 1.111 | 1.057 / 1.143 |
| (b) `sensitivityWallGradient corrected` | 1.311 / 1.177 | 1.245 / 1.332 |
| (a) + (b) | 1.255 / 1.153 | 1.136 / 1.228 |
| (c) `ATCInlets strict` | 1.209 / 1.119 | 1.149 / 1.229 |

**(a) linearUpwind** is the largest single gain on the snapped T-junction and on the manifold, and
the objective path can now select it without touching fvSchemes; it does not bring R2 inside
[0.9, 1.1] (1.143) and R2 still grows with refinement (0.111 -> 0.143). Not made the default
(see below).

**(c) the 1 % inlet rule** moves every ratio by about 1 % (the outlet cells it stops zeroing are
1 % of the domain): not the cause.

**(b) wall gradient / snapped cells, tested on a body-fitted bump.** The same T-junction from
blockMesh alone: no snappy, no bump-to-wall corner; a cos^2 bump pressed into the top wall of each
branch, y_top = 10 mm - a_i bx(x) bz(z), half widths 4 mm in x and z, centred 25 mm from the
junction, a1 = 4 mm, a2 = 2.5 mm (steepest wall 54 to 57 degrees), points moved by -a_i b y/H, so
every mesh has the same topology for every a_i and the finite differences are noise-free (a
1e-5 m step reproduces the 1e-4 m slope to 2e-8 in J). dJ/da_i = sum s (-b_i n_y) A. Split
0.488 / 0.512, J = 3.11e-4 (1 mm) and 2.98e-4 (0.5 mm). Steps +-0.1 mm, central:

| Mesh | cells | dJ/da1 | dJ/da2 | cosine (as V2) |
|---|---|---|---|---|
| 1 mm | 13,000 | 0.4810 | -0.3385 | 0.9997 |
| 0.5 mm | 104,000 | 0.4643 | -0.3230 | 0.99998 |
| 0.33 mm | 351,000 | 0.4648 | -0.3202 | 0.99998 |

| Variant (a1 / a2 ratios) | 1 mm | 0.5 mm | 0.33 mm |
|---|---|---|---|
| as V2 (upwind, uncorrected) | 1.047 / 0.994 | 1.160 / 1.175 | 1.126 / 1.143 |
| (a) linearUpwind | 1.114 / 1.041 | 1.165 / 1.164 | 1.120 / 1.118 |
| (c) strict | 1.035 / 0.982 | 1.152 / 1.165 | |
| (a) + (c) | 1.103 / 1.028 | 1.156 / 1.154 | |
| ATC `standard` | 0.913 / 0.869 | 1.105 / 1.106 | |
| ATC `cancel` | 1.377 / 1.285 | 1.663 / 1.625 | |
| (b) corrected | 0.929 / 0.878 | **1.036 / 1.058** | 1.031 / 1.056 |
| (a) + (b) | 0.982 / 0.914 | **1.032 / 1.042** | **1.021 / 1.030** |

(The finite differences converge: dJ/da2 -0.3385, -0.3230, -0.3202. Iterations 165, 361, 621;
every V1 invariant met on every run, J from phi to 2e-12, Euler 3e-12, adjoint port flux 5e-9.)

- With no snapping and no corner the uncorrected gradient **still does not converge** (1.05 / 0.99
  -> 1.16 / 1.18 -> 1.13 / 1.14), so the snapped cells and the bump-to-wall corner are not what
  keeps R2 away from 1.
- The uncorrected wall gradient is the error on a sloped wall (here up to 57 degrees, cell
  columns vertical): correcting it moves the body-fitted ratios by 10 to 13 %. Upwind + corrected
  gives 1.036 / 1.058 at 0.5 mm with |ratio - 1| halving from 1 mm, then flat at 0.33 mm
  (1.031 / 1.056); linearUpwind + corrected gives 1.032 / 1.042 at 0.5 mm and 1.021 / 1.030 at
  0.33 mm, closing in on 1 (from 1 mm, a1 first moves away, 0.018 -> 0.032). **On the body-fitted
  case, (a) + (b) meets "both in [0.9, 1.1] at 0.5 mm and improving with refinement" from 0.5 mm
  on; (b) alone meets it from 1 to 0.5 mm.**
- On the snapped T-junction the same correction makes every ratio 6 to 9 % worse, and on the
  manifold it helps by 2 %: on snapped polyhedra the cell gradient next to the wall, which the
  correction relies on, is itself poor. It is therefore kept as an option, not a default.
- The ATC's discretisation is a large term: `standard` against `UaGradU` differs by 13 % at 1 mm
  and 5 % at 0.5 mm (both converge to the same continuous term), and dropping it (`cancel`)
  raises the gradient by 30 to 60 %.
- A second-order one-sided wall gradient (two cells along the column, post-processed from the
  written upwind fields, no k correction) gives dJ/da1 0.733 at 1 mm and 0.644 at 0.5 mm (ratios
  1.52, 1.39): without the tangential correction a higher-order normal difference makes it worse,
  which suggests the error on these meshes is mostly the off-normal cell centre rather than the
  order of the difference.

**Decision: no default changes.** None of (a), (b), (c), alone or together, brings both
ratios of the V2 case (the snapped T-junction) inside [0.9, 1.1] at 0.5 mm: the best is (a) at
1.057 / 1.143, and (b), which is what fixes the body-fitted case, makes the snapped one worse
(1.245 / 1.332; with (a) 1.136 / 1.228). So the options stay opt-in, the defaults are the earlier
behaviour, and the case Stokes compiles (golden `adjoint_flow_split`) keeps its meaning. The
measured direction for a fix is a wall gradient that is right on snapped cells as well (a
least-squares wall gradient, or boundary layers so the wall cells are orthogonal), then (a) + (b)
as the objective path's default; that is not done here.

### Verdict

**The flow split's adjoint gradient is not adjoint-verified under ADR-0041's V2/V3 thresholds.**
V2 fails as before (R2 1.241 at 0.5 mm and growing; no variant fixes R2 on the snapped
T-junction). V3 fails as Stokes compiles the case (h2 0.752 < 0.8; h4 inconclusive), and passes
the [0.8, 1.25] band narrowly only with two non-default options that worsen V2. The direction is
right everywhere (cosines 0.994 to 0.99998, no wrong sign), so as ADR-0041 says the adjoint ships,
if at all, `unchecked` and experimental, with finite differences the default; the
`shape.optimise` part of V3 (gradient: adjoint, solves and minutes against 16 and 84) was not run.

## Not done here

- `shape.check` / `shape.optimise` with `gradient: adjoint` on the manifold through Stokes (the
  host has no `adjointFoam` beside the install, and the agent side was being changed).
- A finer snapped T-junction (0.33 mm), and the T-junction with boundary layers, to see whether
  the snapped R2 trend turns.
- A wall gradient that is right on snapped cells (least squares), or ESI's E-SI surface formula
  (roadmap step 4); then (a) + (b) as the objective path's default if V2 and V3 pass with it.
- `startIteration`, parallel runs, turbulent cases.
