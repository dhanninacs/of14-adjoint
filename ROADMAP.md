# Roadmap: merging ESI's adjoint into OpenFOAM 14

ESI-OpenCFD's adjoint optimisation library (`src/optimisation/adjointOptimisation/adjoint` in
OpenFOAM v2xxx; about 85,000 lines in 185 source files) was read at commit
`6e2595746a396d688edd9890f4aebaaa7ff4ecfd` (2026-07-28) of `https://gitlab.com/openfoam/core/openfoam`.
It is ported a piece at a time, each piece measured before the next. The Foundation and ESI code
bases differ in their momentum-transport, boundary-condition and solver-control APIs, so every
piece is adapted, not copied.

| Step | What | Status |
|---|---|---|
| 1 | Fork Foundation's `adjointShapeOptimisationFoam`: porosity optional, inlets found by condition, `shapeSensitivity` on walls from the boundary snGrad, objective logged | done |
| 2 | ESI's adjoint transpose convection: `UaGradU` (Foundation's form, default), `standard`, `cancel`; limiter next to patch types with smoothing; `extraConvection` (`ATC.H`) | done: the limiter next to inlets and outlets makes the turbulent pitzDaily adjoint stable ([tests/atc.md](tests/atc.md)); accuracy against finite differences still to measure |
| 3 | ESI's objectives as adjoint boundary conditions: total-pressure losses (`objectivePtLosses`), force and moment (`objectiveForce`, `objectiveMoment`), flow uniformity (`objectiveUniformityPatch`); the outlet and inlet adjoint conditions they imply (`adjointOutletVelocity`, `adjointOutletPressure`, `adjointInletVelocity` of ESI) | partly done: the flow-rate partition (`objective { type flowRatePartition; }`, ESI's `objectiveFlowRatePartition` with its outlet conditions as `objectiveOutletPressure`/`objectiveOutletVelocity`), laminar, measured against finite differences ([tests/flowsplit.md](tests/flowsplit.md)): right signs, cosine 0.9999, ratios 1.16 / 1.24 at 0.5 mm (high; R2 does not improve with refinement, so not yet accepted). Still to do: run-time-selectable objective classes, `objectivePtLosses` with ESI's outlet source, force/moment/uniformity, inlet-side sources |
| 4 | ESI's E-SI surface sensitivity (`sensitivitySurface`), including the ATC and objective terms, in place of the frozen-turbulence formula | |
| 5 | Differentiated Spalart–Allmaras (`adjointSpalartAllmaras`), so turbulent cases have proper gradients rather than frozen turbulence | |
| 6 | Tests against ESI's own tutorials (`tutorials/incompressible/adjointOptimisationFoam/sensitivityMaps`) run in an ESI build, and against finite differences | |

## Step 1: measured

On a laminar channel with a spherical obstacle (the Stokes test case), the solver's objective and
its `shapeSensitivity` agree with Stokes' own reconstruction from the written fields: the same J to
every printed digit, face ratios 0.97 to 1.04 (median 1.000), integrals within 0.4 %.
