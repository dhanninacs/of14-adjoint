# ATC choices on a turbulent case (measured 2026-11-12, FoamWin)

The OpenFOAM 14 `pitzDaily` tutorial of `adjointShapeOptimisationFoam` (k-epsilon, frozen in the
adjoint, U = 10 m/s, nu = 1e-5), with the porosity off, 1000 SIMPLE iterations, run with
`bin/adjointFoam.exe` under the FoamWin environment. Stable is a bounded Ua (the inlet sets
|Ua| = 1) and a small adjoint continuity error.

| `ATC` in `system/adjointDict` | max \|Ua\| at 1000 | adjoint continuity error (sum local) | |
|---|---|---|---|
| `{ type UaGradU; }` (the default: the Foundation solver) | 1.0e5 | 404 | diverges, identical to the Foundation binary |
| `{ type UaGradU; zeroATCPatchTypes (patch); }` | 6.35 | 0.026 | stable |
| `{ type standard; zeroATCPatchTypes (patch); }` | 30.9 | 0.058 | bounded, large |
| `{ type standard; zeroATCPatchTypes (wall patch); nSmooth 1; }` | 4.41 | 1.6e-4 | stable |
| `{ type cancel; }` | 3.43 | 3.0e-5 | most stable, least accurate |

**What makes the difference is ESI's limiter**: zeroing the adjoint transpose convection in the
cells next to the inlets and outlets (`patch` type) stops the turbulent adjoint diverging, with
either form of the term. This was first seen by accident: an early inlet rule also matched the
outlet's `inletOutlet` condition, and a bisection against the Foundation binary traced the
difference to it (the rule now excludes mixed conditions).

Regression: with the default, the laminar obstacle case of the Stokes tests gives the same
objective and the same wall sensitivity, to every printed digit, as before the ATC port, and the
turbulent case reproduces the Foundation binary's divergence exactly.

Not yet measured: how each choice changes the **accuracy** of the gradient against finite
differences (`cancel` is known to bias it). That is the next test, on a case with a turbulent
finite-difference reference.
