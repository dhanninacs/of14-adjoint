# Provenance

Every source file keeps the copyright and licence header of the code it was taken from, and adds a
line naming this project's changes. All of it is GPL-3.0-or-later.

| Path | From | Copyright |
|---|---|---|
| `applications/adjointFoam/adjointFoam.C`, `createFields.H`, `createPhia.H`, `adjointContinuityErrs.H`, `initAdjointContinuityErrs.H` | OpenFOAM 14 (Foundation), `applications/legacy/incompressible/adjointShapeOptimisationFoam` | OpenFOAM Foundation |
| `applications/adjointFoam/adjointOutletPressure`, `adjointOutletVelocity` | the same | OpenFOAM Foundation |
| `applications/adjointFoam/ATC.H` | ESI-OpenCFD OpenFOAM v2xxx, `src/optimisation/adjointOptimisation/adjoint/ATCModel` (`ATCModel`, `ATCstandard`, `ATCUaGradU`, `cancelATC`, `zeroATCcells/faceCells`), collapsed into one header and adapted to the Foundation API | PCOpt/NTUA (2007-2023), FOSS GP (2013-2023), OpenCFD Ltd. (2019-2021) |
| `applications/adjointFoam/adjointControls.H`, `shapeSensitivity.H` | of14-adjoint | of14-adjoint contributors |

ESI's sources were read from `https://gitlab.com/openfoam/core/openfoam` (the commit is
recorded in ROADMAP.md).
