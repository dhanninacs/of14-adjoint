# Provenance

Every source file keeps the copyright and licence header of the code it was taken from, and adds a
line naming this project's changes. All of it is GPL-3.0-or-later.

| Path | From | Copyright |
|---|---|---|
| `applications/adjointFoam/adjointFoam.C`, `createFields.H`, `createPhia.H`, `adjointContinuityErrs.H`, `initAdjointContinuityErrs.H` | OpenFOAM 14 (Foundation), `applications/legacy/incompressible/adjointShapeOptimisationFoam` | OpenFOAM Foundation |
| `applications/adjointFoam/adjointOutletPressure`, `adjointOutletVelocity` | the same | OpenFOAM Foundation |
| `applications/adjointFoam/ATC.H` | ESI-OpenCFD OpenFOAM v2xxx, `src/optimisation/adjointOptimisation/adjoint/ATCModel` (`ATCModel`, `ATCstandard`, `ATCUaGradU`, `cancelATC`, `zeroATCcells/faceCells`), collapsed into one header and adapted to the Foundation API | PCOpt/NTUA (2007-2023), FOSS GP (2013-2023), OpenCFD Ltd. (2019-2021) |
| `applications/adjointFoam/objectiveOutletPressure`, `objectiveOutletVelocity` | class layout from the Foundation's `adjointOutletPressure`/`adjointOutletVelocity` (above); the boundary values ported from ESI-OpenCFD OpenFOAM v2xxx, `src/optimisation/adjointOptimisation/adjoint/adjointBoundaryConditions/adjointOutletPressure` and `adjointOutletVelocity` (the normal and Robin tangential forms, the `UaGradU` term), with the objective source's sign flipped for this solver's convention (Ua = -u, pa = -q), ESI's explicit deviatoric flux left out, and no tangential source | OpenFOAM Foundation; PCOpt/NTUA (2007-2023), FOSS GP (2013-2023), OpenCFD Ltd. (2019-2023) |
| `applications/adjointFoam/objective.H`, `objectiveUpdate.H` | of14-adjoint, after ESI-OpenCFD's `objectives/incompressible/objectiveFlowRatePartition` (keywords `outletPatches`, `targetFractions`, the patch-index ordering, the default 1/N targets); J is defined as Stokes' flow split (no 1/2, normalised by the summed outlet flux), so dJ/dQ is derived afresh, not copied | of14-adjoint contributors; ESI-OpenCFD as above |
| `applications/adjointFoam/adjointControls.H`, `shapeSensitivity.H` | of14-adjoint | of14-adjoint contributors |

ESI's sources were read from `https://gitlab.com/openfoam/core/openfoam` (the commit is
recorded in ROADMAP.md). The flow-rate objective and its outlet conditions were read at the same
commit, `6e2595746a396d688edd9890f4aebaaa7ff4ecfd`, from the raw files of
`src/optimisation/adjointOptimisation/adjoint/objectives/incompressible/objectiveFlowRatePartition/`
and `src/optimisation/adjointOptimisation/adjoint/adjointBoundaryConditions/adjointOutlet{Pressure,Velocity}/`.
