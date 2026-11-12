# of14-adjoint

A continuous-adjoint shape-sensitivity solver for **OpenFOAM 14** (the OpenFOAM Foundation
release), built for Windows with the FoamWin MinGW-w64 toolchain and on Linux with the ordinary
`wmake`. It starts from the Foundation's own `adjointShapeOptimisationFoam` and brings in pieces of
ESI-OpenCFD's adjoint optimisation library (OpenFOAM v2xxx, `src/optimisation`), which the
Foundation release does not have.

This is an independent community project. It is not part of, approved by or endorsed by the
OpenFOAM Foundation or OpenCFD Ltd. OPENFOAM® is a registered trade mark of OpenCFD Limited.

## What it does

`adjointFoam` solves a steady incompressible flow (SIMPLE) and its continuous adjoint for the
total-pressure loss between inlets and outlets (Othmer 2008), and writes:

- `U`, `p`, `Ua`, `pa` as usual;
- **`shapeSensitivity`** on every wall patch: the change of the objective per unit area and per unit
  displacement of the wall along its normal pointing out of the fluid, from OpenFOAM's own boundary
  normal gradients;
- `postProcessing/adjoint/0/objective.dat`: the objective each iteration.

Compared with the Foundation solver it starts from:

| | Foundation `adjointShapeOptimisationFoam` | `adjointFoam` |
|---|---|---|
| Purpose | topology optimisation by porosity | shape sensitivity (porosity optional, off by default) |
| Adjoint transpose convection | always the full term (diverges on a turbulent pitzDaily without porosity) | `ATC` choice from ESI: `UaGradU` (default, the same term), `standard`, `cancel`, with ESI's limiter next to chosen patch types: `zeroATCPatchTypes (patch)` makes that case stable ([tests/atc.md](tests/atc.md)) |
| Wall sensitivity | none | `shapeSensitivity` written per wall face |
| Inlets | a patch named `inlet` | every non-wall patch with a plain fixed-value U |
| Inlet adjoint velocity | as given | as given (the pressure-loss objective wants `-U_in`) |

Settings: `system/adjointDict`, every entry optional ([etc/adjointDict](etc/adjointDict)). The plan for the ESI merge, in order, is in [ROADMAP.md](ROADMAP.md).

## Build

```bash
# Windows binaries, from Linux or WSL, against a built FoamWin tree:
FOAMWIN_SRC=~/foamwin-clean/OpenFOAM-14 bash build.sh      # -> bin/adjointFoam.exe
```

Copy `bin/adjointFoam.exe` beside the install's other solvers. On Linux, source an OpenFOAM 14
`etc/bashrc` and run `wmake` in `applications/adjointFoam` with `EXE` set to taste.

## Licence and provenance

GPL-3.0-or-later, as OpenFOAM itself ([LICENSE](LICENSE)). Every file keeps the copyright header of
the code it came from; [NOTICE.md](NOTICE.md) lists which file came from where.
