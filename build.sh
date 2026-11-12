#!/bin/bash
# Build the adjoint solver for Windows against a FoamWin (OpenFOAM 14, MinGW-w64) build tree.
#
#   FOAMWIN_SRC=~/foamwin-clean/OpenFOAM-14 bash build.sh
#
# The tree must already be built (its build-all.sh): this links against its libraries. The binary
# goes to bin/adjointFoam.exe; copy it beside the FoamWin install's other solvers to use it.
set -e
here="$(cd "${0%/*}" && pwd)"
src="${FOAMWIN_SRC:-$HOME/foamwin-clean/OpenFOAM-14}"
if [ ! -f "$src/etc/pre-bashrc-mingw_w64.sh" ]; then
    echo "FOAMWIN_SRC ($src) is not a FoamWin OpenFOAM-14 tree" >&2
    exit 2
fi
cd "$src"
# shellcheck disable=SC1091
source etc/pre-bashrc-mingw_w64.sh >/dev/null 2>&1 || true
export WM_THIRD_PARTY_DIR="$(cd .. && pwd)"
export FOAM_EXT_LIBBIN="$WM_THIRD_PARTY_DIR/platforms/linuxmingw-w64GccDPInt32/lib"
echo "WM_OPTIONS=$WM_OPTIONS  FOAM_LIBBIN=$FOAM_LIBBIN"
mkdir -p "$here/bin"
cd "$here/applications/adjointFoam"
wmake
ls -la "$here/bin"
