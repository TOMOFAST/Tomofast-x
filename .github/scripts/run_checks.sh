#!/bin/bash
#
# Runs the Tomofast-x unit tests and the Mansfield example inversion with a packaged
# binary, in a scratch folder. Exits non-zero if anything fails.
#
# Usage: run_checks.sh <path-to-bundle-folder> <executable-name>
#
set -euo pipefail

BUNDLE="$(cd "$1" && pwd)"
EXE_NAME="$2"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

WORK="$(mktemp -d)"
cd "$WORK"

echo "=== Unit tests (OMP_NUM_THREADS=$OMP_NUM_THREADS) ==="
# ftnunit restarts the program after a test that stops it; ftnunit.lst exists until all tests ran.
echo ALL > ftnunit.run
for run in $(seq 1 10); do
  "$BUNDLE/$EXE_NAME" 2>&1 | tee "tests_$run.log"
  [ -f ftnunit.lst ] || break
done
if [ -f ftnunit.lst ]; then
  echo "Unit tests did not complete after $run runs."
  exit 1
fi
if ! grep -h "Number of failed assertions" tests_*.log | grep -Eq ':[[:space:]]+0[[:space:]]*$'; then
  echo "Unit tests reported failures."
  exit 1
fi
rm -f ftnunit.*

echo "=== Example inversion (Mansfield slice) ==="
cp -R "$BUNDLE/parfiles" "$BUNDLE/data" .
"$BUNDLE/$EXE_NAME" -p parfiles/Parfile_mansf_slice.txt 2>&1 | tee example.log
grep -q "THE END" example.log || { echo "Example did not finish."; exit 1; }
[ -s output/mansf_slice/Paraview/grav_final_model3D_full.vtk ] || { echo "Final model missing."; exit 1; }

echo "All checks passed."
