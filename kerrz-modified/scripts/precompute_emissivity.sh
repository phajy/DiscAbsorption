#!/usr/bin/env bash
# Precompute lamppost emissivity FITS files for a grid of spins and heights.
# Untracked helper for local Gradus / library work — not part of the kerrz build.
#
# Usage (from repo root):
#   ./scripts/precompute_emissivity.sh
#
# Edit SPINS / HEIGHTS / NPHOTONS below as needed.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

OUT_DIR="${OUT_DIR:-emissivity}"
# Note: kerrz currently panics for exact a=0 (Carlson RF domain). Use a small
# positive spin (e.g. 0.001) as a Schwarzschild stand-in — same clamp as the C API.
SPINS=(0.001 0.5 0.9 0.998)
HEIGHTS=(3 5 10 20)
NPHOTONS="${NPHOTONS:-3000}"
PHOTON_INDEX="${PHOTON_INDEX:-2.0}"

mkdir -p "$OUT_DIR"

echo "Building kerrz CLI..."
zig build

KERRZ="./zig-out/bin/kerrz"
if [[ ! -x "$KERRZ" ]]; then
  echo "error: expected $KERRZ after zig build" >&2
  exit 1
fi

for spin in "${SPINS[@]}"; do
  for height in "${HEIGHTS[@]}"; do
    out="${OUT_DIR}/emis_a${spin}_h${height}.fits"
    echo "=== spin=${spin} height=${height} -> ${out}"
    "$KERRZ" emissivity \
      --spin "$spin" \
      --lamppost "h:${height},vr:0" \
      --photon-index "$PHOTON_INDEX" \
      --nphotons "$NPHOTONS" \
      -o "$out"
  done
done

echo "Done. Files in ${OUT_DIR}/"
ls -la "$OUT_DIR"
