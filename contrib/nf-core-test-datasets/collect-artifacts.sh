#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${ROOT}/contrib/nf-core-test-datasets/build_output"
DEST="${ROOT}/contrib/nf-core-test-datasets/extension_base"

mkdir -p "${DEST}"

cp "${OUT}/combine/integrate/scvi/scvi_model/model.pt" "${DEST}/model.pt"
cp "${OUT}/finalized/merged.h5ad" "${DEST}/merged.h5ad"
cp "${OUT}/combine/integrate/harmony/harmony_reference.h5ad" "${DEST}/harmony_reference.h5ad"

echo "Collected artifacts into ${DEST}:"
ls -lh "${DEST}"
