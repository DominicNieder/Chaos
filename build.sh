#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

echo "==> syncing glossary..."
bash "$ROOT/code/scripts/sync_glossary.sh"

echo "==> cleaning previous output..."
rm -rf "$ROOT/_site"

echo "==> rendering notes..."
quarto render "$ROOT"

echo "==> done"