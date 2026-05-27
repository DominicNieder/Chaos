#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"

echo "==> syncing glossary..."
bash "$ROOT/code/scripts/sync_glossary.sh"

echo "==> rendering notes..."
quarto render "$ROOT/notes/"

echo "==> done"
