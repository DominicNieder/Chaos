#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GLOSSARY="$ROOT/glossary.json"
MACROS="$ROOT/styles/macros.html"
QMD="$ROOT/glossary.qmd"

# --- styles/macros.html ---
{
  echo "<script>"
  echo "MathJax = { tex: { macros: {"
  jq -r 'to_entries[]
    | select(.value | type == "object")
    | select(.value.macro != null and .value.macro != "")
    | "  \(.value.macro): \"\(.value.latex)\","' "$GLOSSARY" | sed 's/\\/\\\\/g'
  echo "}}};"
  echo "</script>"
} > "$MACROS"

# --- notes/glossary.qmd ---
{
  echo "| Term | Definition |"
  echo "|------|------------|"
  jq -r 'to_entries[]
    | select(.value | type == "object")
    | "| <abbr title=\"\(.value.def)\">\(.key)</abbr> | \(.value.def) |"' "$GLOSSARY"
} > "$QMD"

echo "sync_glossary: macros.html and glossary.qmd updated"
