#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
GLOSSARY="$ROOT/notes/glossary.json"
MACROS="$ROOT/notes/styles/macros.html"
QMD="$ROOT/notes/glossary.qmd"

# --- styles/macros.html ---
{
  echo "<script>"
  echo "MathJax = { tex: { macros: {"
  jq -r 'to_entries[]
    | select(.value.macro != null and .value.macro != "")
    | "  \(.value.macro): \"\(.value.latex)\","' "$GLOSSARY"
  echo "}}};"
  echo "</script>"
} > "$MACROS"

# --- notes/glossary.qmd ---
{
  echo "---"
  echo "title: Glossary"
  echo "---"
  echo ""
  echo "| Term | Definition |"
  echo "|------|------------|"
  jq -r 'to_entries[]
    | "| \(.key) | \(.value.def) |"' "$GLOSSARY"
} > "$QMD"

echo "sync_glossary: macros.html and glossary.qmd updated"
