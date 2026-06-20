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

# --- styles/glossary.html ---
GLOSSARY_JS="$ROOT/styles/glossary.html"
{
  echo "<script>"
  echo "document.addEventListener('DOMContentLoaded', () => {"
  echo "  const terms = {"
  jq -r 'to_entries[]
    | select(.value | type == "object")
    | select(.value.def != null and .value.def != "")
    | "    " + (.key | @json) + ": " + (.value.def | @json) + ","' "$GLOSSARY"
  echo "  };"
  cat << 'EOF'
  const walker = document.createTreeWalker(
    document.querySelector("main") || document.body,
    NodeFilter.SHOW_TEXT
  );
  const sorted = Object.keys(terms).sort((a, b) => b.length - a.length);
  const nodes = [];
  while (walker.nextNode()) nodes.push(walker.currentNode);
  nodes.forEach(node => {
    if (node.parentElement.closest("code, pre, script, style, abbr")) return;
    let html = node.textContent.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
    let text = node.textContent;
    let result = text;
    sorted.forEach(term => {
      const re = new RegExp("\\b(" + term.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + ")\\b", "gi");
      result = result.replace(re, `<abbr title="${terms[term]}">$1</abbr>`);
    });
    if (result !== text) {
      const span = document.createElement("span");
      span.innerHTML = result;
      node.parentElement.replaceChild(span, node);
    }
  });
});
EOF
  echo "</script>"
} > "$GLOSSARY_JS"

echo "sync_glossary: macros.html, glossary.qmd, and glossary.html updated"
