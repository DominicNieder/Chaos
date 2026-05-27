#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TEMPLATES="$ROOT/templates/yy-mm-dd-log"
TOPICS="$ROOT/topics"
INDEX="$ROOT/index.qmd"

# use provided date or today
DATE="${1:-$(date +%y-%m-%d)}"
DIR="$TOPICS/$DATE-log"

# format date: "26-05-27" -> "27 May, 2026"
format_date() {
  local yy mm dd
  IFS='-' read -r yy mm dd <<< "$1"
  date -d "20${yy}-${mm}-${dd}" "+%-d %B, 20%y"
}
DISPLAY="$(format_date "$DATE")"

if [ -d "$DIR" ]; then
  echo "error: $DIR already exists"
  exit 1
fi

# --- create topic directory from templates ---
mkdir -p "$DIR"
for f in _content.qmd _index.qmd _slides.qmd; do
  sed "s/{{DATE}}/$DISPLAY/g" "$TEMPLATES/$f" > "$DIR/${f#_}"
done
echo "created: $DIR"

# --- insert new row after table header separator ---
NEW_ROW="| [$DISPLAY](topics/${DATE}-log/index.qmd) | log entry | [Slides](topics/${DATE}-log/slides.qmd) |"
awk -v row="$NEW_ROW" '
  /^\|---/ && !inserted { print; print row; inserted=1; next }
  { print }
' "$INDEX" > "$INDEX.tmp" && mv "$INDEX.tmp" "$INDEX"

# --- update Latest section ---
awk -v display="$DISPLAY" -v date="$DATE" '
  /^## Latest:/ { print "## Latest: " display; next }
  /\{\{< include topics\/.*-log\/_content\.qmd >\}\}/ {
    print "{{< include topics/" date "-log/_content.qmd >}}"; next
  }
  { print }
' "$INDEX" > "$INDEX.tmp" && mv "$INDEX.tmp" "$INDEX"

echo "updated: $INDEX"
echo "done — open notes/topics/${DATE}-log/_content.qmd to start writing"
