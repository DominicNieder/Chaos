#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TEMPLATES="$ROOT/templates/_topic_name"
TOPICS="$ROOT/topics"
QUARTO="$ROOT/_quarto.yml"

# topic name from argument
if [ $# -eq 0 ]; then
  echo "usage: bash new_topic.sh <topic-name>"
  echo "example: bash new_topic.sh self-organized-criticality"
  exit 1
fi

NAME="$1"
TITLE="${NAME//-/ }"
TITLE="$(echo "$TITLE" | awk '{for(i=1;i<=NF;i++) $i=toupper(substr($i,1,1)) substr($i,2); print}')"
DIR="$TOPICS/$NAME"

if [ -d "$DIR" ]; then
  echo "error: $DIR already exists"
  exit 1
fi

# --- create topic directory from templates ---
mkdir -p "$DIR"
for f in _index.qmd _slides.qmd; do
  sed "s/{{TITLE}}/$TITLE/g" "$TEMPLATES/$f" > "$DIR/${f#_}"
done
sed "s/{{TITLE}}/$TITLE/g" "$TEMPLATES/_content.qmd" > "$DIR/_content.qmd"
echo "created: $DIR"

# --- add to Topics section in _quarto.yml ---
awk -v name="$NAME" '
  /section: "Topics"/ { print; in_topics=1; next }
  in_topics && /contents:/ { print; print "          - topics/" name "/index.qmd"; in_topics=0; next }
  { print }
' "$QUARTO" > "$QUARTO.tmp" && mv "$QUARTO.tmp" "$QUARTO"

echo "updated: $QUARTO"
echo "done — open topics/$NAME/_content.qmd to start writing"
