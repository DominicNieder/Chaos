#!/usr/bin/env bash
set -euo pipefail

if [[ "${RENDERING_SLIDE_DECK:-}" == "1" ]]; then
  exit 0
fi

export RENDERING_SLIDE_DECK=1
quarto render topics/two_trpapped_pointCharges/slides.qmd \
  --to revealjs \
  --output slides-revealjs.html \
  --output-dir _site/topics/two_trpapped_pointCharges
