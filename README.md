# Chaos

A personal knowledge base rendered as a website using [Quarto](https://quarto.org).

Published at: **https://dominicnieder.github.io/Chaos/**

## Structure

```
Chaos/
├── notes/          # Quarto website source
│   ├── index.qmd   # Homepage with topic table
│   ├── topics/     # One folder per topic
│   │   └── topic-name/
│   │       ├── _content.qmd   # Shared content (single source of truth)
│   │       ├── index.qmd      # Web page view
│   │       └── slides.qmd     # Slide deck view
│   └── _site/      # Rendered output (ignored by git)
├── data/
├── lit/
├── pic/
└── to_read/
```

## Notes & website

Each topic lives in `notes/topics/<topic-name>/`. Content is written once in `_content.qmd` and included by both the reading page and the slide deck, so edits propagate to both outputs automatically.

To preview locally:

```bash
cd notes
quarto preview
```

To build the full site:

```bash
cd notes
quarto render
```

The site is published automatically to GitHub Pages from the `_site/` output.
