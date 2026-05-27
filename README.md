# Chaos

I will be studying and exploring the topic of Chaos for at least the next year. Among this I will use various non-linear models and to gain a understanding of chaotic systems. Using analytical and numerically methods as tools in the hope to find patterns and generalized descriptions. In this, the study of the phase space will be crucial in order to gain an understanding of the dynamics.

In order to keep an overview of thoughts and to keep track of my path, notes ([Quarto](https://quarto.org)) will be taken and can be understood here. 

Published here: **https://dominicnieder.github.io/Chaos/**

## Structure

```
Chaos/
├── notes/                  # Quarto website source
│   ├── index.qmd           # Homepage with topic table
│   ├── glossary.json       # Definitions and symbols
│   ├── topics/             # One folder per topic
│   │   └── topic-name/
│   │       ├── _content.qmd   # Shared content
│   │       ├── index.qmd      # Web page view
│   │       └── slides.qmd     # Slide deck view
│   └── _site/              # Rendered output (ignored by git)
├── code/                   # Julia source
│   ├── models/             # System definitions (equations of motion)
│   │   └── henon_heiles.jl
│   ├── analysis/           # Derived quantities
│   │   ├── poincare.jl     # Poincaré section
│   │   └── lyapunov.jl     # Lyapunov exponent
│   ├── scripts/
│   │   ├── simulate.jl     # Run + save trajectory to data/
│   │   └── explore.jl      # Interactive GLMakie dashboard
│   ├── notebooks/
│   │   └── henon_heiles.jl # i.e. Pluto reactive notebook
│   └── Project.toml        # Julia environment
├── figures/                # Generated plots (ignored by git)
├── data/                   # Simulation outputs and parameter sets
│   └── orientation.json    # data orientation notes
├── lit/                    # Literature and bibliography
└── to_read/
```

## Code & visualisation

The Julia codebase is split into three layers:

- **`models/`** — equations of motion; import these everywhere
- **`analysis/`** — Poincaré sections, Lyapunov exponents, etc.
- **`scripts/explore.jl`** — interactive GLMakie dashboard with sliders for energy and initial conditions; plots update live
- **`notebooks/`** — Pluto reactive notebooks for structured exploration

Figures are saved to `figures/` and referenced in notes with relative paths. Data files are saved to `data/` and shared between code and notes.

To run the interactive dashboard:
```julia
julia --project=code code/scripts/explore.jl
```

To open the Pluto notebook:
```julia
julia --project=code -e 'using Pluto; Pluto.run()'
```

