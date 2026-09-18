using CairoMakie, Colors

# ── Font ──────────────────────────────────────────────────────────────
# Set this to match your LaTeX document's math/text font.
#   pdflatex default (Computer Modern) -> "CMU Serif"       (needs cmu-serif installed)
#   lualatex/xelatex + Latin Modern    -> "Latin Modern Roman"
#   fontspec + libertine/newtxmath     -> "Linux Libertine"
# Makie/CairoMakie resolves system font names via the OS font cache,
# so install the matching TTF/OTF (e.g. `cm-unicode` or `lmodern` fonts
# package) on your machine first.
const FONT = "Latin Modern Roman"
const FONT_BOLD = FONT * " Bold"

# ── Colours ───────────────────────────────────────────────────────────
const BG    = RGBf(1.0, 1.0, 1.0)
const INK   = RGBf(0.0, 0.0, 0.0)   # near-black text/axes for print contrast

# Okabe–Ito palette: the standard colorblind-safe categorical set
# (distinguishable under all common CVD types; also fine in greyscale print)
const C_ORANGE = "#E69F00"
const C_SKY    = "#56B4E9"
const C_GREEN  = "#009E73"
const C_YELLOW = "#F0E442"
const C_BLUE   = "#0072B2"
const C_VERMIL = "#D55E00"
const C_PINK   = "#CC79A7"
const C_BLACK  = "#000000"

const OKABE_ITO = [C_BLUE, C_VERMIL, C_GREEN, C_ORANGE,
                    C_PINK, C_SKY, C_BLACK, C_YELLOW]

# Continuous/sequential data: perceptually uniform + colorblind-safe
const CMAP_SEQ  = :viridis   # single-hue trends, densities
# Diverging data (e.g. signed quantities around zero): Fabio Crameri's
# "vik" is colorblind-safe and print-friendly (avoid classic red-blue jet-likes)
const CMAP_DIV  = :vik

const LATEX_THEME = Theme(
    backgroundcolor = BG,
    textcolor       = INK,
    fonts = (regular = FONT, bold = FONT_BOLD, italic = FONT, bold_italic = FONT_BOLD),
    fontsize = 11,             # match \normalsize (adjust to your document class)
    Axis = (
        backgroundcolor = BG,
        xgridcolor = RGBAf(0, 0, 0, 0.10),
        ygridcolor = RGBAf(0, 0, 0, 0.10),
        xticklabelcolor = INK, yticklabelcolor = INK,
        xlabelcolor = INK,     ylabelcolor = INK,
        titlecolor  = INK,
        spinecolor  = INK,
        xtickcolor  = INK, ytickcolor = INK,
    ),
    Legend = (
        backgroundcolor  = BG,
        framecolor       = RGBAf(0, 0, 0, 0.3),
        labelcolor       = INK,
        titlecolor       = INK,
        patchstrokecolor = RGBAf(0, 0, 0, 0.15),
    ),
    Colorbar = (
        ticklabelcolor = INK,
        labelcolor     = INK,
        topspinecolor = INK, bottomspinecolor = INK,
        leftspinecolor = INK, rightspinecolor = INK,
    ),
    palette = (color = OKABE_ITO,),
)

pick_color(i, cmap) = cmap[mod1(i, length(cmap))]