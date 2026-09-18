using GLMakie, Colors, ColorSchemes
using CairoMakie
# =========================================================================
# FONTS
# =========================================================================
const FONT_LATEX = "Latin Modern Roman"    # pdflatex+lmodern / xelatex+fontspec
# const FONT_LATEX = "CMU Serif"           # plain pdflatex (Computer Modern)
const FONT_TYPST = "New Computer Modern"   # Typst default
const FONT_DARK  = "CMU Serif"

const DOC_FONT = FONT_LATEX   # <- switch this one line for LaTeX vs Typst

# =========================================================================
# DARK THEME (screen / day-to-day)
# =========================================================================
const BG_DARK = RGBf(0.102, 0.102, 0.102)
const CREAM   = RGBf(0.973, 0.949, 0.878)

const PALETTE_DARK = ["#f48484", "#7ec8c8", "#e8c97d", "#ae87f7", "#9ed07a",
                       "#f0a56b", "#7ea8e8", "#e88ec5", "#b6bcc6", "#f8f2e0"]
const CMAP_SEQ_DARK = :viridis
const CMAP_DIV_DARK = :vik

const DARK_THEME = Theme(
    backgroundcolor = BG_DARK, textcolor = CREAM,
    fonts = (regular = FONT_DARK, bold = FONT_DARK * " Bold",
             italic = FONT_DARK, bold_italic = FONT_DARK * " Bold"),
    fontsize = 14,
    Axis = (backgroundcolor = BG_DARK,
        xgridcolor = RGBAf(1,1,1,0.08), ygridcolor = RGBAf(1,1,1,0.08),
        xticklabelcolor = CREAM, yticklabelcolor = CREAM,
        xlabelcolor = CREAM, ylabelcolor = CREAM,
        titlecolor = RGBf(0.980, 0.973, 0.965), spinecolor = CREAM),
    Legend = (backgroundcolor = BG_DARK, framecolor = RGBAf(1,1,1,0.15),
        labelcolor = CREAM, titlecolor = CREAM, patchstrokecolor = RGBAf(1,1,1,0.15)),
    Colorbar = (ticklabelcolor = CREAM, labelcolor = CREAM,
        topspinecolor = CREAM, bottomspinecolor = CREAM,
        leftspinecolor = CREAM, rightspinecolor = CREAM),
    palette = (color = PALETTE_DARK,),
)

# =========================================================================
# PRINT THEME (LaTeX / Typst) — white bg, colorblind-safe (Okabe–Ito)
# =========================================================================
const BG_PRINT = RGBf(1.0, 1.0, 1.0)

const PALETTE_PRINT = ["#0072B2", "#D55E00", "#009E73", "#E69F00",
                        "#CC79A7", "#56B4E9", "#000000", "#F0E442"]
const CMAP_SEQ_PRINT = :viridis
const CMAP_DIV_PRINT = :vik

const LATEX_THEME = Theme(
    backgroundcolor = BG_PRINT, textcolor = RGBf(0,0,0),
    fonts = (regular = DOC_FONT, bold = DOC_FONT * " Bold",
             italic = DOC_FONT, bold_italic = DOC_FONT * " Bold"),
    fontsize = 11,
    Axis = (backgroundcolor = BG_PRINT,
        xgridcolor = RGBAf(0,0,0,0.10), ygridcolor = RGBAf(0,0,0,0.10),
        xticklabelcolor = :black, yticklabelcolor = :black,
        xlabelcolor = :black, ylabelcolor = :black,
        titlecolor = :black, spinecolor = :black,
        xtickcolor = :black, ytickcolor = :black),
    Legend = (backgroundcolor = BG_PRINT, framecolor = RGBAf(0,0,0,0.3),
        labelcolor = :black, titlecolor = :black, patchstrokecolor = RGBAf(0,0,0,0.15)),
    Colorbar = (ticklabelcolor = :black, labelcolor = :black,
        topspinecolor = :black, bottomspinecolor = :black,
        leftspinecolor = :black, rightspinecolor = :black),
    palette = (color = PALETTE_PRINT,),
)

# =========================================================================
# LIVE GLOBALS — these are what your scripts actually reference
# (`COLOR_SCHEME`, `INK`, `pick_color`), and they follow whichever
# theme is active. NOT `const`, so `set_style!` can repoint them.
# =========================================================================
COLOR_SCHEME = PALETTE_DARK
INK          = CREAM
CMAP_SEQ     = CMAP_SEQ_DARK
CMAP_DIV     = CMAP_DIV_DARK

"""
    set_style!(:dark)   -- day-to-day screen work
    set_style!(:print)  -- LaTeX/Typst paper figures
"""
function set_style!(name::Symbol)
    global COLOR_SCHEME, INK, CMAP_SEQ, CMAP_DIV
    if name == :dark
        GLMakie.activate!()
        set_theme!(DARK_THEME)
        COLOR_SCHEME, INK = PALETTE_DARK, CREAM
        CMAP_SEQ, CMAP_DIV = CMAP_SEQ_DARK, CMAP_DIV_DARK
    elseif name == :print
        CairoMakie.activate!()
        set_theme!(LATEX_THEME)
        COLOR_SCHEME, INK = PALETTE_PRINT, RGBf(0,0,0)
        CMAP_SEQ, CMAP_DIV = CMAP_SEQ_PRINT, CMAP_DIV_PRINT
    else
        error("unknown style $name, use :dark or :print")
    end
    return name
end

"""
    pick_color(i)                    -- i-th color of active COLOR_SCHEME (cycles)
    pick_color(i, cmap::AbstractVector) -- i-th color of a given categorical palette
    pick_color(t, :seq | :div)       -- t ∈ [0,1] along active sequential/diverging colormap
    pick_color(i, cmap::Symbol)      -- i-th of 8 samples from a named ColorSchemes map
"""
pick_color(i::Integer) = COLOR_SCHEME[mod1(i, length(COLOR_SCHEME))]
pick_color(i::Integer, cmap::AbstractVector) = cmap[mod1(i, length(cmap))]

function pick_color(t::Real, kind::Symbol)
    cmap = kind == :seq ? CMAP_SEQ : kind == :div ? CMAP_DIV : kind
    return get(colorschemes[cmap], clamp(t, 0.0, 1.0))
end
pick_color(i::Integer, cmap::Symbol) = get(colorschemes[cmap], (i - 1) / 7)

set_style!(:print)   # this script's default (matches your existing set_theme!(LATEX_THEME) call)