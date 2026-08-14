const BG    = RGBf(0.102, 0.102, 0.102)
const CREAM = RGBf(0.973, 0.949, 0.878)
const C_RED    = "#f48484"
const C_PURPLE = "#ae87f7"
const C_CREAM  = "#f8f2e0"
const C_TEAL   = "#7ec8c8"
const C_ORANGE = "#f0a56b"
const C_GOLD   = "#e8c97d"
const C_GREEN  = "#9ed07a"
const C_BLUE   = "#7ea8e8"
const C_PINK   = "#e88ec5"
const C_GREY   = "#b6bcc6"

const QUARTO_THEME = Theme(
    backgroundcolor = BG,
    textcolor       = CREAM,
    Axis = (
        backgroundcolor = BG,
        xgridcolor = RGBAf(1,1,1,0.08),
        ygridcolor = RGBAf(1,1,1,0.08),
        xticklabelcolor = CREAM, yticklabelcolor = CREAM,
        xlabelcolor = CREAM,     ylabelcolor = CREAM,
        titlecolor  = RGBf(0.980, 0.973, 0.965),
        spinecolor  = CREAM,
    ),
    Legend = (
        backgroundcolor  = BG,
        framecolor       = RGBAf(1,1,1,0.15),
        labelcolor       = CREAM,
        titlecolor       = CREAM,
        patchstrokecolor = RGBAf(1,1,1,0.15),
    ),
    Colorbar = (
        ticklabelcolor = CREAM,
        labelcolor     = CREAM,
        topspinecolor = CREAM, bottomspinecolor = CREAM,
        leftspinecolor = CREAM, rightspinecolor = CREAM,
    ),
    palette = (color = [C_RED, C_TEAL, C_GOLD, C_PURPLE, C_GREEN,
                        C_ORANGE, C_BLUE, C_PINK, C_GREY, C_CREAM],),
)