const C_RED    = "#f48484"
const C_PURPLE = "#ae87f7"
const C_CREAM  = "#f8f2e0"
const C_TEAL   = "#7ec8c8"

const QUARTO_THEME = Theme(
    backgroundcolor = RGBf(0.102, 0.102, 0.102),   # #1a1a1a
    textcolor       = RGBf(0.973, 0.949, 0.878),   # #f8f2e0
    Axis = (
        backgroundcolor = RGBf(0.102, 0.102, 0.102),
        xgridcolor  = RGBAf(1,1,1,0.08),
        ygridcolor  = RGBAf(1,1,1,0.08),
        xticklabelcolor = RGBf(0.973, 0.949, 0.878),
        yticklabelcolor = RGBf(0.973, 0.949, 0.878),
        xlabelcolor = RGBf(0.973, 0.949, 0.878),
        ylabelcolor = RGBf(0.973, 0.949, 0.878),
        titlecolor  = RGBf(0.980, 0.973, 0.965),   # #f9f8f6
        spinecolor  = RGBf(0.973, 0.949, 0.878),
    ),
    palette = (color = [C_RED, C_PURPLE, C_CREAM, C_TEAL],),

)
