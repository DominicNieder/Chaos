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
    palette = (color = [
        "#f48484",   # link red — first series
        "#ae87f7",   # code purple — second series
        "#f8f2e0",   # cream — third
        "#7ec8c8",   # teal — fourth
    ],),
)
