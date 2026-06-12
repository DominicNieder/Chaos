using GLMakie 
using CairoMakie
include("../models/henon_heiles.jl")
include("../analysis/poincare.jl")
include("../styles/makie_theme.jl")

savepath = joinpath(@__DIR__, "../../figures/henon-heiles/potential.png")


# creating a plot of the potential

with_theme(QUARTO_THEME) do
    r = range(-0.5, 0.5, length = 120)
    epot = [potential(x,y) for x in r, y in r]

    f = Figure(size = (1000, 400), transparency = true)

    a1 = Axis3(f[1, 1], zlabel="E", zlabeloffset=50)
    a2 = Axis(f[1, 2], xlabel="x", ylabel="y")

    colorscale = ReversibleScale(epot -> epot^(1 / 10), epot -> epot^10)
    levels = logrange(0.009, 0.089, 7)

    contour3d!(a1, r, r, epot, linewidth=3, levels=20,  colormap=:hsv, colorscale=colorscale)
    contour!(a2, r, r, epot, labels=true, levels=levels, colormap=:hsv, colorscale=colorscale)
    
    
    mkpath(dirname(savepath))
    save(savepath, f, px_per_unit=2)

end


