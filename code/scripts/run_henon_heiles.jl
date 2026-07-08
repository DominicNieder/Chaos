using GLMakie 
using CairoMakie
include("../models/henon_heiles.jl")
include("../styles/makie_theme.jl")

# "two-in-one" file
savepath_2 = joinpath(@__DIR__, "../../figures/henon-heiles/potential/potential-a$(param[1])-m$(param[2])-w$(param[3]).png")

# contour file
savepath_cont = joinpath(@__DIR__, "../../figures/henon-heiles/potential/potential_contour-a$(param[1])-m$(param[2])-w$(param[3]).png")

# 3D file
savepath_3D = joinpath(@__DIR__, "../../figures/henon-heiles/potential/potential_3D-a$(param[1])-m$(param[2])-w$(param[3]).png")


# creating a plot of the potential

with_theme(QUARTO_THEME) do
    set_colourmap  = :hsv # :hsv
    # require a different rnage of vairables and levels for different parameters \alpha
    if param[1] == 2.0
        r = range(-0.5, 0.5, length = 120)
        levels = logrange(0.009, 0.089, 7)
    elseif param[1] == 1.0
        r = range(-1.0, 1.0, length = 120)
        levels = logrange(5.0*0.009,6.9*0.089,7)
    else 
        r = range(-12, 12, length = 120)
        levels = logrange(4*0.009,4*0.089,21)
    end

    println("---- Plotting the Potential----")
    println("*dark theme")
    println("---\nmodel parameters:\na=$(param[1]),\nm=$(param[2]),\nw=$(param[3])\n---")

    epot = [HenonHeiles.potential(x,y, param) for x in r, y in r]


    colorscale = identity # ReversibleScale(epot -> abs(epot)^(1 / 10), epot -> abs(epot)^10)

    println("plot joint fiugre")
    f = Figure(size = (1000, 400), transparency = true)
    a1 = Axis3(f[1, 1], zlabel="E", zlabeloffset=50)
    a2 = Axis(f[1, 2], xlabel="x", ylabel="y")
    contour3d!(a1, r, r, epot, linewidth=3, levels=20,  colormap=set_colourmap, colorscale=colorscale)
    contour!(a2, r, r, epot, labels=true, levels=levels, colormap=set_colourmap, colorscale=colorscale)

    
    mkpath(dirname(savepath_2))
    save(savepath_2, f, px_per_unit=1)

    println("plotting contour standalone")
    f_cont = Figure(size = (800, 800))
    a2_cont = Axis(f_cont[1, 1], xlabel="x", ylabel="y")
    contour!(a2_cont, r, r, epot, labels=true, levels=levels, colormap=set_colourmap, colorscale=colorscale)

    mkpath(dirname(savepath_cont))
    save(savepath_cont, f_cont, px_per_unit=1)

    println("plotting 3D standalone")
    f_3D =   Figure(size = (800, 800))
    a1_3D = Axis3(f_3D[1, 1], zlabel="E", zlabeloffset=50)
    contour3d!(a1_3D, r, r, epot, linewidth=3, levels=20,  colormap=set_colourmap, colorscale=colorscale)

    mkpath(dirname(savepath_3D))
    save(savepath_3D, f_3D, px_per_unit=1)

end


