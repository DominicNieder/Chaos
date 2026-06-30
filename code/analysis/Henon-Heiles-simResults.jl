import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GLMakie, JLD2
include("../styles/makie_theme.jl")
include("../models/henon_heiles.jl")
using .HenonHeiles

DATA_DIR        = joinpath(@__DIR__, "../../data/henon-heiles/simulation")
FIG_DIR         = joinpath(@__DIR__, "../../figures/henon-heiles/simulation")

FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

EXPL_DATA_DIR = joinpath(@__DIR__, "../../data/henon-heiles/explore")

json = JSON3.read(read(DATA_ORIENTATION, String))

run = json[:"(E,y)100x256_p00.0/"]
file = map(i -> String(i), run.file)


println(typeof(file))
energies = map(i -> Float64(i),run.energy)
println(typeof(energies))

println("data directory: $(isdir(DATA_DIR)),\nexplore directory: $(isdir(EXPL_DATA_DIR)),\nfigures directory: $(isdir(FIG_DIR))")

fsize=(1200,1200)

# --- contour for E=0.1243 ---
e = 0.1243
param = (1,1,1)
n=200
function section_set_plot_lim(e,param,n)
    sec_boundry = HenonHeiles.section_boundary_ranges(e, param, n)
    y_max, py_max = sec_boundry[1], sec_boundry[2]
    boundary = HenonHeiles.section_boundary(y_max, py_max)
    ylim = extrema(y_max)
    ydiff = (ylim[2]-ylim[1])/50
    ylim = ylim[1]-ydiff, ylim[2]+ydiff 
    p_max_value = maximum(py_max)
    pydiff = p_max_value/50
    pylim = (-p_max_value-pydiff, p_max_value+pydiff)
    (boundary, ylim, pylim)
end



file1 = joinpath(DATA_DIR, "E0.1243-T10000.0-py0.0-n128.jld2")
data  = load(file1)
results = data["results"]
resoltuion = length(data["results"])
println(typeof(results), keys(results[1]))
resolution = 128
fname   = "E$(round(e,digits=4))-T$(10000)-py0.0-n128.png"
outfig = joinpath(FIG_DIR, fname)

function plot_phasespaceof!(ax, e, param, data, fig_size)
    b1, ylim1, pylim1 = HenonHeiles.section_set_plot_lim(e,param,120)
    colors = cgrad(:viridis, resolution)[range(0, 1, length=resolution)]
    for i in eachindex(results)
        scatter!(ax, results[i][:sec_y], results[i][:sec_py], color=colors[i], markersize=3)
    end
    ylims!(ax, pylim1...)
    xlims!(ax, ylim1...)

end

with_theme(QUARTO_THEME) do
    fig1 = Figure(size=fig_size)
    ax  = Axis(fig1[1,1], xlabel="y", ylabel="p_y")
    scatter!(ax, b1)
    plot_phasespaceof!(ax, e, param, data, fig_size=(1200,1200))    
    # --- save figure ---
    save_screen = GLMakie.Screen(; visible=false)
    display(save_screen, fig1)
    save(outfig, fig1, px_per_unit=2)
    close(save_screen)
    display(GLMakie.Screen(), fig1)
end