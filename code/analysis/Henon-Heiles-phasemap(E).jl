import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using GLMakie, JLD2, Glob, JSON3
include("../styles/makie_theme.jl")
include("../models/henon_heiles.jl")
using .HenonHeiles

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

DATA_DIR        = joinpath(@__DIR__, "../../data/henon-heiles/simulation")
FIG_DIR         = joinpath(@__DIR__, "../../figures/henon-heiles/simulation")

FIG_ORIENTATION  = joinpath(FIG_DIR, "orientation.json")
DATA_ORIENTATION = joinpath(DATA_DIR, "orientation.json")

EXPL_DATA_DIR = joinpath(@__DIR__, "../../data/henon-heiles/explore")

json = JSON3.read(read(DATA_ORIENTATION, String))

run = json[:"(E,y)100x256_p00.0/"]

file = map(i -> String(i), run.file)
energies = map(i -> Float64(i),run.energy)

full_screen_size= (1920 , 1200)
save_size       = (1200, 1200)
files = sort(glob("E*-T*-py0.0-n*.jld2", DATA_DIR))
energies = [parse(Float64, match(r"E([\d.]+)-", basename(f))[1]) for f in files]


# --- preload everything once, flattened per frame ---
all_y  = Vector{Vector{Float32}}(undef, length(files))
all_py = Vector{Vector{Float32}}(undef, length(files))
all_c  = Vector{Vector{Float32}}(undef, length(files))
bound, ylim, pylim  = section_set_plot_lim(energies[1],param,n)

for (fi, f) in enumerate(files)
    data = load(f, "results")
    y_all, py_all, c_all = Float32[], Float32[], Float32[]
    for (i, d) in enumerate(data)
        append!(y_all, d.sec_y)
        append!(py_all, d.sec_py)
        append!(c_all, fill(Float32(i), length(d.sec_y)))
    end
    all_y[fi], all_py[fi], all_c[fi] = y_all, py_all, c_all
end

with_theme(QUARTO_THEME) do
    fig = Figure(size=full_screen_size)
    ax = Axis(fig[1, 1][2:10,1:10], xlabel="y", ylabel="p_y", aspect=nothing)    
    sl = Slider(fig[1, 1][1,1:10], range=1:length(files), startvalue=1)
    Label(fig[1, 1][2,1:10], @lift("E = $(energies[$(sl.value)])"))

    ys  = Observable(all_y[20])
    pys = Observable(all_py[20])
    cs  = Observable(all_c[20])
    boundary = Observable(bound)

    scatter!(ax, ys, pys, color=cs, colormap=:viridis, markersize=1.5)
    scatter!(ax, boundary)
    on(sl.value) do idx
        ys[]  = all_y[idx]
        pys[] = all_py[idx]
        cs[]  = all_c[idx]
        boundary[], ylim, pylim  = section_set_plot_lim(energies[idx],param,n)
        xlim!(ax, ylim...)
        ylim(ax, pylim...)
    end

    display(GLMakie.Screen(), fig)
end