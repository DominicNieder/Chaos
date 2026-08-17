import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using DynamicalSystems, OrdinaryDiffEq, LinearAlgebra, GLMakie, Random, JSON3, JLD2,
      NonlinearSolve, ADTypes, DataFrames, Dates, ProgressMeter, Printf

include("../styles/makie_theme.jl")
set_theme!(QUARTO_THEME)
include("../models/henon_heiles.jl")
using .HenonHeiles

const CONFIG_DIR = joinpath(@__DIR__, "../sim_config/henon_heiles.json")
const DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/simulation/simn-y256-py0/")
const SAVE_DATA_DIR   = joinpath(@__DIR__, "../../data/henon-heiles/periodic-orbits/")
const FIG_DIR    = joinpath(@__DIR__, "../../figures/henon-heiles/periodic-orbits/")
data_file        = joinpath(DATA_DIR, "E0.1127-T10000.0-py0.0-n256.jld2")

configurations = JSON3.read(read(CONFIG_DIR, String))
cfg            = configurations.explore
param          = [Float64(cfg.a.value), Float64(cfg.m.value), Float64(cfg.w.value)]
dt             = Float64(cfg.dt.value)
x0             = Float64(cfg.x0.value)

# --- background section ---
data = load(data_file, "results")
y_all, py_all = Float32[], Float32[]
for d in data
    append!(y_all,  d.sec_y)
    append!(py_all, d.sec_py)
end

const EPS_OFF = 1e-9
# precompute once, outside the loop
const RGRID  = range(-1.0, 1.0, length = 240)
const EPOT   = [HenonHeiles.potential(x, y, param) for x in RGRID, y in RGRID]
const LEVELS = collect(logrange(5.0*0.009, 6.9*0.089, 7))


### describing orbits
function kind_index(τ)
    abs(τ) < 2      && return 1   # elliptic — stable
    τ ≥  2          && return 2   # hyperbolic — unstable
    τ ≤ -2          && return 3   # inverse hyperbolic — unstable, reflection
    return 4                      # parabolic / marginal
end

const KIND_LS    = [:solid, :dash, :dashdot, :dot]
const KIND_MS    = [:circle, :xcross, :diamond, :utriangle]
const KIND_LABEL = ["elliptic", "hyperbolic", "inverse hyperbolic", "parabolic"]

f = joinpath(SAVE_DATA_DIR, "orbits_E0.01-0.1127_20260816-0742.jld2")   
df, timing, Es = load(f, "df", "timing", "Es")