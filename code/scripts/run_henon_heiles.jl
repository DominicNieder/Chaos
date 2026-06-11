using GLMakie 
include("../models/henon_heiles.jl")
include("../analysis/poincare.jl")


r = range(-0.4, 0.4, length = 120)
data2d = [potential(x,y) for x in r, y in r]

f = Figure(size = (600, 400))


a2 = Axis3(f[1, 1])
contour3d!(a2, -0.4..0.4, -0.4..0.4, data2d, linewidth = 3, levels = 220)
f