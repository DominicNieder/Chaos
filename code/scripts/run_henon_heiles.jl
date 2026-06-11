using GLMakie 
include("../models/henon_heiles.jl")
include("../analysis/poincare.jl")


r = range(-20, 20, length = 40)
data2d = [potential(x,y) for x in r, y in r]

f = Figure(size = (700, 400))


a2 = Axis3(f[1, 1])
contour3d!(a2, -10..10, -10..10, data2d, linewidth = 3, levels = 10)
f