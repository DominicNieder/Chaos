# -*- coding: utf-8 -*-
# Pluto.jl notebook

### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ deps
begin
    using PlutoUI, OrdinaryDiffEq, CairoMakie
    include("../models/henon_heiles.jl")
    include("../analysis/poincare.jl")
end

# ╔═╡ slider
@bind E Slider(0.05:0.01:0.5, default=0.125, show_value=true)

# ╔═╡ solve
begin
    x0, y0 = 0.0, 0.0
    px0 = sqrt(max(0.0, 2E - y0^2 - x0^2))
    sol = solve_trajectory([x0, y0, px0, 0.0], tspan=(0.0, 500.0))
end

# ╔═╡ poincare plot
begin
    xs, pxs = poincare_section(sol)
    fig, ax, _ = scatter(xs, pxs, markersize=3,
        axis=(; title="Poincaré section  E = $E", xlabel="x", ylabel="pₓ"))
    fig
end
