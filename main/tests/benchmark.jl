using BenchmarkTools

# Load simulation code
include("../main.jl")

println("\n" * "="^60)
println("Benchmarking Elastic Wave Solver")
println("="^60)

# Setup
mx, my = 100, 100
λ, λ_d, μ, μ_d, D, D_v1, D_v2 = get_populated_material(Params.λ_tag, Params.μ_tag, Params.D_tag, mx, my)

x_l, x_r = -1.0, 1.0
y_l, y_r = -1.0, 1.0
hx = (x_r - x_l) / (mx - 1)
hy = (y_r - y_l) / (my - 1)
x_values = LinRange(x_l, x_r, mx)
y_values = LinRange(y_l, y_r, my)

grid = Grid(x_l, x_r, y_l, y_r, mx, my, hx, hy, x_values, y_values)
stencils = compute_stencils(mx, my, hx, hy, 4)
material = Material(D, D_v1, D_v2, λ, λ_d, μ, μ_d, Params.ω)
bcs = BoundaryConditions(:traction_free, :traction_free, :traction_free, :traction_free)

v1 = rand(mx*my)
v2 = rand(mx*my)
u1 = rand(mx*my)
u2 = rand(mx*my)

t = 0.0
dt = 0.001

println("\nWarming up...")
RHS(v1, v2, u1, u2, t, grid, stencils, material, bcs)
RK4(v1, v2, u1, u2, t, dt, grid, stencils, material, bcs)

println("\n📊 Benchmarking RHS function:")
println("-" * "="^59)
b1 = @benchmark RHS($v1, $v2, $u1, $u2, $t, $grid, $stencils, $material, $bcs)
display(b1)

println("\n\n📊 Benchmarking RK4 timestep:")
println("-" * "="^59)
b2 = @benchmark RK4($v1, $v2, $u1, $u2, $t, $dt, $grid, $stencils, $material, $bcs)
display(b2)

println("\n\n✅ Summary:")
println("-" * "="^59)
println("RHS median time:  ", BenchmarkTools.prettytime(median(b1.times)))
println("RHS allocations:  ", b1.allocs, " (", BenchmarkTools.prettymemory(b1.memory), ")")
println()
println("RK4 median time:  ", BenchmarkTools.prettytime(median(b2.times)))
println("RK4 allocations:  ", b2.allocs, " (", BenchmarkTools.prettymemory(b2.memory), ")")
println("="^60)
