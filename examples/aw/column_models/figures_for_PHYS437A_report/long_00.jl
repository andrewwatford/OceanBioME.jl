using Oceananigans, CairoMakie
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

file = "../data_files_for_PHYS437A_report/long_00.jld2"

@inline Q₁₀(T) = 1.88 ^ (T / 10) # T in °C
@inline light_limitation(PAR, α, μ₀) = α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)
@inline nutrient_limitation(N, kₙ) = N / (kₙ + N)

kₙ = 3.4151
μ₀ = 0.7894/day
α = 0.2099/day

N = FieldTimeSeries(file, "N")
P = FieldTimeSeries(file, "P")
Z = FieldTimeSeries(file, "Z")
D = FieldTimeSeries(file, "D")
PAR = FieldTimeSeries(file, "PAR")
T = FieldTimeSeries(file, "T")
x, y, z = nodes(N)
times = N.times;
# Create a grid for z and t
z_grid = repeat(z, 1, length(times));      # Shape z into a matrix
t_grid = repeat(times', length(z), 1);     # Shape t into a matrix

fig = Figure(size = (2560, 2304), fontsize = 50)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, times[end] / days), (z[1], 0)))
# Plot the 
axN = Axis(fig[1, 1]; title = "N (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmN = heatmap!(times / days, z, N[1,1,:,:]', colormap = :balance, colorrange = (2.5, 12.5))
Colorbar(fig[1, 2], hmN)
axP = Axis(fig[2, 1]; title = "P (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmP = heatmap!(times / days, z, P[1,1,:,:]', colormap = :balance, colorrange = (0, 0.4))
Colorbar(fig[2, 2], hmP)
axZ = Axis(fig[3, 1]; title = "Z (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmZ = heatmap!(times / days, z, Z[1,1,:,:]', colormap = :balance, colorrange = (0, 0.4))
Colorbar(fig[3, 2], hmZ)
axD = Axis(fig[4, 1]; title = "D (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmD = heatmap!(times / days, z, D[1,1,:,:]', colormap = :balance, colorrange = (0, 0.7))
Colorbar(fig[4, 2], hmD) 

save("long_00.png", fig);