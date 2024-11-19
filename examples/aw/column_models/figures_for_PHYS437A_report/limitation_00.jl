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
PAR = FieldTimeSeries(file, "PAR")
T = FieldTimeSeries(file, "T")
LL = light_limitation.(PAR, α, μ₀ .* Q₁₀.(T))[1,1,:,:]'
NL =  nutrient_limitation.(N, kₙ)[1,1,:,:]'
x, y, z = nodes(N)
times = N.times;
# Create a grid for z and t
z_grid = repeat(z, 1, length(times));      # Shape z into a matrix
t_grid = repeat(times', length(z), 1);     # Shape t into a matrix

# Create the figure
fig = Figure(size = (2400, 1000), fontsize = 50)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((365.0, 730.0), (z[1], 0)))

axLL = Axis(fig[1, 1]; title = "Light Limitation", xlabel = "Time (days)", ylabel = "z (m)", limits = ((365.0, 730.0), (-100, 0)))
hmLL = heatmap!(axLL, times / days, z, LL, colormap = :balance, colorrange = (0.0, 1.0))
Colorbar(fig[1, 2], hmLL)

axNL = Axis(fig[1, 3]; title = "Nutrient Limitation", axis_kwargs...)
hmNL = heatmap!(axNL, times / days, z, NL, colormap = :balance, colorrange = (0.0, 1.0))
Colorbar(fig[1, 4], hmNL)

# axL = Axis(fig[1, 5]; title = "Total Limitation", axis_kwargs...)
# hmL = heatmap!(axL, times / days, z, LL .* NL, colormap = :balance)
# Colorbar(fig[1, 6], hmL)

save("limitation_00.png", fig);