using CairoMakie, NCDatasets, Interpolations
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

# Create the MLD interpolator
mld_ds = NCDataset("../data/nc/mixed_layer_processed.nc")
mld_dat = convert(Array{Float64}, mld_ds["depth_mean"])
push!(mld_dat, mld_dat[1]) # to make looping work
mld_itp_months = interpolate(-mld_dat, BSpline(Cubic(Line(OnGrid()))))
mld_itp(t) = mld_itp_months(mod(t / month, 12) + 1)

# Define diffusivity based on MLD
@inline κₜ(z, t) = 1e-2 * (1 + tanh((z - mld_itp(t)) / 10)) / 2 + 1e-4

t = (0:1:365) * days
z = (0:-1:-300) * meters
z_grid = repeat(z, 1, length(t));      # Shape z into a matrix
t_grid = repeat(t', length(z), 1);     # Shape t into a matrix

# Create the figure
fig = Figure(size = (2400, 1000), fontsize = 50)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, t[end] / days), (z[end], 0)))

axMLD = Axis(fig[1, 1]; title = "Mixed Layer Depth (m)", axis_kwargs...)
lines!(axMLD, t / days, mld_itp.(t), linewidth = 5, color = :black)

axDiff = Axis(fig[1, 2]; title = "Turbulent Diffusivity (m² / s)", axis_kwargs...)
hmDiff = heatmap!(axDiff, t / days, z,  κₜ.(z_grid, t_grid)', colormap = :balance, colorscale = log10)
Colorbar(fig[1, 3], hmDiff)

save("MLD_plot.png", fig);
