using CairoMakie, NCDatasets, Interpolations, Statistics
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

# Surface photosynthetically active radiation (PAR⁰)
@inline PAR⁰(t) = 90 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2;
@inline PAR(z, t) = PAR⁰(t) * exp(0.12 * z);

## Load in the data, interpolate
months_lst = ["january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december"];
itp_times = []
depth_lst = []
mean_temp_lst = []

for (j_month, m) = enumerate(months_lst)
    ds = NCDataset("../data/nc/$m.nc");
    depth = convert(Matrix{Float64}, ds["Depth"])
    temp = ds["Temperature"]
    push!(itp_times, (j_month - 1) * month)
    push!(depth_lst, depth[:, 1])
    push!(mean_temp_lst, mean(temp, dims=2))
end

# Have to add on another set of january to ensure correct looping
push!(itp_times, 12 * month)
push!(mean_temp_lst, mean_temp_lst[1])

temp_arr = reduce(hcat, mean_temp_lst)
depths = depth_lst[1]

# Create a gridded linear interpolation object
temp_itp = interpolate((depths, itp_times), temp_arr, Gridded(Linear()))
temp_function(z, t) = temp_itp(-z, mod(t, year))

t = (0:1:365) * days
z = (0:-0.1:-300) * meters
z_grid = repeat(z, 1, length(t));      # Shape z into a matrix
t_grid = repeat(t', length(z), 1);     # Shape t into a matrix

# Create the figure
fig = Figure(size = (2400, 1000), fontsize = 50)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, t[end] / days), (z[end], 0)))

axTemp = Axis(fig[1, 1]; title = "Temperature (ᵒC)", axis_kwargs...)
hmTemp = heatmap!(axTemp, t / days, z, temp_function.(z_grid, t_grid)', colormap = :balance)
Colorbar(fig[1, 2], hmTemp)

axPAR = Axis(fig[1, 3]; title = "PAR (W / m²)", xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, t[end] / days), (-100, 0)))
hmPAR = heatmap!(axPAR, t / days, z,  PAR.(z_grid, t_grid)', colormap = :balance)
Colorbar(fig[1, 4], hmPAR)

save("aux_plot.png", fig);
