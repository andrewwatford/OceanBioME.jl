using CairoMakie, NCDatasets, Interpolations, Statistics
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

## Load in the data, interpolate
months_lst = ["january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december"];
itp_times = []
depth_lst = []
mean_N_lst = []

for (j_month, m) = enumerate(months_lst)
    ds = NCDataset("../data/nc/$m.nc");
    depth = convert(Matrix{Float64}, ds["Depth"])
    N = ds["Nitrate"]
    push!(itp_times, (j_month - 1) * month)
    push!(depth_lst, depth[:, 1])
    push!(mean_N_lst, mean(N, dims=2))
end

# Have to add on another set of january to ensure correct looping
push!(itp_times, 12 * month)
push!(mean_N_lst, mean_N_lst[1])

N_arr = reduce(hcat, mean_N_lst)
depths = depth_lst[1]

# Create a gridded linear interpolation object
N_itp = interpolate((depths, itp_times), N_arr, Gridded(Linear()))
N_function(z, t) = N_itp(-z, mod(t, year))

t = (0:1:365) * days
z = (0:-1:-300) * meters
z_grid = repeat(z, 1, length(t));      # Shape z into a matrix
t_grid = repeat(t', length(z), 1);     # Shape t into a matrix

# Create the figure
fig = Figure(size = (1200, 1000), fontsize = 50)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, t[end] / days), (z[end], 0)))

axN = Axis(fig[1, 1]; title = "Nitrate Concentration (mmol N / m³)", axis_kwargs...)
hmN = heatmap!(axN, t / days, z, N_function.(z_grid, t_grid)', colormap = :balance)
Colorbar(fig[1, 2], hmN)

save("NWOA_plot.png", fig);
