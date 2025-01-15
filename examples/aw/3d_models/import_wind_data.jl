using Printf, CairoMakie, NCDatasets, Interpolations, GLMakie, Impute, Dates, Statistics
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

lon_min, lon_max = -50.375, -19.625
lon_nodes = Vector(lon_min:0.25:lon_max)
lat_min, lat_max = 24.625, 55.375
lat_nodes = Vector(lat_min:0.25:lat_max)
t_nodes = 0:1:12

ds = NCDataset("./data/copernicus/wind/new.nc");
lat_min_idx, lat_max_idx = findfirst(ds["latitude"] .== lat_min), findfirst(ds["latitude"] .== lat_max)
lon_min_idx, lon_max_idx = findfirst(ds["longitude"] .== lon_min), findfirst(ds["longitude"] .== lon_max)

data_by_month = [[],[],[],[],[],[],[],[],[],[],[],[]];
wind_dat = convert(Array{Union{Missing, Float64}}, ds["eastward_stress"][lon_min_idx:lon_max_idx, lat_min_idx:lat_max_idx, :])
for t = 1:(size(wind_dat)[end])
    wind_dat[:,:,t] = Impute.interp(wind_dat[:,:,t]) |> Impute.locf() |> Impute.nocb();
    push!(data_by_month[Dates.month(ds["time"][t])], wind_dat[:,:,t]);
end

mean_data_by_month = []
for t = 1:12
    push!(mean_data_by_month, mean(data_by_month[t]));
end
push!(mean_data_by_month, mean_data_by_month[1]);
wind_dat = stack(mean_data_by_month);

wind_itp = interpolate((lon_nodes, lat_nodes, t_nodes), wind_dat, Gridded(Linear()));
τ_surface(x, y, t) = wind_itp(x, y, mod(t, year) / month);

# Create an animation on a finer grid
x = lon_min:0.1:lon_max
y = lat_min:0.1:lat_max
x_grid = repeat(x, 1, length(y))
y_grid = repeat(y, 1, length(x))'
t_vals = (0.1:0.1:11) * month
N = length(t_vals)
matrix = τ_surface.(x_grid, y_grid, 0)

fig, ax, hm = heatmap(x, y, matrix, colormap = :balance)
record(fig, "wind_animation.mp4", 1:1:N) do i
    hm[3] = τ_surface.(x_grid, y_grid, t_vals[i]) # update data
    autolimits!(ax) # update limits
end

# Also create an animation for -40 degrees
time = Observable(0.0)
xs = @lift(τ_surface.(-40, y, ($time)*month))
fig = lines(xs, y, axis = (title = @lift("t = $(round($time, digits = 1)) months"), limits = ((-0.5, 0.5), nothing)))
framerate = 10
timestamps = range(0, 11, step=1/framerate)
record(fig, "wind_animation2.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end
