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

# fig = Figure(size = (1000,1000))
# axis_kwargs = (xlabel = "Longitude", ylabel = "Latitude");
# x = -50:0.1:-20;
# y = 25:0.1:40;
# xgrid = repeat(x, 1, length(y));
# ygrid = repeat(y, 1, length(x))';

# for (m, data) = enumerate(mean_data_by_month[1:12])
#     i = 1 + div(m-1, 3);
#     j = 1 + mod(m-1, 3);
#     ax = Axis(fig[i, j]; title = "Month $(m)", axis_kwargs...)
#     hm = heatmap!(x, y, wind_itp.(xgrid, ygrid, m - 1), colormap = :balance, colorrange = (-0.2, 0.2))
#     if m == 12
#         Colorbar(fig[1:4, 4], hm)
#     end
# end

# save("fig.png", fig);


