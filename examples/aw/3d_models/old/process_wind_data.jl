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
    imputed_data = Impute.interp(wind_dat[:,:,t]) |> Impute.locf() |> Impute.nocb();
    zonally_averaged = mean(imputed_data, dims = 1)[1, :]
    push!(data_by_month[Dates.month(ds["time"][t])], zonally_averaged);
end

mean_data_by_month = []
for t = 1:12
    push!(mean_data_by_month, mean(data_by_month[t]));
end
push!(mean_data_by_month, mean_data_by_month[1]);
wind_dat = stack(mean_data_by_month);

wind_itp = interpolate((lat_nodes, t_nodes), wind_dat, Gridded(Linear()));

# Export to file
output_ds = NCDataset("./data/processed/winds.nc", "c")
defDim(output_ds,"lat",size(wind_dat, 1));
defDim(output_ds,"month",size(wind_dat, 2))
τ_x = defVar(output_ds,"τ_x",Float64,("lat","month"))
τ_x[:,:] = wind_dat
τ_x.attrib["units"] = "N/m²"
lat = defVar(output_ds,"lat", Float64,("lat",))
lat[:] = lat_nodes
m = defVar(output_ds,"month", Float64, ("month",))
m[:] = t_nodes
close(output_ds)

# fig = Figure(size = (1000,1000))
# axis_kwargs = (xlabel = "Eastward Stress", ylabel = "Latitude", limits = ((-0.1, 0.1), nothing));
# y = 25:0.1:40;

# for (m, data) = enumerate(mean_data_by_month[1:12])
#     i = 1 + div(m-1, 3);
#     j = 1 + mod(m-1, 3);
#     ax = Axis(fig[i, j]; title = "Month $(m)", axis_kwargs...)
#     lines!(wind_itp.(y, m - 1), y)
# end

# save("zonal_avg_winds.png", fig);
