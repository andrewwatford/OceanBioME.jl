using Printf, CairoMakie, NCDatasets, Interpolations, GLMakie, Impute, Statistics
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

szn_dict = (
            "winter"=>"3",
            "spring"=>"4",
            "summer"=>"5",
            "autumn"=>"6",
);

lon_min, lon_max = -50.5, -19.5
lon_nodes = Vector(lon_min:1:lon_max)
lat_min, lat_max = 24.5, 55.5
lat_nodes = Vector(lat_min:1:lat_max)
t_nodes = [0, 3, 6, 9, 12]

temp_grids = [];
sal_grids = [];

# Extract the surface data
for (szn, num) in szn_dict
    temp_ds = NCDataset("./data/woa/temperature/woa23_decav_t1$(num)_01.nc");
    lat_min_idx, lat_max_idx = findfirst(temp_ds["lat"] .== lat_min), findfirst(temp_ds["lat"] .== lat_max)
    lon_min_idx, lon_max_idx = findfirst(temp_ds["lon"] .== lon_min), findfirst(temp_ds["lon"] .== lon_max)
    temp_dat = convert(Matrix{Union{Missing, Float64}}, temp_ds["t_an"][lon_min_idx:lon_max_idx, lat_min_idx:lat_max_idx, 1, 1])
    temp_dat = Impute.interp(temp_dat) |> Impute.locf() |> Impute.nocb()
    zonally_averaged_temp_dat = mean(temp_dat, dims = 1)[1, :, :]
    push!(temp_grids, zonally_averaged_temp_dat)
    close(temp_ds)

    sal_ds = NCDataset("./data/woa/salinity/woa23_decav_s1$(num)_01.nc");
    lat_min_idx, lat_max_idx = findfirst(sal_ds["lat"] .== lat_min), findfirst(sal_ds["lat"] .== lat_max)
    lon_min_idx, lon_max_idx = findfirst(sal_ds["lon"] .== lon_min), findfirst(sal_ds["lon"] .== lon_max)
    sal_dat = convert(Matrix{Union{Missing, Float64}}, sal_ds["s_an"][lon_min_idx:lon_max_idx, lat_min_idx:lat_max_idx, 1, 1])
    sal_dat = Impute.interp(sal_dat) |> Impute.locf() |> Impute.nocb()
    zonally_averaged_sal_dat = mean(sal_dat, dims = 1)[1, :, :]
    push!(sal_grids, zonally_averaged_sal_dat)
    close(sal_ds)
end

# To get interpolations looping
push!(temp_grids, temp_grids[1]);
push!(sal_grids, sal_grids[1]);

temp_dat = stack(temp_grids)[:,1,:];
sal_dat = stack(sal_grids)[:,1,:];

# Create interpolator for surface data
temp_itp = interpolate((lat_nodes, t_nodes), temp_dat, Gridded(Linear()));
sal_itp = interpolate((lat_nodes, t_nodes), sal_dat, Gridded(Linear()));

# Export to file
output_temp_ds = NCDataset("./data/processed/sst.nc", "c")
defDim(output_temp_ds,"lat",size(temp_dat,1));
defDim(output_temp_ds,"month",size(temp_dat,2));
lat = defVar(output_temp_ds,"lat",Float64,("lat",))
m = defVar(output_temp_ds,"month", Float64, ("month",))
temp = defVar(output_temp_ds,"temp",Float64, ("lat", "month"))
lat[:] = lat_nodes
m[:] = t_nodes
temp[:,:] = temp_dat
temp.attrib["units"] = "ᵒC"
close(output_temp_ds)

output_sal_ds = NCDataset("./data/processed/salinity.nc", "c")
defDim(output_sal_ds,"lat",size(sal_dat,1));
defDim(output_sal_ds,"month",size(sal_dat,2));
lat = defVar(output_sal_ds,"lat",Float64,("lat",))
m = defVar(output_sal_ds,"month", Float64, ("month",))
sal = defVar(output_sal_ds,"salinity",Float64, ("lat", "month"))
lat[:] = lat_nodes
m[:] = t_nodes
sal[:,:] = sal_dat
sal.attrib["units"] = "practical salinity"
close(output_sal_ds)

# # Create some animations on a finer grid
# x = lon_min:0.1:lon_max
# y = lat_min:0.1:lat_max
# x_grid = repeat(x, 1, length(y))
# y_grid = repeat(y, 1, length(x))'
# t_vals = (0.1:0.1:11) * month
# N = length(t_vals)

# S_matrix = S_surface.(x_grid, y_grid, 0)
# T_matrix = T_surface.(x_grid, y_grid, 0)

# fig, ax, hm = heatmap(x, y, S_matrix, colormap = :balance, colorrange = (()))
# record(fig, "salinity_animation.mp4", 1:1:N) do i
#     hm[3] = S_surface.(x_grid, y_grid, t_vals[i]) # update data
#     autolimits!(ax) # update limits
# end

# fig, ax, hm = heatmap(x, y, T_matrix, colormap = :balance)
# record(fig, "temp_animation.mp4", 1:1:N) do i
#     hm[3] = S_surface.(x_grid, y_grid, t_vals[i]) # update data
#     autolimits!(ax) # update limits
# end