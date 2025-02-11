using Printf, CairoMakie, NCDatasets, Interpolations, Impute, Statistics
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

nitrate_grids = [];
depth_nodes = [];

# Extract the surface data
for (szn, num) in szn_dict
    nitrate_ds = NCDataset("./data/woa/nitrate/woa23_all_n1$(num)_01.nc");
    lat_min_idx, lat_max_idx = findfirst(nitrate_ds["lat"] .== lat_min), findfirst(nitrate_ds["lat"] .== lat_max)
    lon_min_idx, lon_max_idx = findfirst(nitrate_ds["lon"] .== lon_min), findfirst(nitrate_ds["lon"] .== lon_max)
    nitrate_dat = nitrate_ds["n_an"][lon_min_idx:lon_max_idx, lat_min_idx:lat_max_idx, :, 1]
    nitrate_dat_imputed = []
    for z in range(1,size(nitrate_dat, 3))
        local_nitrate_dat = nitrate_dat[:,:,z]
        local_nitrate_dat_imputed = Impute.interp(local_nitrate_dat) |> Impute.locf() |> Impute.nocb()
        push!(nitrate_dat_imputed, local_nitrate_dat_imputed)
    end
    nitrate_dat_imputed = stack(nitrate_dat_imputed)
    zonally_averaged_nitrate_dat = mean(nitrate_dat_imputed, dims = 1)[1, :, :]
    push!(nitrate_grids, zonally_averaged_nitrate_dat)
    push!(depth_nodes, deepcopy(nitrate_ds["depth"][:]))
    close(nitrate_ds)
end

# To get interpolations looping
push!(nitrate_grids, nitrate_grids[1]);
nitrate_dat = stack(nitrate_grids);

# Create interpolator for surface data
depth_nodes = depth_nodes[1]
nitrate_itp = interpolate((lat_nodes, depth_nodes, t_nodes), nitrate_dat, Gridded(Linear()));

# Export to file
output_nitrate_ds = NCDataset("./data/processed/nitrate.nc", "c")
defDim(output_nitrate_ds,"lat",size(nitrate_dat,1));
defDim(output_nitrate_ds,"depth",size(nitrate_dat,2))
defDim(output_nitrate_ds,"month",size(nitrate_dat,3));
lat = defVar(output_nitrate_ds,"lat",Float64,("lat",))
depth = defVar(output_nitrate_ds,"depth",Float64,("depth",))
m = defVar(output_nitrate_ds,"month", Float64, ("month",))
nitrate = defVar(output_nitrate_ds,"nitrate",Float64, ("lat", "depth","month"))
lat[:] = lat_nodes
depth[:] = depth_nodes
m[:] = t_nodes
nitrate[:,:,:] = nitrate_dat
nitrate.attrib["units"] = "mmolN/kg"
close(output_nitrate_ds)