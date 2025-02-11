using Printf, CairoMakie, NCDatasets, Statistics
using Oceananigans
using Oceananigans.Units
using Printf

const year = years = 365days;
const month = months = year / 12;

const Lλ = 20
const Lφ = 15
const Lz = 4000

const λ_start = -45
const φ_start = 25

const Nλ = 200
const Nφ = 200
const Nz = 200
const halo_size = 6

const α = 2e-4  # [ᵒC⁻¹] thermal expansion
const β = 7.4e-4 # [(g/kg)⁻¹] haline contraction

### Load in and trim the datasets ###
sim_ds = NCDataset("./spinup.nc");
woa_temp_ds = NCDataset("./data/woa/temperature/woa23_decav_t13_01.nc");
woa_sal_ds = NCDataset("./data/woa/salinity/woa23_decav_s13_01.nc");

### Create the grids and predefine computations ###
# Simulation #
chebychev_spaced_z_faces(k) = -Lz + Lz * sin(π * (k - 1) / 2 / Nz);
sim_grid = LatitudeLongitudeGrid(CPU();
                           size = (Nλ, Nφ, Nz),
                           halo = (halo_size, halo_size, halo_size),
                      longitude = (λ_start,  λ_start + Lλ),
                       latitude = (φ_start,  φ_start + Lφ),
                              z = chebychev_spaced_z_faces,
                       topology = (Bounded, Bounded, Bounded))

simT = Field{Center, Center, Center}(sim_grid);
avgsimT = Field(Average(simT, dims = (1,2)));
simS = Field{Center, Center, Center}(sim_grid);
avgsimS = Field(Average(simS, dims = (1,2)));
simB = α*avgsimT - β*avgsimS;
simN² = ∂z(simB) * (hour / 2pi)^2
# WOA #
lon_min, lon_max = λ_start - 0.5, λ_start + Lλ + 0.5;
lat_min, lat_max = φ_start - 0.5, φ_start + Lφ + 0.5;
lat_min_idx, lat_max_idx = findfirst(woa_temp_ds["lat"] .== lat_min), findfirst(woa_temp_ds["lat"] .== lat_max)
lon_min_idx, lon_max_idx = findfirst(woa_temp_ds["lon"] .== lon_min), findfirst(woa_temp_ds["lon"] .== lon_max)
depth_nodes = woa_temp_ds["depth"]
woa_grid = LatitudeLongitudeGrid(CPU();
            size = (Lλ+2, Lφ+2,101),
            halo = (halo_size, halo_size, halo_size),
            longitude = (λ_start,  λ_start + Lλ),
            latitude = (φ_start,  φ_start + Lφ),
            z = depth_nodes,
            topology = (Bounded, Bounded, Bounded))
woaT = Field{Center, Center, Center}(woa_grid);
avgwoaT = Field(Average(woaT, dims = (1,2)));
woaS = Field{Center, Center, Center}(woa_grid);
avgwoaS = Field(Average(woaS, dims = (1,2)));
woaB = α * avgwoaT - β * avgwoaS
woaN² = ∂z(woaB) * (hour / 2pi)^2

### Perform computations ###
set!(simT, sim_ds["T"][:,:,:,end])
simT_data = compute!(Field(avgsimT))
set!(simS, sim_ds["S"][:,:,:,end])
simS_data = compute!(Field(avgsimS))
set!(woaT, replace!(woa_temp_ds["t_an"][lon_min_idx:lon_max_idx,lat_min_idx:lat_max_idx,1:end-1,1], missing => 0))
woaT_data = compute!(Field(avgwoaT))
set!(woaS, replace!(woa_sal_ds["s_an"][lon_min_idx:lon_max_idx,lat_min_idx:lat_max_idx,1:end-1,1], missing => 0))
woaS_data = compute!(Field(avgwoaS))

### Plot temperature on same axes ###
fig = Figure(size = (400, 600))
axis_kwargs = (xlabel = "Temperature (ᵒC)", ylabel = "z (m)")
ax = Axis(fig[1,1]; title = "Temperature profile", axis_kwargs...)
# plot simulated data #
lines!(ax, simS_data[1:size(simS_data, 3)], sim_ds["zC"])
# plot WOA data #
lines!(ax, woaS_data[1,1,1:size(woaS_data, 3)], -depth_nodes[2:end]; linestyle = :dash)

save("./figures/T_plot.png", fig)

# close(sim_ds)
close(woa_temp_ds)
close(woa_sal_ds)




## Load in the simulated data
xF = sim_ds["xF"];
yF = sim_ds["yF"];
zF = sim_ds["zF"];

# Pre-define the computations
simT = Field{Center, Center, Center}(grid);
simT = Average(simT, dims = (1,2));
simS = Field{Center, Center, Center}(grid);
simS = Average(simS, dims = (1,2));
simB = α * simT - β * simS
simN² = ∂z(simB) * (hour / 2pi)^2

### Compute the simulation data
set!(T, sim_ds["T"][:,:,:,end-3])
set!(S, sim_ds["S"][:,:,:,end-3])
datJ = compute!(Field(averageN²))

set!(T, sim_ds["T"][:,:,:,end-2])
set!(S, sim_ds["S"][:,:,:,end-2])
datM = compute!(Field(averageN²))

set!(T, sim_ds["T"][:,:,:,end-1])
set!(S, sim_ds["S"][:,:,:,end-1])
datS = compute!(Field(averageN²))

simData = [datJ, datM, datS]

########### Part II: WOA DATA and plotting ###########

# Set up the figure
fig = Figure(size = (1200, 600))
axis_kwargs = (xlabel = "N² (cph²)", ylabel = "z (m)", limits = ((-10, 10), (-2000, 0)))

lon_min, lon_max = λ_start - 0.5, λ_start + Lλ + 0.5;
lat_min, lat_max = φ_start - 0.5, φ_start + Lφ + 0.5;

szn_dict = ["3", "4", "5"]

for (axNum, num) in enumerate(szn_dict)
    temp_ds = NCDataset("./data/woa/temperature/woa23_decav_t1$(num)_01.nc");
    sal_ds = NCDataset("./data/woa/salinity/woa23_decav_s1$(num)_01.nc");
    lat_min_idx, lat_max_idx = findfirst(temp_ds["lat"] .== lat_min), findfirst(temp_ds["lat"] .== lat_max)
    lon_min_idx, lon_max_idx = findfirst(temp_ds["lon"] .== lon_min), findfirst(temp_ds["lon"] .== lon_max)
    depth_nodes = temp_ds["depth"]

    WOAgrid = LatitudeLongitudeGrid(CPU();
            size = (Lλ+2, Lφ+2,101),
            halo = (halo_size, halo_size, halo_size),
            longitude = (λ_start,  λ_start + Lλ),
            latitude = (φ_start,  φ_start + Lφ),
            z = depth_nodes,
            topology = (Bounded, Bounded, Bounded))
    
    T = Field{Center, Center, Center}(WOAgrid);
    S = Field{Center, Center, Center}(WOAgrid);
    b = α * T - β * S
    N² = ∂z(b) * (hour / 2pi)^2
    averageN² = Average(N², dims = (1,2))

    set!(T, replace!(temp_ds["t_an"][lon_min_idx:lon_max_idx,lat_min_idx:lat_max_idx,1:end-1,1], missing => 0))
    set!(S, replace!(sal_ds["s_an"][lon_min_idx:lon_max_idx,lat_min_idx:lat_max_idx,1:end-1,1], missing => 0))
    dat = compute!(Field(averageN²))

    ax = Axis(fig[1,axNum]; title = "Month $((axNum-1)*4+1)", axis_kwargs...)
    plotSim = lines!(simData[axNum][1:size(zF,1)], zF)
    plotWOA = lines!(dat[1:102], depth_nodes)
    
    close(temp_ds)
    close(sal_ds)
end

save("./figures/N_plot_4.png", fig)

# close(sim_ds)