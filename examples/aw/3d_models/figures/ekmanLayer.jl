using Oceananigans
using NCDatasets
using CairoMakie
using CUDA
using Statistics

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

const ρ₀ = 1035 # [kg m⁻³] reference density

ds = NCDataset("./data/runs/NPZD3_10y.nc");
# wind_ds = NCDataset("./data/processed/winds.nc")

## Create the grid for averaging purposes
xF, xC = ds["xF"][:], ds["xC"][:];
yF, yC = ds["yF"][:], ds["yC"][:];
zF, zC = ds["zF"][:], ds["zC"][:];
t = ds["time"][:];
Nλ, Nφ, Nz, Nt = size(xC, 1), size(yC, 1), size(zC, 1), size(t, 1);

grid = LatitudeLongitudeGrid(GPU(), Float32;
                           size = (Nλ, Nφ, Nz),
                           halo = (6, 6, 6),
                      longitude = xF,
                       latitude = yF,
                              z = zF,
                       topology = (Bounded, Bounded, Bounded));

## Predefine the relevant computations
V = Field{Center, Face, Center}(grid);
mᵥ = CumulativeIntegral(V, dims=3);

## Compute
v = CuArray(ds["v"][:,:,:,:]);
vAvg = mean(v, dims=4);
set!(V, vAvg);
mᵥ_field = compute!(Field(mᵥ));
mᵥAvg_field = compute!(Field(Average(mᵥ_field, dims=1)));

## Figure
fig = Figure();
ax = Axis(fig[1,1]; title = "Cumulative meridional velocity [m²/s]", xlabel = "Latitude [ᵒN]", ylabel = "z [m]", limits = (nothing, (-1000, 0)));
hm = heatmap!(ax, yF, zC, -mᵥAvg_field, colormap = :balance);
Colorbar(fig[1,2], hm);

save("./figures/ekmanTransport.png", fig)
