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

ds = NCDataset("./examples/aw/3d_models/data/runs/NPZD3_10y.nc");

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
vField = Field{Center, Face, Center}(grid);
vBtField = Field(Average(vField, dims=3));
vEkmanField = Field(Average(Field(vField - vBtField), dims=1));

## Compute
vDat = CuArray(ds["v"][:,:,:,:]);
vAvg = mean(vDat, dims=4);
set!(vField, vAvg);
CUDA.@allowscalar vEkman = Array(compute!(vEkmanField));

## Figure
fig = Figure();
ax = Axis(fig[1,1]; title = "Average ageostrophic meridional velocity [m/s]", xlabel = "Latitude [ᵒN]", ylabel = "z [m]", limits = (nothing, (-100, 0)));
hm = heatmap!(ax, yF, zC, vEkman[1,1:size(yF,1), 1:Nz], colormap = :balance);
Colorbar(fig[1,2], hm);

save("./examples/aw/3d_models/figures/ekmanTransport.png", fig)
