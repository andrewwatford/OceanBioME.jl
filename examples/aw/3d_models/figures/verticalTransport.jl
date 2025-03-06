using Oceananigans
using NCDatasets
using CairoMakie
using CUDA
using Statistics

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

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
wField = Field{Center, Center, Face}(grid);
NField = Field{Center, Center, Center}(grid);
fluxField = Field(wField * NField);

## Compute
w = CuArray(ds["w"][:,:,:,:]);
N = CuArray(ds["N"][:,:,:,:]);
w = mean(w, dims=4);
N = mean(N, dims=4);
set!(wField, w);
set!(NField, N);
flux = compute!(fluxField);
CUDA.@allowscalar fluxUpper = Array(flux[1:Nλ,1:Nφ,end] - flux[1:Nλ,1:Nφ,172]);

### Figures
fig = Figure();
ax = Axis(fig[1,1]; title = "Vertical nutrient flux into top 100m [mmol N / m² s]", xlabel = "Longitude [ᵒE]", ylabel = "Latitude [ᵒN]");
co = contourf!(ax, xC, yC, fluxUpper);
Colorbar(fig[1,2], co);
save("./examples/aw/3d_models/figures/verticalTransport.png", fig)