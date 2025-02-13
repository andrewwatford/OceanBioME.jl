using Oceananigans
using NCDatasets
using CairoMakie
using CUDA

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

ds = NCDataset("./data/runs/spinup0_100y.nc");

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

### Computations to get vector field to play nice ###
surfaceU = ds["u"][:,:,end,end];
surfaceU = (surfaceU[1:end-1, :] + surfaceU[2:end, :]) / 2;
surfaceV = ds["v"][:,:,end,end];
surfaceV = (surfaceV[:, 1:end-1] + surfaceV[:, 2:end]) / 2;

### Plot stuff ###
fig = Figure(); # Play with figure size?
xlabel = "Longitude [ᵒE]";
ylabel = "Latitude [ᵒN]";
surfaceAx = Axis(fig[1,1]; title = "Surface velocity [m/s] after spin-up", xlabel = xlabel, ylabel = ylabel);
arrows!(xC[1:10:end], yC[1:10:end], surfaceU[1:10:end,1:10:end], surfaceV[1:10:end,1:10:end], arrowsize = 8, lengthscale = 8);
save("./figures/surfaceFlow.png", fig);
close(ds)