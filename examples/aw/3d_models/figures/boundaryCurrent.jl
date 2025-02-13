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

### Step 1: Identify the boundary current and its extent ###
fig = Figure(); # Play with figure size?
xlabel = "Longitude [ᵒE]";
ylabel = "Latitude [ᵒN]";
zlabel = "Depth [m]";
xlim = (-45, -43);
zlim = (-750, 0);
colorrange = (-0.05, 0.25)
## (a) Surface view
surfaceV = ds["v"][:,:,end,end];
surfaceAx = Axis(fig[1,1]; title = "Surface v [m/s] after spin-up", limits = (xlim, nothing), xlabel = xlabel, ylabel = ylabel);
surfaceHm = heatmap!(surfaceAx, xC, yF, surfaceV, colormap = :balance, colorrange = colorrange);
yIdx = 43
hlines!(yF[yIdx], linestyle = :dash, color = :black)
## (b) Vertical profile at peak
profileV = ds["v"][:,yIdx,:,end];
profileAx = Axis(fig[1,2]; title = "Vertical profile at $(yF[yIdx])ᵒN", limits = (xlim, zlim), xlabel = xlabel, ylabel = zlabel);
profileHm = heatmap!(profileAx, xC, zC, profileV, colormap = :balance, colorrange = colorrange);
Colorbar(fig[1,3], profileHm);
save("./figures/boundaryCurrent.png", fig);
close(ds)