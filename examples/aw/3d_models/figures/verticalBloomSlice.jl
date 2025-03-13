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

tI = 39:41;
tIdx = 40;
yIdx = 180;
xlabel = "Longitude [ᵒE]";
ylabel = "Latitude [ᵒN]";
zlabel = "Depth [m]";
tracerName = "Z";
tracerMax = maximum(abs, ds[tracerName][:,:,:,tI]);
dayNum = (tIdx-1)*10
tracerDat = ds[tracerName][:,:,:,tIdx];
fig = Figure(size = (800, 500));
ax1 = Axis(fig[1,1:2]; title = "Surface $tracerName [mmol N / m³] at $(dayNum) days", xlabel = xlabel, ylabel = ylabel);
hm1 = heatmap!(ax1, xC, yC, tracerDat[:,:,end], colormap = :balance, colorrange = (0, tracerMax));
hlines!(yF[yIdx], linestyle = :dash, color = :black);
ax2 = Axis(fig[1,3]; title = "Vertical profile at $(yF[yIdx])ᵒN", xlabel = xlabel, ylabel = zlabel, limits = (nothing, (-50, 0)));
hm2 = heatmap!(ax2, xC, zC, tracerDat[:,yIdx,:], colormap = :balance, colorrange = (0, tracerMax));
Colorbar(fig[1,4], hm2);
save("./examples/aw/3d_models/figures/verticalBloomSlice$(tracerName)$(dayNum).png", fig);


