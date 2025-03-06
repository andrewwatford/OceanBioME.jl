using Oceananigans
using NCDatasets
using CairoMakie
using GLMakie
using CUDA

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

## Labels
labels = ["N", "P", "Z", "D"]
xlabel = "Longitude [ᵒE]";
ylabel = "Latitude [ᵒN]";
## Create the figure
index = Observable(1);
data = ds["P"][:,:,end,:];
maxval = maximum(abs, data);
ax_tpl = (; title = @lift("Surface P [mmolN/m³] at $(($(index)-1)*10) days"), xlabel = xlabel, ylabel = ylabel);
fig, ax, hm = heatmap(xC, yC, data[:,:,1], colormap = :balance, colorrange = (0, maxval), axis = ax_tpl);
Colorbar(fig[1, 2], hm);
video = VideoStream(fig, format = "mp4", framerate = 6);
for i = 1:1:366
	index[] = i;
	hm[3] = data[:,:,i];
	hlines!(yF[189], linestyle = :dash, color = :black);
	recordframe!(video);
	yield();
end
save("./examples/aw/3d_models/figures/P_animation.mp4", video);
close(ds);
