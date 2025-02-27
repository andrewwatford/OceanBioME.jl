using Oceananigans
using NCDatasets
using CairoMakie
using GLMakie
using CUDA

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

ds = NCDataset("./data/runs/NPZD2_20y.nc");

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
## Create the figures
for (j, label) in enumerate(labels)
	fig = Figure()
	hmAx = Axis(fig[1,1]; title = "Surface $(label) [mmolN/m³]", xlabel = xlabel, ylabel = ylabel)
	hm = heatmap!(hmAx, xC, yC, ds[label][:,:,end,end], colormap = :balance)
	Colorbar(fig[1,2], hm);
	save("./figures/surface$(label).png", fig)
end

for (j, label) in enumerate(labels)
	time = Observable(0.0)
	data = ds[label][:,:,end,:]
	maxval = maximum(abs, data)
	ax_tpl = (; title = @lift("Surface $(label) [mmolN/m³] at month $time"), xlabel = xlabel, ylabel = ylabel)
	fig, ax, hm = heatmap(xC, yC, data[:,:,1], colormap = :balance, colorrange = (0, maxval), axis = ax_tpl)
	Colorbar(fig[1, 2], hm)
	CairoMakie.record(fig, "./figures/surface$(label).mp4", 1:1:12, framerate = 2) do i
		hm[3] = data[:,:,i]
		time[] = i
		autolimits!(ax)
	end
end

close(ds);
