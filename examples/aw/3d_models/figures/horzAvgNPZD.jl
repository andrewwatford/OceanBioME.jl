using Oceananigans
using NCDatasets
using CairoMakie
using CUDA

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

ds = NCDataset("./data/runs/NPZD3_10y.nc");

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

## Predefine computations
F = Field{Center, Center, Center}(grid);
horzAvgF = Average(F, dims = (1,2))

## Perform computations
labels = ["N", "P", "Z", "D"]
timeseries = []
for (j, label) in enumerate(labels)
	dat = CuArray(ds[label][:,:,:,:]);
	horzAvgTimeSeries = Array{Float64}(undef, Nt, Nz);
	CUDA.allowscalar() do
		for ti in 1:Nt
			set!(F, dat[:,:,:,ti]);
			horzAvgTimeSeries[ti, :] = compute!(Field(horzAvgF))[1,1,1:Nz];
		end
	end
	push!(timeseries, horzAvgTimeSeries)
end

## Create the figures
for (j, label) in enumerate(labels)
	fig = Figure()
	hmAx = Axis(fig[1,1]; title = "Horizontally averaged $(label) [mmolN/m³] over 1 year period", xlabel = "t [mo]", ylabel = "z [m]", limits = ((0, 12), (-100, 0)))
	hm = heatmap!(hmAx, t[1:37] / month, zC, timeseries[j][1:37,:], colormap = :balance)
	Colorbar(fig[1,2], hm);
	save("./figures/horzAvgWeakYear$(label).png", fig)
end
for (j, label) in enumerate(labels)
	fig = Figure()
	hmAx = Axis(fig[1,1]; title = "Horizontally averaged $(label) [mmolN/m³] over 1 year period", xlabel = "t [mo]", ylabel = "z [m]", limits = ((12, 24), (-100, 0)))
	hm = heatmap!(hmAx, t[37:74] / month, zC, timeseries[j][37:74,:], colormap = :balance)
	Colorbar(fig[1,2], hm);
	save("./figures/horzAvgStrongYear$(label).png", fig)
end
using Statistics
for (j, label) in enumerate(labels)
	fig = Figure()
	ax = Axis(fig[1,1]; title = "Average $(label) [mmolN/m³] over 10 year period", xlabel = "t [y]", ylabel = "$(label) [mmolN/m³]")
	plot = lines!(t / year, mean(timeseries[j], dims=2)[:,1])
	save("./figures/horzAvgTS$(label).png", fig)
end
close(ds);
