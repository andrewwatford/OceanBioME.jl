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

## Predefine computations
T = Field{Center, Center, Center}(grid);
zonAvgT = Average(T, dims = (1,));
merAvgT = Average(T, dims = (2,));
horzAvgT = Average(T, dims = (1,2));

## Compute the horizontally averaged profile for each time
datT = CuArray(ds["T"][:,:,:,:]);
horzAvgTimeSeries = Array{Float64}(undef, Nt, Nz);
CUDA.allowscalar() do
	for ti in 1:Nt
		set!(T, datT[:,:,:,ti]);
		horzAvgTimeSeries[ti, :] = compute!(Field(horzAvgT))[1,1,1:Nz];
	end
end

## Create the heatmap
fig = Figure(size=(1000, 600));
hmAx = Axis(fig[1,1:2], title = "Horizontally averaged temperature [ᵒC] over spin-up period", xlabel = "Time [y]", ylabel = "z [m]", limits = (nothing, (-1000, 0)));
hm = heatmap!(hmAx, t / year, zC, horzAvgTimeSeries, colormap = :balance);
Colorbar(fig[1,3], hm);
## Create profiles for the last three data points
zAx = Axis(fig[1,4], title = "Vertical profiles at end of spin-up", xlabel = "T [ᵒC]", ylabel = "z [m]", limits = ((18, 25), (-100, 0)));
timeLabels = ["Jan", "May", "Sep"]
timeIdx = [3, 2, 1]
for (j, label) in enumerate(timeLabels)
	idx = timeIdx[j];
	plt = lines!(zAx, horzAvgTimeSeries[(Nt-idx),:], zC[1:Nz], label = label)
end
axislegend()
save("./figures/horzAvgTemp.png", fig);

close(ds);