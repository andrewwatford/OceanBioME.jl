using Oceananigans
using NCDatasets
using CairoMakie
using CUDA

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

const α = 2e-4  # [ᵒC⁻¹] thermal expansion
const β = 7.4e-4 # [(g/kg)⁻¹] haline contraction
const g = 9.81 # [m/s²]

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
b = Field{Center, Center, Center}(grid);
horzAvgBuoy = Average(b, dims = (1,2));

## Compute the horizontally averaged profile for each time
datS = Array(ds["S"][:,:,:,:]);
datT = Array(ds["T"][:,:,:,:]);
datBuoy = CuArray(g * (α * datT - β * datS));
horzAvgTimeSeries = Array{Float64}(undef, Nt, Nz);
CUDA.allowscalar() do
	for ti in 1:Nt
		set!(b, datBuoy[:,:,:,ti]);
		horzAvgTimeSeries[ti, :] = compute!(Field(horzAvgBuoy))[1,1,1:Nz];
	end
end

## Create the heatmap
fig = Figure(size=(1000, 600));
hmAx = Axis(fig[1,1:2], title = "Horizontally averaged buoyancy [m/s²] over spin-up period", xlabel = "Time [y]", ylabel = "z [m]", limits = (nothing, (-1000, 0)));
hm = heatmap!(hmAx, t / year, zC, horzAvgTimeSeries, colormap = :balance);
Colorbar(fig[1,3], hm);
## Create profiles for the last three data points
zAx = Axis(fig[1,4], title = "Vertical profiles at end of spin-up", xlabel = "b [m/s²]", ylabel = "z [m]", limits = ((-0.230, -0.220), (-100, 0)));
timeLabels = ["Jan", "May", "Sep"]
timeIdx = [3, 2, 1]
for (j, label) in enumerate(timeLabels)
	idx = timeIdx[j];
	plt = lines!(zAx, horzAvgTimeSeries[(Nt-idx),:], zC[1:Nz], label = label)
end
axislegend()
save("./figures/horzAvgBuoy.png", fig);

close(ds)