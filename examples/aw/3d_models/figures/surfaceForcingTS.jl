using NCDatasets
using CairoMakie

## Create the base figure
fig = Figure();

### Temperature plots ###
## Load in dataset, unpack
dsT = NCDataset("./examples/aw/3d_models/data/processed/sst.nc");
labelsT = ["Jan", "Apr", "Jul", "Oct"];
y = dsT["lat"][:];
t = dsT["month"][:];
T = dsT["temp"][:,:];

## Create axis and plot
ax = Axis(fig[1, 1]; title = "Surface temperature by season", xlabel = "Surface temperature [ᵒC]", ylabel = "Latitude [ᵒN]", limits = ((15, 30), (25, 40)))
for (j, label) in enumerate(labelsT)
    plt = lines!(ax, T[:,j], y, label = label)
end
axislegend()
close(dsT)

### Salinity plots ###
## Load in dataset, unpack
dsS = NCDataset("./examples/aw/3d_models/data/processed/salinity.nc");
labelsS = ["Jan", "Apr", "Jul", "Oct"];
y = dsS["lat"][:];
t = dsS["month"][:];
S = dsS["salinity"][:,:];

## Create figure and plots
ax = Axis(fig[1, 2]; title = "Surface salinity by season", xlabel = "Surface salinity [g/kg]", ylabel = "Latitude [ᵒN]", limits = ((35.5, 37.5), (25, 40)))
for (j, label) in enumerate(labelsS)
    plt = lines!(ax, S[:,j], y, label = label)
end
axislegend()
close(dsS)

## Save temperature and salinity figure
save("./examples/aw/3d_models/figures/surfaceForcingTS.png", fig)

