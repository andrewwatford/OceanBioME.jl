using NCDatasets
using CairoMakie

## Create the base figure
fig = Figure();

### Temperature plots ###
## Load in dataset, unpack
dsT = NCDataset("./data/processed/sst.nc");
labelsT = ["Jan", "Apr", "Jul", "Oct"];
y = dsT["lat"][:];
t = dsT["month"][:];
T = dsT["temp"][:,:];

## Create axis and plot
ax = Axis(fig[1, 1]; title = "Surface temperature by season", xlabel = "Surface temperature [ᵒC]", ylabel = "Latitude [ᵒN]")
for (j, label) in enumerate(labelsT)
    plt = lines!(ax, T[:,j], y, label = label)
end
axislegend()
close(dsT)

### Salinity plots ###
## Load in dataset, unpack
dsS = NCDataset("./data/processed/salinity.nc");
labelsS = ["Jan", "Apr", "Jul", "Oct"];
y = dsS["lat"][:];
t = dsS["month"][:];
S = dsS["salinity"][:,:];

## Create figure and plots
ax = Axis(fig[1, 2]; title = "Surface salinity by season", xlabel = "Surface salinity [psu]", ylabel = "Latitude [ᵒN]")
for (j, label) in enumerate(labelsS)
    plt = lines!(ax, S[:,j], y, label = label)
end
axislegend()
close(dsS)

## Save temperature and salinity figure
save("./figures/surfaceForcingTS.png", fig)

