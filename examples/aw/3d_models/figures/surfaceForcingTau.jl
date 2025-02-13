using NCDatasets
using CairoMakie

## Load in dataset, unpack
dsTau = NCDataset("./data/processed/winds.nc");
labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
y = dsTau["lat"][:];
t = dsTau["month"][:];
τ_x = dsTau["τ_x"][:,:];

## Create the base figure
fig = Figure(size=(800,800));

## Create axes and plot
for (j, label) in enumerate(labels)
    r, c = fld(j-1, 4) + 1, mod(j-1, 4) + 1
    ax = Axis(fig[r, c]; title = label, xlabel = "τₓ [N/m²]", ylabel = "Latitude [ᵒN]", limits = ((-0.25, 0.25), nothing))
    plot = lines!(ax, τ_x[:,j], y)
end
close(dsTau)

## Put title up
Label(fig[0,:], text = "Surface wind stress by month", fontsize = 24, font = :bold)

## Save figure
save("./figures/surfaceForcingTau.png", fig)

