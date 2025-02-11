using Printf, CairoMakie, NCDatasets, Statistics

ds = NCDataset("single_gyre_gpu_fields.nc");
xC = Vector(ds["xC"])
yC = Vector(ds["yC"])
zF = Vector(ds["zF"]);
dz = zF[2:end] .- zF[1:end-1]
ω = Array(ds["ω"]);
ω_BT = sum(ω .* reshape(dz, (1, 1, :, 1)), dims=3)[:,:,1,:]

max_abs = maximum(abs, ω_BT[:,:,end])
colorrange = (-max_abs, max_abs)

fig = Figure()
ax = Axis(fig[1,1])
hm = heatmap!(xC, yC, ω_BT[:,:,end], colormap = :balance, colorrange = colorrange)
Colorbar(fig[1,2], hm)

save("./figures/BT_vorticity.png", fig)

close(ds)
