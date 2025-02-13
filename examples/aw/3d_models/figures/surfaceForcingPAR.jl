using NCDatasets
using CairoMakie

const hour = 60 * 60;
const day = 24hour;
const year = 365day;
const month = year / 12;

const Lφ = 15
const φ_start = 25

### Define PAR consistent with NPZD_single_gyre_gpu.jl !!!!! ###
@inline PAR⁰(t) = 90 * (1 - cos((t + 15day) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200day) / 50day)^2)));
@inline PAR(y, t) = PAR⁰(t) * (1 - 0.1 * (y - φ_start) / Lφ);

t = (0:0.1:365)* day
y = (φ_start:0.01:(φ_start+Lφ))
tgrid = repeat(t, 1, length(y))
ygrid = repeat(y, 1, length(t))'
par = PAR.(ygrid, tgrid)

fig = Figure();
hmAx = Axis(fig[1,1], title = "Photosynthetically active radiation [W/m²]", xlabel = "Time [d]", ylabel = "Latitude [ᵒN]");
hm = heatmap!(hmAx, t / day, y, par, colormap = :balance);
Colorbar(fig[1,2], hm);

## Save figure
save("./figures/surfaceForcingPAR.png", fig)

