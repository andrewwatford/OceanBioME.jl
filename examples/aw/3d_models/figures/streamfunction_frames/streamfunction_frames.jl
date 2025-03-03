using Oceananigans
using NCDatasets
using CairoMakie
using CUDA
using Statistics

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

grid = LatitudeLongitudeGrid(GPU(), Float32;
                           size = (Nλ, Nφ, Nz),
                           halo = (6, 6, 6),
                      longitude = xF,
                       latitude = yF,
                              z = zF,
                       topology = (Bounded, Bounded, Bounded));

## Predefine the relevant computations
uField = Field{Face, Center, Center}(grid);
vField = Field{Center, Face, Center}(grid);
uCumField = Field(CumulativeIntegral(uField, dims=2));
vCumField = Field(CumulativeIntegral(vField, dims=1));
uBtFieldCum = Field(Integral(uCumField, dims=3));
vBtFieldCum = Field(Integral(vCumField, dims=3));

## Compute for first year
for i = 1:1:36
    print("Now on iteration $(i)\n")
    u = CuArray(ds["u"][:,:,:,i]);
    v = CuArray(ds["v"][:,:,:,i]);
    set!(uField, u);
    set!(vField, v);
    uBtCum, vBtCum = compute!(uBtFieldCum), compute!(vBtFieldCum);
    # Put data onto center grid
    CUDA.@allowscalar uBtCumCenter = (uBtCum[1:Nλ,1:Nφ,1] + uBtCum[2:Nλ+1,1:Nφ,1]) / 2;
    CUDA.@allowscalar vBtCumCenter = (vBtCum[1:Nλ,1:Nφ,1] + vBtCum[1:Nλ,2:Nφ+1,1]) / 2;
    # Compute streamfunction
    CUDA.@allowscalar psi1 = Array(vBtCumCenter[:,1] .- uBtCumCenter[:,:]);
    CUDA.@allowscalar psi2 = Array(vBtCumCenter[:,:] .- uBtCumCenter[1,:]);
    psiAvg = 0.5 * (psi1 + psi2);
    psiDiff = psi2 - psi1

    ### Figures
    fig = Figure();
    ax = Axis(fig[1,1]; title = "Barotropic streamfunction (averaged) [m²/s]", xlabel = "Longitude [ᵒE]", ylabel = "Latitude [ᵒN]");
    co = contourf!(ax, xC, yC, psiAvg);
    Colorbar(fig[1,2], co);
    save("./examples/aw/3d_models/figures/streamfunction_frames/streamfunction_frame$(i).png", fig)
end