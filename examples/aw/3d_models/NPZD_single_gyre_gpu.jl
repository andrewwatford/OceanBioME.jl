using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: FunctionField
using Printf
using Interpolations, NCDatasets
using OceanBioME
using CUDA

const year = years = 365days;
const month = months = year / 12;

const Lλ = 20
const Lφ = 15
const Lz = 4000

const λ_start = -45
const φ_start = 25

const Nλ = 200
const Nφ = 200
const Nz = 200
const halo_size = 6

const Δt = 5minutes

const ρ₀ = 1035 # [kg m⁻³] reference density
const α = 2e-4  # [ᵒC⁻¹] thermal expansion
const β = 7.4e-4 # [(g/kg)⁻¹] haline contraction

const νh = 1e3 # [m²/s]
const νz = 3e-3
const κh = 1e3
const κz = 1e-5

## Create the grid
chebychev_spaced_z_faces(k) = -Lz + Lz * sin(π * (k - 1) / 2 / Nz);
grid = LatitudeLongitudeGrid(GPU();
                           size = (Nλ, Nφ, Nz),
                           halo = (halo_size, halo_size, halo_size),
                      longitude = (λ_start,  λ_start + Lλ),
                       latitude = (φ_start,  φ_start + Lφ),
                              z = chebychev_spaced_z_faces,
                       topology = (Bounded, Bounded, Bounded))

## Boundary conditions
@inline current_time_index(time)                 = mod(unsafe_trunc(Int32, time / month), 12) + 1;
@inline next_time_index(time)                    = mod(unsafe_trunc(Int32, time / month) + 1, 12) + 1;
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / month, 1) * (u₂ - u₁);

# Wind stress
wind_ds = NCDataset("./data/processed/winds.nc", share=true);
wind_itp = interpolate((wind_ds["lat"], wind_ds["month"]), wind_ds["τ_x"], Gridded(Linear()));
const τ_on_grid = CuArray(wind_itp.(repeat(grid.φᵃᶜᵃ[1:Nφ], 1, 12), (0:11)'));
close(wind_ds);

@inline function τ_surface(i, j, grid, clock, fields, p)
    time = clock.time
    n₁ = current_time_index(time)
    n₂ = next_time_index(time)
    τ₁ = p[j, n₁]
    τ₂ = p[j, n₂]
    τ = cyclic_interpolate(τ₁, τ₂, time);
    return - τ / ρ₀
end

# Temperature
temp_ds = NCDataset("./data/processed/sst.nc", share=true);
temp_itp = interpolate((temp_ds["lat"], temp_ds["month"]), temp_ds["temp"], Gridded(Linear()));
const T_on_grid = CuArray(temp_itp.(repeat(grid.φᵃᶜᵃ[1:Nφ], 1, 12), (0:11)'));
close(temp_ds);

@inline function T_surface(i, j, grid, clock, fields, p)
    time = clock.time
    n₁ = current_time_index(time)
    n₂ = next_time_index(time)
    T₁ = p[j, n₁]
    T₂ = p[j, n₂]
    T = cyclic_interpolate(T₁, T₂, time);
    return T
end

# Salinity
sal_ds = NCDataset("./data/processed/salinity.nc", share=true);
sal_itp = interpolate((sal_ds["lat"], sal_ds["month"]), sal_ds["salinity"], Gridded(Linear()));
const S_on_grid = CuArray(sal_itp.(repeat(grid.φᵃᶜᵃ[1:Nφ], 1, 12), (0:11)'));
close(sal_ds);

@inline function S_surface(i, j, grid, clock, fields, p)
    time = clock.time
    n₁ = current_time_index(time)
    n₂ = next_time_index(time)
    S₁ = p[j, n₁]
    S₂ = p[j, n₂]
    S = cyclic_interpolate(S₁, S₂, time);
    return S
end

u_surface_bc = FluxBoundaryCondition(τ_surface, discrete_form = true, parameters = τ_on_grid)
T_surface_bc = ValueBoundaryCondition(T_surface, discrete_form = true, parameters = T_on_grid)
S_surface_bc = ValueBoundaryCondition(S_surface, discrete_form = true, parameters = S_on_grid)

u_bcs = FieldBoundaryConditions(top = u_surface_bc)
T_bcs = FieldBoundaryConditions(top = T_surface_bc)
S_bcs = FieldBoundaryConditions(top = S_surface_bc)

## Diffusive closures
horizontal_diffusive_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)
vertical_diffusive_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                                         ν = νz, κ = κz)

## Buoyancy
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = α,
                                                                   haline_contraction = β))

## Defining BGC setup
@inline PAR⁰(t) = 90 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2)));
@inline PAR(x, y, t) = PAR⁰(t) * (1 - 0.1 * (y - φ_start) / Lφ);
biogeochemistry = NPZD(; grid, surface_photosynthetically_active_radiation = PAR, scale_negatives=true)

model = HydrostaticFreeSurfaceModel(; grid = grid, buoyancy,
                        momentum_advection = WENOVectorInvariant(),
                        tracer_advection = WENO(),
                        timestepper = :SplitRungeKutta3,
                        coriolis = HydrostaticSphericalCoriolis(),
                        biogeochemistry = biogeochemistry,
                        closure = (vertical_diffusive_closure, horizontal_diffusive_closure),
                        tracers = (:N, :P, :Z, :D, :T, :S, ),
                        boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs,))

### Here, I am setting the initial conditions based on the end of the previous run
# I also set initial values for the NPZD variables.
prev_ds = NCDataset("./spinup.nc")
prevT = convert(CuArray{Float64}, prev_ds["T"][:,:,:,end])
prevS = convert(CuArray{Float64}, prev_ds["S"][:,:,:,end])
U, V, W =  convert(CuArray{Float64}, prev_ds["u"][:,:,:,end]), convert(CuArray{Float64}, prev_ds["v"][:,:,:,end]), convert(CuArray{Float64}, prev_ds["w"][:,:,:,end])
set!(model; T = prevT, S = prevS, u = U, v = V, w = W, N = 10, P = 0.1, Z = 0.1, D = 0)
close(prev_ds)

function progress(sim)
    umax = maximum(abs, sim.model.velocities.u)
    vmax = maximum(abs, sim.model.velocities.v)
    Tmax = maximum(abs, sim.model.tracers.T)
    Smax = maximum(abs, sim.model.tracers.S)

 @info @sprintf("Iter: %d, time: %.2e, Δt: %s, max|u|: %.2e, max|v|: %.2e, max|T|: %.2e, max|S|: %.2e",
       iteration(sim), time(sim), prettytime(sim.Δt), umax, vmax, Tmax, Smax)

    return nothing
end

simulation = Simulation(model; Δt=Δt, stop_time=10years)
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:p] = Callback(progress, TimeInterval(10days))
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

u, v, w = model.velocities
ω = ∂x(v) - ∂y(u)
N, P, Z, D, T, S = model.tracers
outputs = (; u, v, w, T, S, ω, N, P, Z, D);

simulation.output_writers[:fields] = NetCDFOutputWriter(
        model, outputs;
        filename = "NPZD_spinup.nc",
        schedule = TimeInterval(1month),
        array_type = Array{Float32},
        overwrite_existing = true)

@time begin
    run!(simulation)
end
