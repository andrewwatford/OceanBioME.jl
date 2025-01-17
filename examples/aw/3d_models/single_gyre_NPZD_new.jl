using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
using Oceananigans.Units
using OceanBioME
using Printf, Statistics, CairoMakie
using CUDA: @allowscalar

include("import_wind_data.jl")
include("import_WOA_data.jl")

const Lλ = 30
const Lφ = 15
const Lz = 4000

const λ_start = -50
const φ_start = 25

const Nλ = 20
const Nφ = 20
const Nz = 20

const Δt = 4minutes

const ρ₀ = 1035 # [kg m⁻³] reference density
const α = 2e-4
const β = 7.4e-4

const νh = 1e5 # [m²/s]
const νz = 3e-3
const κh = 1000
const κz = 1e-5

chebychev_spaced_z_faces(k) = -Lz + Lz * sin(π * (k - 1) / 2 / Nz);

@inline function τ_surface(i, j, grid, clock, fields)
    t = clock.time
    x = grid.λᶜᵃᵃ[i]
    y = grid.φᵃᶜᵃ[j]
    τ = wind_itp(x, y, mod(t, year) / month);
    return τ / ρ₀;
end

@inline function T_surface(i, j, grid, clock, fields)
    t = clock.time
    x = grid.λᶜᵃᵃ[i]
    y = grid.φᵃᶜᵃ[j]
    return temp_itp(x, y, mod(t, year) / month);
end

@inline function S_surface(i, j, grid, clock, fields)
    t = clock.time
    x = grid.λᶜᵃᵃ[i]
    y = grid.φᵃᶜᵃ[j]
    return sal_itp(x, y, mod(t, year) / month);
end

grid = LatitudeLongitudeGrid(CPU(); # was GPU()
                           size = (Nλ, Nφ, Nz),
                           halo = (6, 6, 6),
                      longitude = (λ_start,  λ_start + Lλ),
                       latitude = (φ_start,  φ_start + Lφ),
                              z = chebychev_spaced_z_faces,
                       topology = (Bounded, Bounded, Bounded))

u_surface_bc = FluxBoundaryCondition(τ_surface, discrete_form = true)
T_surface_bc = ValueBoundaryCondition(T_surface, discrete_form = true)
S_surface_bc = ValueBoundaryCondition(S_surface, discrete_form = true)

u_bcs = FieldBoundaryConditions(top = u_surface_bc)
T_bcs = FieldBoundaryConditions(top = T_surface_bc)
S_bcs = FieldBoundaryConditions(top = S_surface_bc)

horizontal_diffusive_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)
vertical_diffusive_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                                         ν = νz, κ = κz)

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = α,
                                                                   haline_contraction = β))

model = HydrostaticFreeSurfaceModel(; grid = grid, buoyancy,
                        momentum_advection = WENOVectorInvariant(),
                          tracer_advection = WENO(),
                                  coriolis = HydrostaticSphericalCoriolis(),
                                   closure = (vertical_diffusive_closure, horizontal_diffusive_closure),
                                   tracers = (:T, :S, ),
                       boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs,))

set!(model; )

function progress(sim)
    umax = maximum(abs, sim.model.velocities.u)
    vmax = maximum(abs, sim.model.velocities.v)
    Tmax = maximum(abs, sim.model.tracers.T)
    Smax = maximum(abs, sim.model.tracers.S)

 @info @sprintf("Iter: %d, time: %.2e, Δt: %s, max|u|: %.2e, max|v|: %.2e, max|T|: %.2e, max|S|: %.2e",
       iteration(sim), time(sim), prettytime(sim.Δt), umax, vmax, Tmax, Smax)

    return nothing
end

simulation = Simulation(model; Δt=Δt, stop_time= 400days)
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:p] = Callback(progress, TimeInterval(1hours))
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))

u, v, w = model.velocities
T, S = model.tracers
outputs = (; u, v, w, T, S)

simulation.output_writers[:fields] = NetCDFOutputWriter(
        model, outputs;
        filename = "single_gyre_fields.nc",
        schedule = TimeInterval(6hours),
        array_type = Array{Float32},
            overwrite_existing = true)

@time begin
    run!(simulation)
end