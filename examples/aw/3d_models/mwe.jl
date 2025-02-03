using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
using Oceananigans.Units
using Printf, Statistics
using CUDA: @allowscalar

const year = years = 365days;
const month = months = year / 12;

const Lλ = 20
const Lφ = 15
const Lz = 4000

const λ_start = -45
const φ_start = 25

const Nλ = 50
const Nφ = 50
const Nz = 3

const Δt = 4minutes

const ρ₀ = 1035 # [kg m⁻³] reference density
const α = 2e-4  # [ᵒC⁻¹] thermal expansion
const β = 7.4e-4 # [(g/kg)⁻¹] haline contraction

const νh = 1e3 # [m²/s]
const νz = 3e-3
const κh = 1e3
const κz = 1e-5

chebychev_spaced_z_faces(k) = -Lz + Lz * sin(π * (k - 1) / 2 / Nz);

@inline function τ_surface(i, j, grid, clock, fields)
    y = grid.φᵃᶜᵃ[j]
    τ = 0.01 .* (y .- 33)
    return - τ ./ ρ₀;
end

@inline function T_surface(i, j, grid, clock, fields)
    t = clock.time
    y = grid.φᵃᶜᵃ[j]
    return 0.25 .* (y .- 33) .+ 20
end

@inline function S_surface(i, j, grid, clock, fields)
    t = clock.time
    y = grid.φᵃᶜᵃ[j]
    return 0.075 .* (y .- 35) .+ 36.7
end

grid = LatitudeLongitudeGrid(GPU();
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

set!(model; T = 12, S = 36)

function progress(sim)
    umax = maximum(abs, sim.model.velocities.u)
    vmax = maximum(abs, sim.model.velocities.v)
    Tmax = maximum(abs, sim.model.tracers.T)
    Smax = maximum(abs, sim.model.tracers.S)

 @info @sprintf("Iter: %d, time: %.2e, Δt: %s, max|u|: %.2e, max|v|: %.2e, max|T|: %.2e, max|S|: %.2e",
       iteration(sim), time(sim), prettytime(sim.Δt), umax, vmax, Tmax, Smax)

    return nothing
end

simulation = Simulation(model; Δt=Δt, stop_time=10days)
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:p] = Callback(progress, TimeInterval(6hours))
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))

u, v, w = model.velocities
ω = ∂x(v) - ∂y(u)
T, S = model.tracers
outputs = (; u, v, w, T, S, ω)

simulation.output_writers[:fields] = NetCDFOutputWriter(
        model, outputs;
        filename = "single_gyre_fields.nc",
        schedule = TimeInterval(6hours),
        array_type = Array{Float32},
            overwrite_existing = true)

@time begin
    run!(simulation)
end
