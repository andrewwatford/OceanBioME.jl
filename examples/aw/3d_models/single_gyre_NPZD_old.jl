using Oceananigans
using Oceananigans.Grids: λnodes, φnodes
using Oceananigans.Units
using OceanBioME
using Printf, Statistics, CairoMakie
using CUDA: @allowscalar

include("FunctionsLibrary.jl")
using .FunctionsLibrary

const Lλ = 25
const Lφ = 5
const Lz = 4000
const λ_start = -60
const φ_start = 25
const λ0 = -25
const φ0 = 20
const Nλ = 10
const Nφ = 10
const Nz = 10

const Δt = 250
const stop_time = 1days

const ρ₀ = 1035 # [kg m⁻³] reference density

const νh = 100000    
const νz = 3e-3
const κh = 1000
const κz = 1e-5

chebychev_spaced_z_faces(k) = -Lz + Lz * sin(π * (k - 1) / 2 / Nz);

grid = LatitudeLongitudeGrid(CPU(); # was GPU()
                           size = (Nλ, Nφ, Nz),
                           halo = (6, 6, 6),
                      longitude = (λ_start,  λ_start + Lλ),
                       latitude = (φ_start,  φ_start + Lφ),
                              z = chebychev_spaced_z_faces,
                       topology = (Bounded, Bounded, Bounded))

parameters = (Lφ = Lφ,
              Lz = Lz,
      T_seasonal = 360days, 
    τ_amp_winter = 0.1 / ρ₀,                   
    τ_amp_summer = 0.08 / ρ₀,
        τ_period = 1.9,
        τ_offset = -0.81,
         sst_amp = 14.5,
      sst_period = 0.6,
      sst_offset = 0.0,
  sal_amp_winter = 5.07e-8,
  sal_amp_summer = 4.44e-8,
      sal_period = 0.9,
      sal_offset = 3.1,
  sst_min_winter = 11,
  sst_min_summer = 13.8)

u_surface_bc = FluxBoundaryCondition(τ_surface,  discrete_form = true, parameters = parameters)
T_surface_bc = ValueBoundaryCondition(T_surface, discrete_form = true, parameters = parameters)
S_surface_bc = ValueBoundaryCondition(S_surface, discrete_form = true, parameters = parameters)

u_bcs = FieldBoundaryConditions(top = u_surface_bc)
T_bcs = FieldBoundaryConditions(top = T_surface_bc)
S_bcs = FieldBoundaryConditions(top = S_surface_bc)

horizontal_diffusive_closure = HorizontalScalarDiffusivity(ν = νh, κ = κh)
vertical_diffusive_closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                                         ν = νz, κ = κz)

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4,
                                                                   haline_contraction = 7.4e-4))
biogeochemistry = NPZD(; grid, scale_negatives = true)

model = HydrostaticFreeSurfaceModel(; grid = grid, buoyancy, biogeochemistry, 
                        momentum_advection = WENOVectorInvariant(),
                          tracer_advection = WENO(),
                                  coriolis = HydrostaticSphericalCoriolis(),
                                   closure = (vertical_diffusive_closure, horizontal_diffusive_closure),
                                   tracers = (:S, ),
                       boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs,))

set!(model; P = 0.1, Z = 0.1, N = 10.0, D = 0)

function progress(sim)
    umax = maximum(abs, sim.model.velocities.u)
    vmax = maximum(abs, sim.model.velocities.v)
    Tmax = maximum(abs, sim.model.tracers.T)
    Smax = maximum(abs, sim.model.tracers.S)

 @info @sprintf("Iter: %d, time: %.2e, Δt: %s, max|u|: %.2e, max|v|: %.2e, max|T|: %.2e, max|S|: %.2e",
       iteration(sim), time(sim), prettytime(sim.Δt), umax, vmax, Tmax, Smax)

    return nothing
end

simulation = Simulation(model; Δt=Δt, stop_time= 100days)
wizard = TimeStepWizard(cfl=0.2, max_change=1.1, max_Δt=20minutes)
simulation.callbacks[:p] = Callback(progress, TimeInterval(1hours))
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(2))

u, v, w = model.velocities
S, N, P, Z, D, T = model.tracers
outputs = (; u, v, T, S, N, P, Z, D)

simulation.output_writers[:fields] = NetCDFOutputWriter(
        model, outputs;
        filename = "single_gyre_fields.nc",
        schedule = TimeInterval(6hours),
        array_type = Array{Float32},
            overwrite_existing = true)

@time begin
    run!(simulation)
end

plot_forcing_profiles(grid, model, parameters, τ_surface, T_surface, S_surface)
plot_streamlines(grid)
plot_velocity(grid)
plot_velocity_BT(grid)
