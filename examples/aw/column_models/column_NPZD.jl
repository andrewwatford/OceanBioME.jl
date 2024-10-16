using OceanBioME, Oceananigans, Printf
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

using CairoMakie

const year = years = 365days;

Nz = 60
Lz = 300meters
grid = RectilinearGrid(size = Nz, extent = Lz, topology = (Flat, Flat, Bounded));

# Inline functions for diffusivity
@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))
@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))
@inline MLD_corrected(t) = 0.67 * MLD(t) - 19

@inline κₜ(z, t) = 1e-2 * (1 + tanh((z - MLD_corrected(t)) / 10)) / 2 + 1e-4
# @inline κₜ(z, t) = ifelse(z > MLD_corrected(t), 1, 1e-3) * MLD_corrected(t).^2 / day # This is what is actually used in the Kuhn et al. paper

# Temperature
@inline temp(z, t) = 2.4 * cos(t * 2π / year + 50days) + 14

# Surface photosynthetically active radiation (PAR⁰)
@inline PAR⁰(t) = 80 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2;
@inline PAR(z, t) = PAR⁰(t) * exp(0.12 * z);

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time));

npzd_params = (
    initial_photosynthetic_slope = 0.2099/day,
    base_maximum_growth = 0.7894/day,
    nutrient_half_saturation = 3.4151,
    base_respiration_rate = 0.0062/day, 
    phyto_base_mortality_rate = 0.0109/day, 
    maximum_grazing_rate = 2.0811/day,
    grazing_half_saturation = 0.5373, # may be off, they supply this as the square in Table 2...
    assimulation_efficiency = 0.8781,
    base_excretion_rate = 0.0133/day,
    zoo_base_mortality_rate = 0.3998/day,
    remineralization_rate = 0.1401/day,
    sinking_speeds = (; P = 0.1117/day, D = 5.0400/day),
    scale_negatives = true,
    open_bottom = false
); # From Kuhn et al. 2015, Table 2, column NA5

file = "NPZD_column_model.jld2"
    
clock = Clock(; time = 0.0) 
T = FunctionField{Center, Center, Center}(temp, grid; clock)
PAR_field = FunctionField{Center, Center, Center}(PAR, grid; clock)

biogeochemistry = NPZD(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(PAR_field), npzd_params...)
print("Simulation with Biology model ", biogeochemistry, "\n")

model = NonhydrostaticModel(; grid,
                          clock = clock,
                          closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                          biogeochemistry = biogeochemistry,
                          auxiliary_fields = (; T))

set!(model; P = 0.03, Z = 0.03, N = 10.0)

simulation = Simulation(model, Δt = 3minutes, stop_time = 365days)
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields),
                                                    filename = file,
                                                    schedule = TimeInterval(1day),
                                                    overwrite_existing = true)
run!(simulation)

# Read in the model outputs
N = FieldTimeSeries(file, "N")
P = FieldTimeSeries(file, "P")
Z = FieldTimeSeries(file, "Z")
D = FieldTimeSeries(file, "D")
x, y, z = nodes(P)
times = P.times;

# Create a grid for z and t
z_grid = repeat(z, 1, length(times));      # Shape z into a matrix
t_grid = repeat(times', length(z), 1);     # Shape t into a matrix

# Plot the auxiliary fields over the time period
fig_aux = Figure(size = (1000, 800), fontsize = 20)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, times[end] / days), (-200meters, 0)))

# Mixed layer depth and diffusivity
axMLD = Axis(fig_aux[1, 1]; title = "Mixed layer depth (m)", xlabel = "Time (days)",
            ylabel = "MLD (m)", limits = ((0, times[end] / days), nothing))
linesMLD = lines!(axMLD, times / days, MLD_corrected.(times))

axDiff = Axis(fig_aux[1, 3]; title = "Scalar diffusivity (m² / s)", axis_kwargs...)
hmDiff = heatmap!(times / days, z,  κₜ.(z_grid', t_grid'), colormap = :batlow, colorscale = log10)
Colorbar(fig_aux[1, 4], hmDiff)

# Attenuated PAR
axAPAR = Axis(fig_aux[2, 1]; title = "PAR (W / m²)", xlabel = "Time (days)", ylabel = "z (m)", 
                limits = ((0, times[end] / days), (-100meters, 0)))
hmAPAR = heatmap!(times / days, z, PAR.(z_grid', t_grid'), colormap = :batlow)
Colorbar(fig_aux[2, 2], hmAPAR)

# Temperature
axTemp = Axis(fig_aux[2, 3]; title = "Temperature (ᵒC)", axis_kwargs...)
hmTemp = heatmap!(times / days, z, temp.(z_grid', t_grid'), colormap = :batlow)
Colorbar(fig_aux[2, 4], hmTemp)

save("aux_plot.png", fig_aux);

# Plot the BGC fields over the time period
fig_bgc = Figure(size = (1000, 600), fontsize = 20)

axP = Axis(fig_bgc[1, 1]; title = "Phytoplankton (mmol N / m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig_bgc[1, 2], hmP)

axZ = Axis(fig_bgc[1, 3]; title = "Zooplankton (mmol N / m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig_bgc[1, 4], hmZ)

axN = Axis(fig_bgc[2, 1]; title = "Nutrient (mmol N / m³)", axis_kwargs...)
hmN = heatmap!(times / days, z, interior(N, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig_bgc[2, 2], hmN)

axD = Axis(fig_bgc[2, 3]; title = "Detritus (mmol N / m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(D, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig_bgc[2, 4], hmD)

save("bgc_plot.png", fig_bgc);

# Plot total nitrogen content
Ntot = interior(P, 1, 1, :, :)' + interior(N, 1, 1, :, :)' + interior(Z, 1, 1, :, :)' + interior(D, 1, 1, :, :)';

fig_N = Figure(size = (1000, 600), fontsize = 20)
axN = Axis(fig_N[1, 1]; title = "Total Nitrate (mmol N / m³)", xlabel = "Time (days)", ylabel = "z (m)", 
            limits = ((0.0, times[end] / days), (-Lz, 0)))
hmN = heatmap!(times / days, z, Ntot, colormap = :batlow)
Colorbar(fig_N[1, 2], hmN)

save("N_plot.png", fig_N);