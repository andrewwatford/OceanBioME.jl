using OceanBioME, Oceananigans, Printf, CairoMakie, Statistics, NCDatasets, Interpolations
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

# Define functions
# diffusivity
@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))
@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))
@inline MLD_corrected(t) = 0.67 * MLD(t) - 19

@inline κₜ(z, t) = 1e-2 * (1 + tanh((z - MLD_corrected(t)) / 10)) / 2 + 1e-4
# @inline κₜ(z, t) = ifelse(z > MLD_corrected(t), 1, 1e-3) * MLD_corrected(t).^2 / day # This is what is actually used in the Kuhn et al. paper

# Surface photosynthetically active radiation (PAR⁰)
@inline PAR⁰(t) = 80 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2;
@inline PAR(z, t) = PAR⁰(t) * exp(0.12 * z);

# Progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time));

# NPZD parameters
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
    open_bottom = true
); # From Kuhn et al. 2015, Table 2, column NA5

## Load in the data, interpolate
months_lst = ["january", "february", "march", "april", "may", "june", "july", "august", "september", "october", "november", "december"];
itp_times = []
depth_lst = []
mean_temp_lst = []
mean_nitrate_lst = []

for (j_month, m) = enumerate(months_lst)
    ds = NCDataset("./data/nc/$m.nc");
    depth = convert(Matrix{Float64}, ds["Depth"])
    temp = ds["Temperature"]
    nitrate = ds["Nitrate"]
    push!(itp_times, (j_month - 1) * month)
    push!(depth_lst, depth[:, 1])
    push!(mean_temp_lst, mean(temp, dims=2))
    push!(mean_nitrate_lst, mean(nitrate, dims=2))
end

# Have to add on another set of january to ensure correct looping
push!(itp_times, 12 * month)
push!(mean_temp_lst, mean_temp_lst[1])
push!(mean_nitrate_lst, mean_nitrate_lst[1])

nitrate_arr = reduce(hcat, mean_nitrate_lst)
temp_arr = reduce(hcat, mean_temp_lst)
depths = depth_lst[1]

# Create a gridded linear interpolation object
temp_itp = interpolate((depths, itp_times), temp_arr, Gridded(Linear()))
nitrate_itp = interpolate((depths, itp_times), nitrate_arr, Gridded(Linear()))

temp_function(z, t) = temp_itp(-z, mod(t, year))
nitrate_function(z, t) = nitrate_itp(-z, mod(t, year))

# Create the grid and clock
Nz = 100
Lz = 300meters
grid = RectilinearGrid(size = Nz, extent = Lz, topology = (Oceananigans.Flat, Oceananigans.Flat, Bounded));
clock = Clock(; time = 0.0)

file = "NPZD_forced_column_model.jld2"

# Create the background temperature and PAR fields
T = FunctionField{Center, Center, Center}(temp_function, grid; clock)
PAR_field = FunctionField{Center, Center, Center}(PAR, grid; clock)

# Create BGC model
biogeochemistry = NPZD(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(PAR_field), npzd_params...)

# Set up the forcing for the Nutrient (linear gradient from -250m to -10m)
nitrate_nudging_mask(z) = ifelse(z > -10, 0, ifelse(z > -250, (z + 10) / (-7200), 1/30))
nitrate_nudging = Relaxation(rate = 1 / days, mask = nitrate_nudging_mask, target = nitrate_function)
# nitrate_nudging = Relaxation(; rate = 0.001, target = nitrate_function)

# Create model, set ICs
model = NonhydrostaticModel(; grid,
                          clock = clock,
                          closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                          biogeochemistry = biogeochemistry,
                          auxiliary_fields = (; T),
                          forcing = (; N=nitrate_nudging))

set!(model; P = 0.2, Z = 0.4, N = 8.0, D = 0.7)

# Set up and run sim
simulation = Simulation(model, Δt = 1minutes, stop_time = 365days)
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields),
                                                    filename = file,
                                                    schedule = TimeInterval(1day),
                                                    overwrite_existing = true)
print("Starting simulation with Biology model ", biogeochemistry, "\n")
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
fig_aux = Figure(size = (1800, 500), fontsize = 24)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, times[end] / days), (-Lz, 0)))

# Mixed layer depth and diffusivity
axDiff = Axis(fig_aux[1, 1]; title = "Turbulent diffusivity (m² / s)", axis_kwargs...)
hmDiff = heatmap!(times / days, z,  κₜ.(z_grid, t_grid)', colormap = :balance, colorscale = log10)
Colorbar(fig_aux[1, 2], hmDiff)

# Attenuated PAR
axAPAR = Axis(fig_aux[1, 3]; title = "PAR (W / m²)", xlabel = "Time (days)", ylabel = "z (m)", 
                limits = ((0, times[end] / days), (-100meters, 0)))
hmAPAR = heatmap!(times / days, z, PAR.(z_grid, t_grid)', colormap = :balance)
Colorbar(fig_aux[1, 4], hmAPAR)

# Temperature
axTemp = Axis(fig_aux[1, 5]; title = "Temperature (ᵒC)", axis_kwargs...)
hmTemp = heatmap!(times / days, z, temp_function.(z_grid, t_grid)', colormap = :balance)
Colorbar(fig_aux[1, 6], hmTemp)

save("aux_plot.png", fig_aux);

# Plot the BGC fields over the time period
fig_bgc = Figure(size = (2400, 500), fontsize = 24)

axP = Axis(fig_bgc[1, 1]; title = "Phytoplankton (mmol N / m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap = :balance)
Colorbar(fig_bgc[1, 2], hmP)

axZ = Axis(fig_bgc[1, 3]; title = "Zooplankton (mmol N / m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)', colormap = :balance)
Colorbar(fig_bgc[1, 4], hmZ)

axN = Axis(fig_bgc[1, 5]; title = "Nutrient (mmol N / m³)", axis_kwargs...)
hmN = heatmap!(times / days, z, interior(N, 1, 1, :, :)', colormap = :balance)
Colorbar(fig_bgc[1, 6], hmN)

axD = Axis(fig_bgc[1, 7]; title = "Detritus (mmol N / m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(D, 1, 1, :, :)', colormap = :balance)
Colorbar(fig_bgc[1, 8], hmD)

save("bgc_plot.png", fig_bgc);