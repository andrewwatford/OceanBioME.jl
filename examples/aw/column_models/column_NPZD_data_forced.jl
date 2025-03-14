using OceanBioME, Oceananigans, Printf, CairoMakie, Statistics, NCDatasets, Interpolations
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

# Create the MLD interpolator
mld_ds = NCDataset("./data/nc/mixed_layer_processed.nc")
mld_dat = convert(Array{Float64}, mld_ds["depth_mean"])
push!(mld_dat, mld_dat[1]) # to make looping work
mld_itp_months = interpolate(-mld_dat, BSpline(Cubic(Line(OnGrid()))))
mld_itp(t) = mld_itp_months(mod(t / month, 12) + 1)

# Define diffusivity based on MLD
@inline κₜ(z, t) = 1e-2 * (1 + tanh((z - mld_itp(t)) / 10)) / 2 + 1e-4

# Surface photosynthetically active radiation (PAR⁰)
@inline PAR⁰(t) = 90 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2)));
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
    grazing_half_saturation = 0.5373, # may be off, they supply this as the square in Table 2...but changing it doesn't affect much
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
Nz = 60
Lz = 300meters
grid = RectilinearGrid(size = Nz, extent = Lz, topology = (Oceananigans.Flat, Oceananigans.Flat, Bounded));
clock = Clock(; time = 0.0)

file = "debug.jld2"

# Create the background temperature and PAR fields
T_field = FunctionField{Center, Center, Center}(temp_function, grid; clock)
PAR_field = FunctionField{Center, Center, Center}(PAR, grid; clock)

# Create BGC model
biogeochemistry = NPZD(; grid, light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(PAR_field), npzd_params...)

# Set up the forcing for the Nutrient (linear gradient from -250m to -10m)
nitrate_nudging_mask(z) = ifelse(z > -10, 0, ifelse(z > -250, (z + 10) / (-7200), 1/30))
nitrate_nudging = Relaxation(rate = 1 / days, mask = nitrate_nudging_mask, target = nitrate_function)

# Create model, set ICs
model = NonhydrostaticModel(; grid,
                          clock = clock,
                          closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                          biogeochemistry = biogeochemistry,
                          auxiliary_fields = (; T = T_field),
                          forcing = (; N = nitrate_nudging))

init_N_dist(z) = nitrate_function(z, 0)
set!(model; P = 0.1, Z = 0.1, N = init_N_dist, D = 0)

# Set up and run sim
dt = 10minutes
simulation = Simulation(model, Δt = dt, stop_time = 200days)
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields),
                                                    filename = file,
                                                    schedule = TimeInterval(10day),
                                                    overwrite_existing = true)
print("Starting simulation with Biology model ", biogeochemistry, "\n")
run!(simulation)

# Read in the model outputs
N = FieldTimeSeries(file, "N")
P = FieldTimeSeries(file, "P")
Z = FieldTimeSeries(file, "Z")
D = FieldTimeSeries(file, "D")
PAR_dat = FieldTimeSeries(file, "PAR")
T_dat = FieldTimeSeries(file, "T")
x, y, z = nodes(PAR_dat)
times = PAR_dat.times;

# Create a grid for z and t
z_grid = repeat(z, 1, length(times));      # Shape z into a matrix
t_grid = repeat(times', length(z), 1);     # Shape t into a matrix

# Plot the auxiliary fields over the time period
fig = Figure(size = (2400, 1000), fontsize = 24)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, times[end] / days), (z[1], 0)))

# Mixed layer depth and diffusivity
axDiff = Axis(fig[1, 1]; title = "Turbulent diffusivity (m² / s)", axis_kwargs...)
hmDiff = heatmap!(times / days, z,  κₜ.(z_grid, t_grid)', colormap = :balance, colorscale = log10)
Colorbar(fig[1, 2], hmDiff)

# Attenuated PAR
axAPAR = Axis(fig[1, 3]; title = "PAR (W / m²)", xlabel = "Time (days)", ylabel = "z (m)", 
                limits = ((0, times[end] / days), (-100meters, 0)))
hmAPAR = heatmap!(times / days, z, interior(PAR_dat, 1, 1, :, :)', colormap = :balance)
Colorbar(fig[1, 4], hmAPAR)

# Temperature
axTemp = Axis(fig[1, 5]; title = "Temperature (ᵒC)", axis_kwargs...)
hmTemp = heatmap!(times / days, z, interior(T_dat, 1, 1, :, :)', colormap = :balance)
Colorbar(fig[1, 6], hmTemp)

# Light limitation
@inline Q₁₀(T) = 1.88 ^ (T / 10) # T in °C
@inline light_limitation(PAR, α, μ₀) = α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)
ll_dat = light_limitation.(PAR_dat, npzd_params.initial_photosynthetic_slope, npzd_params.base_maximum_growth .* Q₁₀.(T_dat))
axLL = Axis(fig[1, 7]; title = "Light Limitation", xlabel = "Time (days)", ylabel = "z (m)", 
            limits = ((0, times[end] / days), (-100meters, 0)))
hmLL = heatmap!(times / days, z, ll_dat[1,1,:,:]', colormap = :balance)
Colorbar(fig[1, 8], hmLL)

# Plot the BGC fields over the time period
axP = Axis(fig[2, 1]; title = "Phytoplankton (mmol N / m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap = :balance)
Colorbar(fig[2, 2], hmP)

axZ = Axis(fig[2, 3]; title = "Zooplankton (mmol N / m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)', colormap = :balance)
Colorbar(fig[2, 4], hmZ)

axN = Axis(fig[2, 5]; title = "Nutrient (mmol N / m³)", axis_kwargs...)
hmN = heatmap!(times / days, z, interior(N, 1, 1, :, :)', colormap = :balance)
Colorbar(fig[2, 6], hmN)

axD = Axis(fig[2, 7]; title = "Detritus (mmol N / m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(D, 1, 1, :, :)', colormap = :balance)
Colorbar(fig[2, 8], hmD)

save("plot.png", fig);