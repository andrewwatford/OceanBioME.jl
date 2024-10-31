using Oceananigans, CairoMakie
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

files = ["no_T_dependent.jld2", "growth_T_dependent.jld2", "death_T_dependent.jld2", "all_T_dependent.jld2"]
cases = ["Control (T = 12ᵒC)", "T-dependent Growth", "T-dependent Mortality", "All T-dependent"]
data = []

for (j, file) = enumerate(files)
    # Read in the model outputs
    N = FieldTimeSeries(file, "N")
    P = FieldTimeSeries(file, "P")
    Z = FieldTimeSeries(file, "Z")
    D = FieldTimeSeries(file, "D")
    data_tuple = (; N = N, P = P, Z = Z, D = D)
    push!(data, data_tuple)
end

controlN = data[1][:N]
controlP = data[1][:P]
controlZ = data[1][:Z]
controlD = data[1][:D]
x, y, z = nodes(controlN)
times = controlN.times;
# Create a grid for z and t
z_grid = repeat(z, 1, length(times));      # Shape z into a matrix
t_grid = repeat(times', length(z), 1);     # Shape t into a matrix

fig = Figure(size = (2400, 2000), fontsize = 24)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0.0, times[end] / days), (z[1], 0)))
# Plot the control
axControlN = Axis(fig[1, 1]; title = "N (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmControlN = heatmap!(times / days, z, controlN[1,1,:,:]', colormap = :balance)
Colorbar(fig[1, 2], hmControlN)
axControlP = Axis(fig[2, 1]; title = "P (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmControlP = heatmap!(times / days, z, controlP[1,1,:,:]', colormap = :balance)
Colorbar(fig[2, 2], hmControlP)
axControlZ = Axis(fig[3, 1]; title = "Z (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmControlZ = heatmap!(times / days, z, controlZ[1,1,:,:]', colormap = :balance)
Colorbar(fig[3, 2], hmControlZ)
axControlD = Axis(fig[4, 1]; title = "D (mmol N / m³) for Control (T = 12ᵒC)", axis_kwargs...)
hmControlD = heatmap!(times / days, z, controlD[1,1,:,:]', colormap = :balance)
Colorbar(fig[4, 2], hmControlD)

for (j, data_tuple) = enumerate(data)
    if j != 1
        N = data_tuple[:N]
        P = data_tuple[:P]
        Z = data_tuple[:Z]
        D = data_tuple[:D]
        # Plot the residuals
        label = cases[j]
        axN = Axis(fig[1, 1+j]; title = "ΔN (mmol N / m³) for $label", axis_kwargs...)
        hmN = heatmap!(times / days, z, N[1,1,:,:]' .- controlN[1,1,:,:]', colormap = :balance, colorrange = (-0.75, 0.75))
        axP = Axis(fig[2, 1+j]; title = "ΔP (mmol N / m³) for $label", axis_kwargs...)
        hmP = heatmap!(times / days, z, P[1,1,:,:]' .- controlP[1,1,:,:]', colormap = :balance, colorrange = (-0.07, 0.07))
        axZ = Axis(fig[3, 1+j]; title = "ΔZ (mmol N / m³) for $label", axis_kwargs...)
        hmZ = heatmap!(times / days, z, Z[1,1,:,:]' .- controlZ[1,1,:,:]', colormap = :balance, colorrange = (-0.05, 0.05))
        axD = Axis(fig[4, 1+j]; title = "ΔD (mmol N / m³) for $label", axis_kwargs...)
        hmD = heatmap!(times / days, z, D[1,1,:,:]' .- controlD[1,1,:,:]', colormap = :balance, colorrange = (-0.2, 0.2))
        if j == 4
            Colorbar(fig[1, 6], hmN)
            Colorbar(fig[2, 6], hmP)
            Colorbar(fig[3, 6], hmZ)
            Colorbar(fig[4, 6], hmD)
        end
    end
end
    
save("multiple_plots.png", fig);