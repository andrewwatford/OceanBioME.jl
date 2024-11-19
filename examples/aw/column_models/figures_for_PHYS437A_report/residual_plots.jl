using Oceananigans, CairoMakie
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

files = ["../data_files_for_PHYS437A_report/long_00.jld2",
        "../data_files_for_PHYS437A_report/long_01.jld2", 
        "../data_files_for_PHYS437A_report/long_10.jld2", 
        "../data_files_for_PHYS437A_report/long_11.jld2"]
cases = ["Control (T = 12ᵒC)", "T-dependent Mortality", "T-dependent Growth", "All T-dependent"]
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

fig = Figure(size = (2560, 2304), fontsize = 50)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((365.0, 730.0), (z[1], 0)))
for (j, data_tuple) = enumerate(data)
    if j != 1
        N = data_tuple[:N]
        P = data_tuple[:P]
        Z = data_tuple[:Z]
        D = data_tuple[:D]
        # Plot the residuals
        label = cases[j]
        axN = Axis(fig[1, j-1]; title = "ΔN for Experiment $(j-1)", axis_kwargs...)
        hmN = heatmap!(times / days, z, N[1,1,:,:]' .- controlN[1,1,:,:]', colormap = :balance, colorrange = (-0.75, 0.75))
        axP = Axis(fig[2, j-1]; title = "ΔP for Experiment $(j-1)", axis_kwargs...)
        hmP = heatmap!(times / days, z, P[1,1,:,:]' .- controlP[1,1,:,:]', colormap = :balance, colorrange = (-0.07, 0.07))
        axZ = Axis(fig[3, j-1]; title = "ΔZ for Experiment $(j-1)", axis_kwargs...)
        hmZ = heatmap!(times / days, z, Z[1,1,:,:]' .- controlZ[1,1,:,:]', colormap = :balance, colorrange = (-0.05, 0.05))
        axD = Axis(fig[4, j-1]; title = "ΔD for Experiment $(j-1)", axis_kwargs...)
        hmD = heatmap!(times / days, z, D[1,1,:,:]' .- controlD[1,1,:,:]', colormap = :balance, colorrange = (-0.2, 0.2))
        if j == 4
            Colorbar(fig[1, 4], hmN)
            Colorbar(fig[2, 4], hmP)
            Colorbar(fig[3, 4], hmZ)
            Colorbar(fig[4, 4], hmD)
        end
    end
end
    
save("residuals.png", fig);