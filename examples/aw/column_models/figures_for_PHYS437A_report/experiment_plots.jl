using OceanBioME, Oceananigans, Printf, CairoMakie, Statistics, NCDatasets, Interpolations
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

const year = years = 365days;
const month = months = year / 12;

# Create the MLD interpolator
mld_ds = NCDataset("../data/nc/mixed_layer_processed.nc")
mld_dat = convert(Array{Float64}, mld_ds["depth_mean"])
push!(mld_dat, mld_dat[1]) # to make looping work
mld_itp_months = interpolate(-mld_dat, BSpline(Cubic(Line(OnGrid()))))
mld_itp(t) = mld_itp_months(mod(t / month, 12) + 1)

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
        axN = Axis(fig[1, j-1]; title = "N for Experiment $(j-1)", axis_kwargs...)
        hmN = heatmap!(times / days, z, N[1,1,:,:]', colormap = :balance, colorrange = (2.5, 12.5))
        lines!(axN, times / days, mld_itp.(times), color = :black, linewidth = 5)
        axP = Axis(fig[2, j-1]; title = "P for Experiment $(j-1)", axis_kwargs...)
        hmP = heatmap!(times / days, z, P[1,1,:,:]', colormap = :balance, colorrange = (0, 0.4))
        lines!(axP, times / days, mld_itp.(times), color = :black, linewidth = 5)
        axZ = Axis(fig[3, j-1]; title = "Z for Experiment $(j-1)", axis_kwargs...)
        hmZ = heatmap!(times / days, z, Z[1,1,:,:]', colormap = :balance, colorrange = (0, 0.4))
        lines!(axZ, times / days, mld_itp.(times), color = :black, linewidth = 5)
        axD = Axis(fig[4, j-1]; title = "D for Experiment $(j-1)", axis_kwargs...)
        hmD = heatmap!(times / days, z, D[1,1,:,:]', colormap = :balance, colorrange = (0, 0.7))
        lines!(axD, times / days, mld_itp.(times), color = :black, linewidth = 5)
        if j == 4
            Colorbar(fig[1, 4], hmN)
            Colorbar(fig[2, 4], hmP)
            Colorbar(fig[3, 4], hmZ)
            Colorbar(fig[4, 4], hmD)
        end
    end
end
    
save("experiments_MLD_overlaid.png", fig);