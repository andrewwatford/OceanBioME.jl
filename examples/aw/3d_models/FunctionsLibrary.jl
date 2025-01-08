module FunctionsLibrary
    export τ_surface, T_surface, S_surface
    export plot_forcing_profiles, plot_final_profiles, plot_final_profiles_NPZD
end

using NCDatasets

@inline function τ_surface(i, j, grid, clock, fields, p)
        t = clock.time
    τ_amp = p.τ_amp_winter + (p.τ_amp_summer - p.τ_amp_winter) * abs(sin(π* t/p.T_seasonal))
    return -  τ_amp * cos(p.τ_period * π * (grid.φᵃᶜᵃ[j] - φ0) / p.Lφ + p.τ_offset)
end

@inline function T_surface(i, j, grid, clock, fields, p)
          t = clock.time
    sst_min = p.sst_min_winter + (p.sst_min_summer - p.sst_min_winter) * abs(sin(π* t/p.T_seasonal))
    return p.sst_amp * cos(p.sst_period*π * (grid.φᵃᶜᵃ[j] - φ0) / p.Lφ + p.sst_offset) + sst_min
end

@inline function S_surface(i, j, grid, clock, fields, p)
          t = clock.time
    sal_amp =  p.sal_amp_winter + (p.sal_amp_summer - p.sal_amp_winter) * abs(sin(π* t/p.T_seasonal))
    return - sal_amp * cos(p.sal_period*π * (grid.φᵃᶜᵃ[j] - φ0) / p.Lφ + p.sal_offset)
end

@inline function ζ(u, v, Δx, Δy)

    vx = (v[2:end, 1:end] .- v[1:end-1,1:end] )./Δx    # Nλ-1, Nφ+1
    uy = (u[1:end, 2:end] .- u[1:end, 1:end-1])./Δy    # Nλ+1, Nφ-1
    
    avevx = (vx[1:end-1, 1:end-1] .+ vx[1:end-1, 2:end] .+ vx[2:end, 1:end-1] .+ vx[2:end, 2:end])/4  # Nλ-2, Nφ
    aveuy = (uy[1:end-1, 1:end-1] .+ uy[1:end-1, 2:end] .+ uy[2:end, 1:end-1] .+ uy[2:end, 2:end])/4  # Nλ,   Nφ-2
    
    ζz = avevx[1:end,2:end-1] .- aveuy[2:end-1,1:end]
    
    return ζz
end

function plot_forcing_profiles(grid, model, parameters, τ_surface, T_surface, S_surface)

    Nφ  = grid.Ny
    lat = zeros(Nφ)
    winds_summer, sst_summer, sal_summer = zeros(Nφ), zeros(Nφ), zeros(Nφ)
    winds_winter, sst_winter, sal_winter = zeros(Nφ), zeros(Nφ), zeros(Nφ)

    model.clock.time = 0.
    for j = 1:Nφ
        lat[j] = grid.φᵃᶜᵃ[j]
        winds_winter[j] = τ_surface(1, j, grid, model.clock, 0, parameters)
        sst_winter[j]   = T_surface(1, j, grid, model.clock, 0, parameters)
        sal_winter[j]   = S_surface(1, j, grid, model.clock, 0, parameters)
    end

    model.clock.time = 180*days
    for j = 1:Nφ
        winds_summer[j] = τ_surface(1, j, grid, model.clock, 0, parameters)
        sst_summer[j]   = T_surface(1, j, grid, model.clock, 0, parameters)
        sal_summer[j]   = S_surface(1, j, grid, model.clock, 0, parameters)
    end

    model.clock.time = 0.

    fig = Figure(size = (1400, 700))

    ax_winds = Axis(fig[1, 1];  title = "winds")
    ax_sst   = Axis(fig[1, 2];  title = "SST")
    ax_sal   = Axis(fig[1, 3];  title = "salinity")

    lines!(ax_winds, winds_winter, lat, color=:blue, label="winter")
    lines!(ax_sst,   sst_winter,   lat, color=:blue, label="winter")
    lines!(ax_sal,   sal_winter,   lat, color=:blue, label="winter")
    lines!(ax_winds, winds_summer, lat, color=:red,  label="summer")
    lines!(ax_sst,   sst_summer,   lat, color=:red,  label="summer")
    lines!(ax_sal,   sal_summer,   lat, color=:red,  label="summer")

    axislegend(ax_winds)
    axislegend(ax_sst)
    axislegend(ax_sal)
    save("winds_SST_salinity.png", fig)

end

function plot_final_profiles(grid, index)

    Nλ, Nφ = grid.Nx, grid.Ny

    λc = λnodes(grid, Center())
    φc = φnodes(grid, Center())
    λf = λnodes(grid, Face())
    φf = φnodes(grid, Face())
    z = znodes(grid, Center())
    
    ds = NCDataset(simulation.output_writers[:fields].filepath, "r")
    
    fig = Figure(size = (800, 600))
    
    axis_kwargs_h = (xlabel = "longitude",
                     ylabel = "latitude",
                   limits = ((λf[1], λf[Nλ+1]), (φf[1], φf[Nφ+1])),)
    
    axis_kwargs_v = (xlabel = "latitude",
                   ylabel = "depth",
                   limits = ((φc[1], φc[Nφ]), (-Lz, 0)),)
                       
    Tfinal_v = ds["T"][30,     1:Nφ,   1:Nz, index]
    Sfinal_v = ds["S"][30,     1:Nφ,   1:Nz, index]
    ufinal_h = ds["u"][1:Nλ+1, 1:Nφ,     Nz, index]
    vfinal_h = ds["v"][1:Nλ,   1:Nφ+1,   Nz, index]
    ufinal_v = ds["u"][30,     1:Nφ,   1:Nz, index]
    vfinal_v = ds["v"][30,     1:Nφ+1, 1:Nz, index]
    
    close(ds)
    
    ax_T = Axis(fig[2,1]; title = "temperature", axis_kwargs_v...)               
    T_heatmap = heatmap!(ax_T, φc, z, Tfinal_v, colormap=:heat)
    Colorbar(fig[2,2], T_heatmap; label = "C")
    
    ax_S = Axis(fig[2,3]; title = "salinity", axis_kwargs_v...)               
    S_heatmap = heatmap!(ax_S, φc, z, Sfinal_v, colormap=:heat)
    Colorbar(fig[2,4], S_heatmap; label = "S")
    
    ax_u = Axis(fig[3,1]; title = "x-velocity", axis_kwargs_h...)               
    u_heatmap = heatmap!(ax_u, λf, φc, ufinal_h, colormap=:balance)
    Colorbar(fig[3,2], u_heatmap; label = "m/s")
    
    ax_v = Axis(fig[3,3]; title = "y-velocity", axis_kwargs_h...)               
    v_heatmap = heatmap!(ax_v, λc, φf, vfinal_h, colormap=:balance)
    Colorbar(fig[3,4], v_heatmap; label = "m/s")
    
    ax_u_v = Axis(fig[4,1]; title = "x-velocity (vertical)", axis_kwargs_v...)               
    u_heatmap_v = heatmap!(ax_u_v, φc, z, ufinal_v, colormap=:balance)
    Colorbar(fig[4,2], u_heatmap_v; label = "m/s")
    
    ax_v_v = Axis(fig[4,3]; title = "y-velocity (vertical)", axis_kwargs_v...)               
    v_heatmap_v = heatmap!(ax_v_v, φf, z, vfinal_v, colormap=:balance)
    Colorbar(fig[4,4], v_heatmap_v; label = "m/s")
    
    title = "Vertical profiles of T and S"
    fig[1, :] = Label(fig, title, fontsize=20, tellwidth=false)
    
    save("T_S_u_v_final.png", fig)
end

function plot_final_profiles_NPZD(grid, index)

    Nλ, Nφ = grid.Nx, grid.Ny

    λc = λnodes(grid, Center())
    φc = φnodes(grid, Center())
    λf = λnodes(grid, Face())
    φf = φnodes(grid, Face())
    z = znodes(grid, Center())
    
    ds = NCDataset(simulation.output_writers[:fields].filepath, "r")
    
    fig = Figure(size = (800, 600))
    
    axis_kwargs_h = (xlabel = "longitude",
                     ylabel = "latitude",
                     limits = ((λf[1], λf[Nλ+1]), (φf[1], φf[Nφ+1])),)
    
    axis_kwargs_v = (xlabel = "latitude",
                     ylabel = "depth",
                     limits = ((φc[1], φc[Nφ]), (-Lz, 0)),)
                       
    N_v = ds["N"][Int(Nλ/2), 1:Nφ,   1:Nz, index]
    P_v = ds["P"][Int(Nλ/2), 1:Nφ,   1:Nz, index]
    Z_v = ds["Z"][Int(Nλ/2), 1:Nφ,   1:Nz, index]
    D_v = ds["D"][Int(Nλ/2), 1:Nφ, 1:Nz, index]
    N_h = ds["N"][1:Nλ,      1:Nφ,     Nz, index]
    P_h = ds["P"][1:Nλ,      1:Nφ,   Nz, index]
    Z_h = ds["Z"][1:Nλ,      1:Nφ,     Nz, index]
    D_h = ds["D"][1:Nλ,      1:Nφ,   Nz, index]
    
    close(ds)
    
    ax_N = Axis(fig[2,1]; title = "Nutrient (v)", axis_kwargs_v...)               
    N_hm = heatmap!(ax_N, φc, z, N_v, colormap=:heat)
    Colorbar(fig[2,2], N_hm; label = "")
    
    ax_P = Axis(fig[2,3]; title = "Phytoplankton (v)", axis_kwargs_v...)               
    P_hm = heatmap!(ax_P, φc, z, P_v, colormap=:heat)
    Colorbar(fig[2,4], P_hm; label = "")
    
    ax_Z = Axis(fig[3,1]; title = "Zooplankton (v)", axis_kwargs_v...)               
    Z_hm = heatmap!(ax_Z, φc, z, Z_v, colormap=:heat)
    Colorbar(fig[3,2], Z_hm; label = "")
    
    ax_D = Axis(fig[3,3]; title = "Detritus (v)", axis_kwargs_v...)               
    D_hm = heatmap!(ax_D, φc, z, D_v, colormap=:heat)
    Colorbar(fig[3,4], D_hm; label = "")

    ax_N_h = Axis(fig[4,1]; title = "Nutrient (h)", axis_kwargs_h...)               
    N_hm_h = heatmap!(ax_N_h, λf, φc, N_h, colormap=:heat)
    Colorbar(fig[4,2], N_hm_h; label = "")
    
    ax_P_h = Axis(fig[4,3]; title = "Phytoplankton (h)", axis_kwargs_h...)               
    P_hm_h = heatmap!(ax_P_h, λc, φf, P_h, colormap=:heat)
    Colorbar(fig[4,4], P_hm_h; label = "")
    
    ax_Z_h = Axis(fig[5,1]; title = "Zooplankton (h)", axis_kwargs_h...)               
    Z_hm_h = heatmap!(ax_Z_h, λf, φc, Z_h, colormap=:heat)
    Colorbar(fig[5,2], Z_hm_h; label = "")
    
    ax_D_h = Axis(fig[5,3]; title = "Detritus (h)", axis_kwargs_h...)               
    D_hm_h = heatmap!(ax_D_h, λc, φf, D_h, colormap=:heat)
    Colorbar(fig[5,4], D_hm_h; label = "")
     
    title = "Vertical and Horizontal profiles of N, P, Z, D"
    fig[1, :] = Label(fig, title, fontsize=20, tellwidth=false)
    
    save("NPZD_final.png", fig)
end

function plot_streamlines(grid)

    Nλ, Nφ = grid.Nx, grid.Ny

    Δλ = grid.Δλᶠᵃᵃ
    Δφ = grid.Δφᵃᶠᵃ
    # Δz = grid.Δzᵃᵃᶠ

    ds = NCDataset(simulation.output_writers[:fields].filepath, "r")

    λc = λnodes(grid, Center())
    φc = φnodes(grid, Center())
    λf = λnodes(grid, Face())
    φf = φnodes(grid, Face())
    z = znodes(grid, Center())

    times = ds["time"][:]
    
    fig = Figure(size = (800, 800))

    axis_kwargs_h = (xlabel = "longitude",
                     ylabel = "latitude",
                     limits = ((λf[1], λf[Nλ+1]), (φf[1], φf[Nφ+1])),)

    n = Observable(1)

    u = @lift ds["u"][: , :, :, $n]
    v = @lift ds["v"][: , :, :, $n]
    
    U = @lift mean($u, dims=3)                         
    V = @lift mean($v, dims=3)                      
    ω = @lift ζ($U, $V, Δλ, Δφ)

    psi =  zeros(Nλ, Nφ)

    on(n) do dummy_var
        psi =  zeros(Nλ, Nφ)
    end

    psi_partial = @lift Δφ ./ 2 .* ($U[:, 2:end] .+ $U[:, 1:end-1])
    psi         = @lift cumsum($psi_partial, dims=2)
 
    psi = @lift $psi[1:end-1,:]
    max_ω = @lift max(maximum($ω), 1e-3)
    lim_ω =  @lift [- $max_ω ,$max_ω]
    φp = φc[1:end-1]

    title =  @lift @sprintf("BT ω and streamfunction @ t = %.2f days", times[$n]/(3600*24))
    fig[1, 1:2] = Label(fig, title, fontsize=24, tellwidth=false)

    ax_ω = Axis(fig[2,1]; title = "ω", axis_kwargs_h...)
    ω_heatmap = heatmap!(ax_ω,λc, φc, ω, colormap=:curl, colorrange=lim_ω)
    Colorbar(fig[2,2], ω_heatmap; label = L"s^{-1}", labelsize = 20)
    ct_psi = contour!(λc, φp, psi, color = :black, linewidth = 2, levels=12)

    frames = 1:length(times)
    video1=VideoStream(fig,format="mp4", framerate=12)

    for i=1:frames[end]
        recordframe!(video1)
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[]=i
    end
    save("ZetaPsi.mp4", video1)
    close(ds)
end

function plot_velocity(grid)
    Nλ, Nφ = grid.Nx, grid.Ny
    Δλ = grid.Δλᶠᵃᵃ
    ds = NCDataset(simulation.output_writers[:fields].filepath, "r")
    
    λc, φc = λnodes(grid, Center()), φnodes(grid, Center())
    λf, φf = λnodes(grid, Face()),   φnodes(grid, Face())
    z = znodes(grid, Center())
    
    times = ds["time"][:]
        
    fig = Figure(size = (800, 800))
    
    axis_kwargs_h = (xlabel = "longitude",
                     ylabel = "latitude",
                     limits = ((λc[1], λc[Nλ]), (φc[1], φc[Nφ])),)
    
    n = Observable(1)
    u = @lift ds["u"][1:Nλ , 1:Nφ, end-3, $n]
    v = @lift ds["v"][1:Nλ , 1:Nφ, end-3, $n]
    T = @lift ds["T"][1:Nλ , 1:Nφ, end-3, $n]
    
    title =  @lift @sprintf("t = %.2f days", times[$n]/(3600*24))
    fig[1, 1:2] = Label(fig, title, fontsize=24, tellwidth=false)
    
    ax_u = Axis(fig[2,1]; title = "zonal velocity", axis_kwargs_h...)
    u_heatmap = heatmap!(ax_u, λf, φc, u, colorrange = (-0.2, 0.2), colormap=:balance)
    Colorbar(fig[2,2], u_heatmap)
    
    ax_v = Axis(fig[3,1]; title = "meridional velocity", axis_kwargs_h...)
    v_heatmap = heatmap!(ax_v, λc, φf, v, colorrange = (-0.1, 0.1), colormap=:balance)
    Colorbar(fig[3,2], v_heatmap)

    frames = 1:length(times)
    video1=VideoStream(fig,format="mp4", framerate=12)
    
    for i=1:frames[end]
        recordframe!(video1)
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[]=i
    end
    close(ds)
    save("velocities.mp4", video1)
    
end
    
function plot_velocity_BT(grid)
    Nλ, Nφ = grid.Nx, grid.Ny
    Δλ = grid.Δλᶠᵃᵃ
    ds = NCDataset(simulation.output_writers[:fields].filepath, "r")
    
    λc, φc = λnodes(grid, Center()), φnodes(grid, Center())
    λf, φf = λnodes(grid, Face()),   φnodes(grid, Face())
    z = znodes(grid, Center())

    times = ds["time"][:]
        
    Z = zeros(Float32, Nλ, Nφ, Nz)
    for i=1:Nλ, j=1:Nφ
        Z[i,j,:] = z
    end

    fig = Figure(size = (800, 800))
    
    axis_kwargs_h = (xlabel = "longitude",
                     ylabel = "latitude",
                     limits = ((λc[1], λc[Nλ]), (φc[1], φc[Nφ])),)
    
    n = Observable(1)
    u = @lift ds["u"][1:Nλ , 1:Nφ, 1:Nz, $n]
    v = @lift ds["v"][1:Nλ , 1:Nφ, 1:Nz, $n]
    
    on(n) do dummy_var
        u_BT = zeros(Float32, Nλ, Nφ)
        v_BT = zeros(Float32, Nλ, Nφ)
    end

    u_BT = @lift sum(Z.*$u, dims=3)[:,:,1]/grid.Lz
    v_BT = @lift sum(Z.*$v, dims=3)[:,:,1]/grid.Lz

    title =  @lift @sprintf("t = %.2f days", times[$n]/(3600*24))
    fig[1, 1:2] = Label(fig, title, fontsize=24, tellwidth=false)
    
    ax_u = Axis(fig[2,1]; title = "zonal BT velocity", axis_kwargs_h...)
    u_heatmap = heatmap!(ax_u, λf, φc, u_BT, colorrange = (-0.2, 0.2), colormap=:balance)
    Colorbar(fig[2,2], u_heatmap)
    
    ax_v = Axis(fig[3,1]; title = "meridional BT velocity", axis_kwargs_h...)
    v_heatmap = heatmap!(ax_v, λc, φf, v_BT, colorrange = (-0.1, 0.1), colormap=:balance)
    Colorbar(fig[3,2], v_heatmap)

    frames = 1:length(times)
    video1=VideoStream(fig,format="mp4", framerate=12)
    
    for i=1:frames[end]
        recordframe!(video1)
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[]=i
    end
    close(ds)
    save("velocities_BT.mp4", video1)
    
end
    