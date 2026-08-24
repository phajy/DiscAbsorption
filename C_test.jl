using Plots

a = 0.998
θ = 30.
r = 5.0
h = 10.
Γ = 2.0

include("Kerrz_Lineprof_Ccall.jl")

@time begin
    # Example use
    
    metric = make_metric(1.0, a)
    x_obs = KrzFourVector(0.0, 1e6, deg2rad(θ), 0.0)
    ring  = KrzRingCorona(h, r)   # height, radius 
    
    r_in  = metric.isco + 0.1
    r_out = 400.0
    
    r_grid = logspace(r_in, r_out, 1000)
    g_grid = collect(range(0.0, 2.0; length = 501))  # 100 bins
    
    flux = build_lineprofile(metric, x_obs, ring, r_grid, g_grid,
    emissivity_num_traces = 50000,
    tf_max_points = 100,
    n_threads = 16,
    gamma_index = Γ,      # photon index Γ of illuminating spectrum
    beaming_exponent = 3.0)
    
    dg = sum(diff(g_grid))/length(diff(g_grid))
    total = sum(flux) * dg
    flux ./= total
    g_centers = (g_grid[1:end-1] .+ g_grid[2:end]) ./ 2
    
    plot(g_centers,flux,xlabel="energy",ylabel="flux",label="C-call")
end

include("Kerrz/Kerrz-lineprof.jl")

@time begin
    model_kerrz = RingCoronaLineKerrz(;K = FitParam(1.0),
    r = FitParam(r),
    h = FitParam(h),
    Γ = FitParam(Γ),
    R_in = FitParam(r_in),
    R_out = FitParam(r_out), 
    θ = FitParam(θ),
    a = FitParam(a))
    
    spec_kerrz = invokemodel(g_grid, model_kerrz)
    flux_kerrz = spec_kerrz.parent[:, 1]
    total_kerrz = sum(flux_kerrz) * dg
    flux_kerrz ./= total_kerrz
    plot!(g_centers,flux_kerrz,label="Kerrz File")
end
    
plot!()
    