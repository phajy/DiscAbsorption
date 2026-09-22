using Plots

a = 0.998
θ = 70.
r = 5.0
h = 10.
Γ = 2.0

include("Kerrz_Lineprof_Ccall.jl")


metric = make_metric(1.0, a)
x_obs = KrzFourVector(0.0, 1e6, deg2rad(θ), 0.0)
ring  = KrzRingCorona(h, r)   # height, radius 

r_in  = metric.isco + 0.1
r_out = 400.0

r_grid = logspace(r_in, r_out, 1000)
g_grid = collect(range(0.0, 2.0; length = 1000))  # 100 bins

@time begin
    flux = build_lineprofile(metric, x_obs, ring, r_grid, g_grid,
                            emissivity_num_traces = 1000,
                            tf_max_points = 100,
                            n_threads = 4,
                            gamma_index = Γ,      # photon index Γ of illuminating spectrum
                            beaming_exponent = 3.0)
    
    dg = sum(diff(g_grid))/length(diff(g_grid))
    total = sum(flux) * dg
    flux ./= total
    g_centers = (g_grid[1:end-1] .+ g_grid[2:end]) ./ 2
    
    plot(g_centers,flux,xlabel="energy",ylabel="flux",label="C-call")
end

include("KerrRingLinetest.jl")

@time begin
    model_kerrz = RingCoronaLineKerrz_mod(;K = FitParam(1.0),
    r = FitParam(r),
    h = FitParam(h),
    Γ = FitParam(Γ),
    R_in = FitParam(r_in),
    R_out = FitParam(r_out), 
    θ = FitParam(θ),
    a = FitParam(a)
    
    )
    
    spec_kerrz = invokemodel(g_grid, model_kerrz)
    flux_kerrz = spec_kerrz.parent[:, 1]
    total_kerrz = sum(flux_kerrz) * dg
    flux_kerrz ./= total_kerrz
    plot!(g_centers,flux_kerrz,label="Kerrz File")
end
plot!()
##
@time begin
    
    g_grid_ = collect(logrange(0.1, 1000.0; length = 1000))  # 100 bins
    
    
    model_kerrz = FullModelRingKerrz_mod(;K = FitParam(1.0),
    r = FitParam(r),
    h = FitParam(7.0),
    Γ = FitParam(Γ),
    R_in = FitParam(r_in),
    R_out = FitParam(r_out), 
    θ = FitParam(θ),
    a = FitParam(a)
    
    )
    
    spec_kerrz = invokemodel(g_grid_, model_kerrz)

    plot(g_grid_[1:end-1],spec_kerrz.parent[:,1],xlims=(1.0,100.0),xscale=:log10,yscale=:log10, legend=false)
    
end
    