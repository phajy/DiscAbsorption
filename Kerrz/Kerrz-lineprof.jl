using Plots, SpectralFitting, XSPECModels, Relxill, CFITSIO, Base.Threads, Gradus
include("../gradus-lamp-post.jl")

cur_dir = pwd()
kerrz = "/Users/er19801/kerrz/kerrz-0.1.12-65305f2f7efade22ec09597417524d8afab01676-aarch64-macos-none"
Threads.nthreads() = 20

struct RingCoronaLineKerrz{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Ring Radius"
    r::T
    "Ring Height"
    h::T
    "Photon Index"
    Γ::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
end

function RingCoronaLineKerrz(;K = FitParam(1.0),
    r = FitParam(2.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    h = FitParam(2.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 3., frozen = false),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    RingCoronaLineKerrz(K, r, h, Γ, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::RingCoronaLineKerrz)
    ID = Threads.threadid()
    emisivity_out_file = "emisivity_$(ID)_temp.dat"
    lineprof_out_file = "lineprof_$(ID)_temp.dat"

    g_domain = copy(domain)
    domain_size = length(g_domain)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(KerrMetric(a = model.a))
    else
        R_In = model.R_in
    end 
    if model.R_out < 0 
        R_Out = abs(model.R_out) * Gradus.isco(KerrMetric(a = model.a))
    else
        R_Out = model.R_out
    end 

    run(`$kerrz emissivity --photon-index $(model.Γ) --nthreads $(Threads.nthreads()) --ring-like h:$(model.h),x:$(model.r) --output-file $cur_dir/$emisivity_out_file`)
    
    run(`$kerrz lineprof  --nradii 100 --nangles 200 --spin $(model.a) --incl $(model.θ) --ng $domain_size --rin $R_In --rout $R_Out --emissivity-profile  $cur_dir/$emisivity_out_file --output $cur_dir/$lineprof_out_file`)

    rm("$cur_dir/$emisivity_out_file")
    
    lineprof = parse.(Float64,reduce(hcat, split.(readlines(lineprof_out_file),", ")))
    
    rm("$cur_dir/$lineprof_out_file")

    output .= lineprof[2,:][1:end-1]./(sum(lineprof[2,:][1:end-1].*lineprof[1,:][1:end-1]))
end

energies = collect(range(0,2,1000))

@time begin
spec_Kerrz = invokemodel(energies,RingCoronaLineKerrz(;K = FitParam(1.0),
    r = FitParam(4.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    h = FitParam(6.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    Γ = FitParam(0.0,lower_limit = 1, upper_limit = 3., frozen = false),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(60.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998)))
plot(energies[1:end-1],spec_Kerrz,label="Kerrz")
end

@time begin
    spec_Gradus = invokemodel(energies,RingCoronaLine(;K = FitParam(1.0),
    r = FitParam(4.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    h = FitParam(6.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    E = FitParam(1.0,lower_limit = 1, upper_limit = 3., frozen = false),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(60.,lower_limit=7,upper_limit=85),
    K = FitParam(1.0),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998)))
    plot!(energies[1:end-1],spec_Gradus,label="Gradus")
end