using Plots, SpectralFitting, CFITSIO, Base.Threads, Gradus

Threads.nthreads() = 1

#ring corona line

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
    r = FitParam(2.,lower_limit = 1.5, upper_limit = 10., frozen = false),
    h = FitParam(2.,lower_limit = 1.5, upper_limit = 50., frozen = false),
    Γ = FitParam(2.3,lower_limit = 1.0, upper_limit = 3., frozen = false),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    RingCoronaLineKerrz(K, r, h, Γ, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::RingCoronaLineKerrz)
    cur_dir = pwd()
    #kerrz = "/Users/er19801/kerrz/kerrz-0.1.12-65305f2f7efade22ec09597417524d8afab01676-aarch64-macos-none"
    kerrz = "/data/typhon2/DariusM/kerrz/zig-out/bin/kerrz"
    ID = Threads.threadid()
    emisivity_out_file = "Kerrz/Table/emsvty_g$(model.Γ)_h$(model.h)_r$(model.r).dat"
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

    #println("Writing Emissivity")
    #run(`$kerrz emissivity --velocity corotate --photon-index $(model.Γ) --nthreads $(Threads.nthreads()) --ring-like h:$(model.h),x:$(model.r) --output $cur_dir/$emisivity_out_file`)
    
    println("Writing LineProf")
    run(`$kerrz lineprof  --nradii 100 --nangles 200 --spin $(model.a) --incl $(model.θ) --ng $domain_size --rin $R_In --rout $R_Out --emissivity-profile  $cur_dir/$emisivity_out_file --output $cur_dir/$lineprof_out_file`)

    


    #rm("$cur_dir/$emisivity_out_file")
    
    lineprof = parse.(Float64,reduce(hcat, split.(readlines(lineprof_out_file),", ")))
    
    rm("$cur_dir/$lineprof_out_file")

    output .= lineprof[2,:][1:end-1]./(sum(lineprof[2,:][1:end-1].*lineprof[1,:][1:end-1]))
end

#ring corona full 

struct FullModelRingKerrz{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Ring Radius"
    r::T
    "Ring Height"
    h::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
    "Photon Index"
    Γ::T
    "Iron Abundance"
    A_Fe::T
    "Ionisation Parameter"
    logXi::T
    "density"
    density::T
end

function FullModelRingKerrz(;K = FitParam(1.0),
    r = FitParam(5.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    h = FitParam(5.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    R_in = FitParam(0.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(Inf, lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85, frozen = false),
    a = FitParam(0.7,lower_limit=0.0,upper_limit=0.998, frozen = false),
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 3., frozen = false),
    A_Fe = FitParam(1.0,lower_limit = 0.1, upper_limit = 100., frozen = true),
    logXi = FitParam(3.0,lower_limit= 2., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false))
    FullModelRingKerrz(K,r,h,R_in,R_out,θ,a,Γ,A_Fe,logXi,density)
end

function SpectralFitting.invoke!(output, domain, model::FullModelRingKerrz)
    convmodel = RingCoronaLineKerrz(
    K = FitParam(1.0),
    r = FitParam(model.r),
    h = FitParam(model.h),
    Γ = FitParam(model.Γ),
    R_in = FitParam(model.R_in),
    R_out = FitParam(model.R_out), 
    θ = FitParam(model.θ),
    a = FitParam(model.a))
    
    specmodel = XillverD5(
    K = FitParam(model.K),
    Γ = FitParam(model.Γ),
    A_Fe = FitParam(model.A_Fe),
    logXi = FitParam(model.logXi),
    density = FitParam(model.density), 
    inclination = FitParam(model.θ))
        
    convolution_model = AsConvolution(convmodel)
    Fmodel = convolution_model(specmodel)
    output .= invokemodel(domain,Fmodel)
end

#Lamp Post line

struct LPCoronaLineKerrz{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
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

function LPCoronaLineKerrz(;K = FitParam(1.0),
    h = FitParam(2.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 3., frozen = false),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    LPCoronaLineKerrz(K, h, Γ, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::LPCoronaLineKerrz)
    cur_dir = pwd()
    #kerrz = "/Users/er19801/kerrz/kerrz-0.1.12-65305f2f7efade22ec09597417524d8afab01676-aarch64-macos-none"
    kerrz = "kerrzcli"
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

    #println("$kerrz emissivity --velocity corotate --spin $(model.a) --photon-index $(model.Γ) --output-file $cur_dir/$emisivity_out_file --nphotons 500000 --nthreads $(Threads.nthreads()) --lamppost h:$(model.h),vr:0.0")
    
    #println("$kerrz lineprof  --nradii 100 --nangles 200 --spin $(model.a) --incl $(model.θ) --ng $domain_size --rin $R_In --rout $R_Out --emissivity-profile  $cur_dir/$emisivity_out_file --output $cur_dir/$lineprof_out_file")


    run(`$kerrz emissivity --velocity corotate --spin $(model.a) --photon-index $(model.Γ) --output-file $cur_dir/$emisivity_out_file --nphotons 3000 --nthreads $(Threads.nthreads()) --lamppost h:$(model.h),vr:0.0`)
    
    run(`$kerrz lineprof  --nradii 100 --nangles 200 --spin $(model.a) --incl $(model.θ) --ng $domain_size --rin $R_In --rout $R_Out --emissivity-profile  $cur_dir/$emisivity_out_file --output $cur_dir/$lineprof_out_file`)

    rm("$cur_dir/$emisivity_out_file")
    
    lineprof = parse.(Float64,reduce(hcat, split.(readlines(lineprof_out_file),", ")))
    
    rm("$cur_dir/$lineprof_out_file")

    output .= lineprof[2,:][1:end-1]./(sum(lineprof[2,:][1:end-1].*lineprof[1,:][1:end-1]))
end

#Lamp Post full

struct FullModelLPKerrz{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Corona Height"
    h::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
    "Photon Index"
    Γ::T
    "Iron Abundance"
    A_Fe::T
    "Ionisation Parameter"
    logXi::T
    "density"
    density::T
end

function FullModelLPKerrz(;K = FitParam(1.0),
    h = FitParam(5.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    R_in = FitParam(0.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(Inf, lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85, frozen = false),
    a = FitParam(0.7,lower_limit=0.0,upper_limit=0.998, frozen = false),
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 3., frozen = false),
    A_Fe = FitParam(1.0,lower_limit = 0.1, upper_limit = 100., frozen = true),
    logXi = FitParam(3.0,lower_limit= 2., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false))
    FullModelLPKerrz(K,h,R_in,R_out,θ,a,Γ,A_Fe,logXi,density)
end

function SpectralFitting.invoke!(output, domain, model::FullModelLPKerrz)
    convmodel = LPCoronaLineKerrz(
    K = FitParam(1.0),
    h = FitParam(model.h),
    Γ = FitParam(model.Γ),
    R_in = FitParam(model.R_in),
    R_out = FitParam(model.R_out), 
    θ = FitParam(model.θ),
    a = FitParam(model.a))
    
    specmodel = XillverD5(
    K = FitParam(model.K),
    Γ = FitParam(model.Γ),
    A_Fe = FitParam(model.A_Fe),
    logXi = FitParam(model.logXi),
    density = FitParam(model.density), 
    inclination = FitParam(model.θ))
        
    convolution_model = AsConvolution(convmodel)
    Fmodel = convolution_model(specmodel)
    output .= invokemodel(domain,Fmodel)
end

#= #plotting energies
line_energies = collect(range(0,2,1000))
spec_energies = collect(logrange(0.1,100,1000))

# Ring corona
@time begin
spec_Kerrz = invokemodel(spec_energies,FullModelRingKerrz(
    r = FitParam(1.5,lower_limit = 1.5, upper_limit = 10., frozen = false),
    h = FitParam(5.5,lower_limit = 1.5, upper_limit = 50., frozen = false),
    R_in = FitParam(-1.,lower_limit= 1. ,upper_limit=100, frozen = true),
    R_out = FitParam(400., lower_limit=400. ,upper_limit=600., frozen = true), 
    θ = FitParam(10.,lower_limit=10.,upper_limit=85., frozen = false),
    a = FitParam(0.1,lower_limit=0.0,upper_limit=0.998, frozen = false),
    Γ = FitParam(1.1,lower_limit = 1., upper_limit = 4., frozen = false),
    A_Fe = FitParam(1.5,lower_limit = 0.5, upper_limit = 10., frozen = false),
    logXi = FitParam(0.0,lower_limit= 0.0, upper_limit = 4.0,frozen = false),
    density = FitParam(15.0, lower_limit=15., upper_limit=19., frozen = false)
    ))
plot(spec_energies[1:end-1],spec_Kerrz,xscale=:log10,yscale=:log10,label="Kerrz",title="Ring Corona")
end =#