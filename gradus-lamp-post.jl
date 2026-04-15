using Gradus, Plots, SpectralFitting

struct LampPost{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Corona Height"
    h::T
    "Line energy"
    E::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
end

function LampPost(;K = FitParam(1.0),
    h = FitParam(2.,lower_limit = 1.5, upper_limit = 10., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    LampPost(K, h, E, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::LampPost)
    g_domain = copy(domain)
    
    m = KerrMetric(;a = model.a)
    x_obs = SVector(0.0, 1e3, deg2rad(model.θ), 0.0)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    d = ThinDisc(R_In, model.R_out)

    mode = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, mode)

    data = lineprofile(m, x_obs, d, profile ;bins = g_domain, method = TransferFunctionMethod(), numrₑ = 10)
    output .= data[2][1:end-1]
end

println("LampPost Loaded")


struct CutoffPL{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Photon Index"
    Γ::T
    "Energy Cutoff"
    β::T
end

function CutoffPL(;
    K = FitParam(1.0),
    Γ = FitParam(2.0),
    β = FitParam(100.0))
    CutoffPL{typeof(K)}(K,Γ,β)
end

function SpectralFitting.invoke!(output, domain, model::CutoffPL)
    let Γ = model.Γ, β = model.β
        SpectralFitting.integration_kernel!(output, domain) do E, δE
            δE*E^(-1*Γ)*exp(-1*(E / β))
        end
    end
end

println("CutoffPl Loaded")

struct FullModel{T} <: AbstractSpectralModel{T,Additive}
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

function FullModel(;K = FitParam(1.0),
    h = FitParam(5.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    R_in = FitParam(0.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(Inf, lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85, frozen = false),
    a = FitParam(0.7,lower_limit=0.0,upper_limit=0.998, frozen = false),
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 3., frozen = false),
    A_Fe = FitParam(1.0,lower_limit = 0.1, upper_limit = 100., frozen = true),
    logXi = FitParam(3.0,lower_limit= 2., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false))
    FullModel(K,h,R_in,R_out,θ,a,Γ,A_Fe,logXi,density)
end

function SpectralFitting.invoke!(output, domain, model::FullModel)
    convmodel = LampPost(
    K = FitParam(1.0),
    h = FitParam(model.h),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
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

println("FullModel Loaded")
##
struct RingCoronaLine{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Ring Radius"
    r::T
    "Ring Height"
    h::T
    "Line energy"
    E::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Spin"
    a::T
end

function RingCoronaLine(;K = FitParam(1.0),
    r = FitParam(2.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    h = FitParam(2.,lower_limit = 1.5, upper_limit = 100., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(0.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    RingCoronaLine(K, r, h, E, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::RingCoronaLine)
    g_domain = copy(domain)
    
    m = KerrMetric(;a = model.a)
    x_obs = SVector(0.0, 1e3, deg2rad(model.θ), 0.0)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    d = ThinDisc(R_In, model.R_out)

    mode = RingCorona(r=model.r,h = model.h)
    profile = emissivity_profile(m, d, mode)

    data = lineprofile(m, x_obs, d, profile ;bins = g_domain, method = TransferFunctionMethod(), numrₑ = 30)
    output .= data[2][1:end-1]
end

println("RingCorona Loaded")

struct FullModelRing{T} <: AbstractSpectralModel{T,Additive}
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

function FullModelRing(;K = FitParam(1.0),
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
    FullModelRing(K,r,h,R_in,R_out,θ,a,Γ,A_Fe,logXi,density)
end

function SpectralFitting.invoke!(output, domain, model::FullModelRing)
    convmodel = RingCoronaLine(
    K = FitParam(1.0),
    r = FitParam(model.r),
    h = FitParam(model.h),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
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

println("FullModelRing Loaded")

energies = collect(logrange(2.5,90.0,900))
ring_spec = invokemodel(energies,FullModelRing(r=FitParam(1.5),h=FitParam(1.5)))
plot(energies[1:end-1],ring_spec,xscale=:log10,yscale=:log10)