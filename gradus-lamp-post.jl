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

    data = lineprofile(m, x_obs, d, profile ;bins = g_domain, method = TransferFunctionMethod(), numrₑ = 50)
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

function invoke!(output, domain, model::CutoffPL)
    let Γ = model.Γ, β = model.β
        integration_kernel!(output, domain) do E, δE
            δE*E^(-1*Γ)*exp(-1*(E / β))
        end
    end
end

println("CutoffPl Loaded")

energies=collect(logrange(3, 79, 100))
model = CutoffPL()
output=invokemodel(energies,model)
plot(energies, output, xscale=:log10, yscale=:log10)

#= struct Cutoff{T} <: AbstractSpectralModel{T,Multiplicative}
    "Energy Cutoff"
    β::T
end

function Cutoff(;
    β = FitParam(100.0))
    Cutoff(;β)
end

function invoke!(output, domain, model::Cutoff)
    let β = model.β
        E = domain[1:(end-1)]
        output .= exp.(E./β)
    end
end
println("Cutoff Loaded")
##

 =#