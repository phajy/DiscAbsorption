using Gradus, Plots, SpectralFitting
struct Laor{T} <: AbstractSpectralModel{T,Additive}
    "Line energy"
    E::T
    "Emissivity index"
    a::T
    "Inner Radius"
    R_in::T
    "Outer Radius"
    R_out::T
    "Inclination"
    θ::T
    "Normalisation"
    K::T
end

# add a default keyword constructor
function Laor(;K = FitParam(1.0),
    E = FitParam(6.4,lower_limit = 1., upper_limit = 10., frozen = true),
    a = FitParam(3.,frozen=true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,upper_limit=90))
    Laor(E, a, R_in, R_out, θ, K)
end

function SpectralFitting.invoke!(output, domain, model::Laor)
    g_domain = copy(domain) ./ model.E
    
    m = KerrMetric(;a = 0.998)
    x_obs = SVector(0.0, 1e3, deg2rad(model.θ), 0.0)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    d = ThinDisc(R_In, model.R_out)
    emissivity(r) = r^-model.a
    data = lineprofile(g_domain, emissivity, m, x_obs, d, ; method = BinningMethod())
    output .= data[2][1:end-1]
end

model = Laor()

energies = collect(logrange(0.1,70,1000))
spec = invokemodel(energies,model)
plot(spec,energies[1:end-1])