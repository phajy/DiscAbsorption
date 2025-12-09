using Gradus, Plots, SpectralFitting

struct LampPost{T} <: AbstractSpectralModel{T,Additive}
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
    "Normalisation"
    K::T
end

# add a default keyword constructor
function LampPost(;K = FitParam(1.0),
    h = FitParam(2.,lower_limit = 2.0, upper_limit = 20., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    LampPost(h, E, R_in, R_out, θ, a, K)
end

function SpectralFitting.invoke!(output, domain, model::LampPost)
    g_domain = copy(domain) ./ model.E
    
    m = KerrMetric(;a = 0.998)
    x_obs = SVector(0.0, 1e3, deg2rad(model.θ), 0.0)

    if model.R_in < 0 
        R_In = abs(model.R_in) * Gradus.isco(m)
    else
        R_In = model.R_in
    end 

    d = ThinDisc(R_In, model.R_out)

    mode = LampPostModel(h = model.h)
    profile = emissivity_profile(m, d, mode)

    data = lineprofile(m, x_obs, d, profile ;bins = g_domain, method = TransferFunctionMethod(), numrₑ = 100)
    output .= data[2][1:end-1]
end

#numrₑ = number of transfer functions 
#plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 50.0)

model = LampPost() 


energies = collect(logrange(1.,80.,500))
spec = invokemodel(energies,model)
plot(energies[1:end-1],spec) 

line_convolution = AsConvolution(GaussianLine())
