using Gradus, Plots, SpectralFitting, Warmabs, Colors, XSPECModels, CFITSIO, Relxill

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
    h = FitParam(2.,lower_limit = 1, upper_limit = 100., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    LampPost(K, h, E, R_in, R_out, θ, a)
end

function SpectralFitting.invoke!(output, domain, model::LampPost)
    g_domain = copy(domain)
    
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
println("LampPost Loaded")
#= convmodel = LampPost(    
    h = FitParam(1.5,lower_limit = 1, upper_limit = 20., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85, frozen = true),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998, frozen = true))

specmodel = XillverD5(
    Γ = FitParam(1.5,lower_limit = 1, upper_limit = 20., frozen = false),
    A_Fe = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = false),
    logXi = FitParam(3.,lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false), 
    inclination = FitParam(30.,lower_limit=7,upper_limit=85, frozen = true))

N=10
color = range(colorant"red", stop=colorant"blue", length=N)
energies = collect(logrange(0.1,100.0,1000))
hrange = collect(logrange(1.5,50.0,N))
plot(xlims=(0,1.5))
for i in eachindex(hrange)
    s = hrange[i]
    convmodel.h = s
    spec = invokemodel(energies,convmodel)
    S = round(s,digits = ndigits(N))
    global p = plot!(energies[1:end-1],spec,color=color[i],label = "h = $S")
end
display(p)

unconvspec = invokemodel(energies,specmodel)
area = sum(energies[1:end-1].*unconvspec)
plot(energies[1:end-1],unconvspec./area,yscale=:log10,xscale=:log10,xlims=(0.5,70),xlabel="Energy (keV)",color="black",legend=false)

for i in eachindex(hrange)
    s = hrange[i]
    convmodel.h = s
    convolution_model = AsConvolution(convmodel)
    reflecmodel = convolution_model(specmodel)
    spec = invokemodel(energies,reflecmodel)
    area = sum(energies[1:end-1].*spec)
    S = round(s,digits = ndigits(N)+1)
    global p = plot!(energies[1:end-1],spec./area,yscale=:log10,xscale=:log10,color=color[i],label = "h = $S")
end
display(p)
plot!(legend=false) =#