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
    h = FitParam(2.,lower_limit = 1, upper_limit = 20., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998))
    LampPost(h, E, R_in, R_out, θ, a, K)
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

#numrₑ = number of transfer functions 
#plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 50.0)
using Colors
plot()
color = range(colorant"red", stop=colorant"blue", length=10)
Hrange = collect(logrange(1.5,10.,10))
energies = collect(logrange(0.1,1.5,500))
for i in eachindex(Hrange)
        model = LampPost(h=FitParam(Hrange[i])) 
        spec = invokemodel(energies,model)
        digs = round(Hrange[i],digits=2)
        global p = plot!(energies[1:end-1],spec,label = "h = $digs",color = color[i] )
end
display(p)
##

N=10
color = range(colorant"red", stop=colorant"blue", length=N)

specmodel = XillverD5()

convmodel = GaussianLine(μ = FitParam(1.))
#convmodel = LampPost()
energies = collect(logrange(0.1,50.,500))

plot(legend=false,energies[1:end-1],invokemodel(energies,specmodel)./sum(energies[1:end-1].*invokemodel(energies,specmodel)),yscale=:log10,xscale=:log10,color="black",label="unconvolved")

sigrange = collect(logrange(0.01,0.1,N))
for i in eachindex(sigrange)
    s = sigrange[i]
    convmodel.σ = s
    convolution_model = AsConvolution(convmodel)
    reflecmodel = convolution_model(XillverD5())
    spec = invokemodel(energies,reflecmodel)
    area = sum(energies[1:end-1].*spec)
    S = round(s,digits = ndigits(N)+1)
    global p = plot!(energies[1:end-1],spec./area,yscale=:log10,xscale=:log10,color=color[i],label = "σ = $S")
end
display(p)

