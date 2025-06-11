# define our own compound model
# we want to model XS_Relconv()(XS_WarmAbsorber()*XillverD5())+XS_Relconv()(XillverD5())+PowerLaw()

using Plots
using SpectralFitting
using XSPECModels
using Relxill

function ISCO(a::Float64)
    Z_1 = 1.0+((1-a^2)^(1/3))*((1+a)^(1/3)+(1-a)^(1/3))
    Z_2 = (3*a^2+Z_1^2)^(1/2)
    if a >= 0 
        r = 3+Z_2-((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    else
        r = 3+Z_2+((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    end
end

struct DiscAbsModel{T} <: AbstractSpectralModel{T,Additive}
    # xillver
    "Normalisation."
    K::T
    Γ::T
    A_Fe::T
    logξ::T
    density::T
    θ::T
    # relconv
    index::T
    a::T
    r_abs::T
end

function DiscAbsModel(;
    # xillver
    K = FitParam(1.0e-3, frozen = false, lower_limit = 0.0, upper_limit = 1.0),
    Γ = FitParam(2.0, frozen = false, lower_limit = 1.0, upper_limit = 5.0),
    A_Fe = FitParam(1.0, frozen = true, lower_limit = 0.0, upper_limit = 10.0),
    logξ = FitParam(1.0, frozen = false, lower_limit = 0.0, upper_limit = 4.0),
    density = FitParam(17.0, frozen = true, lower_limit = 15.0, upper_limit = 19.0),
    θ = FitParam(30.0, frozen = false, lower_limit = 4.0, upper_limit = 86.0),
    # relconv
    index = FitParam(3.0, frozen = true, lower_limit = 0.0, upper_limit = 10.0),
    a = FitParam(0.0, frozen = 0.0, lower_limit = 0.0, upper_limit = 0.998),
    # note r_abs in units of r_ISCO
    r_abs = FitParam(3.0, frozen = false, lower_limit = 1.0, upper_limit = 10.0),
    )
    DiscAbsModel(K, Γ, A_Fe, logξ, density, θ, index, a, r_abs)
end

function SpectralFitting.invoke!(output, domain, model::DiscAbsModel)

    # extend domain so we can do the convolution
    # this will extend the domain from 10^-1 keV to 10^1.5 keV which should be fine
    Δ = 0.005
    our_low_bins = collect(-1.0:Δ:log10(domain[1])-Δ)
    our_low_bins = 10 .^ our_low_bins
    our_high_bins = collect(log10(domain[end])+Δ:Δ:1.5)
    our_high_bins = 10 .^ our_high_bins
    our_domain = vcat(our_low_bins, domain, our_high_bins)
    our_output = zeros(length(our_domain)-1)

    # setup first xillver model
    m1 = XillverD5(K = model.K, Γ = model.Γ, A_Fe = model.A_Fe, logXi = model.logξ, density = model.density, inclination = model.θ)

    # setup first relconv model
    m2_r_in = ISCO(model.a)
    m2_r_out = model.r_abs * ISCO(model.a)
    m2_r_break = 0.5*(m2_r_in + m2_r_out)
    m2 = XS_Relconv(model.index, model.index, m2_r_break, model.a, model.θ, m2_r_in, m2_r_out, 0.0)
    
    # evalueat first model which is relconv(xillver) for the inner disc
    invokemodel!(our_output, our_domain, m1)
    # convolution only works if there is something to convolve with (it crashes without the following condition)
    if maximum(our_output) > 0.0
        invokemodel!(our_output, our_domain, m2)
    end

    # return the result from the origin domain excluding the extended bins
    output .= our_output[length(our_low_bins)+1:length(our_low_bins)+length(output)]
end

energy = collect(range(1.0, 10.0, 150))
m = invokemodel(energy, DiscAbsModel())
plot(energy[1:end-1],m./diff(energy),xlim=(2.0, 8.0),xlabel="Energy (keV)")
