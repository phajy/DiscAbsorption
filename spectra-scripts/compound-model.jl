# define our own compound model

using Plots
using SpectralFitting
using XSPECModels
using Relxill

struct OurTestModel{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    incl::T
end

function OurTestModel(;
    K = FitParam(1.0),
    incl = FitParam(30.0, frozen=false, lower_limit=0.0, upper_limit=90.0)
    )
    OurTestModel(K, incl)
end

function SpectralFitting.invoke!(output, domain, model::OurTestModel)

    # extend domain so we can do the convolution
    # this will extend the domain from 10^-1 keV to 10^2 keV which should be fine
    Δ = 0.005
    our_low_bins = collect(-1.0:Δ:log10(domain[1])-Δ)
    our_low_bins = 10 .^ our_low_bins
    our_high_bins = collect(log10(domain[end])+Δ:Δ:2.0)
    our_high_bins = 10 .^ our_high_bins
    our_domain = vcat(our_low_bins, domain, our_high_bins)
    our_output = zeros(length(our_domain)-1)

    m1 = XillverD5(K = model.K, Γ = 2.0, A_Fe = 1.0, logXi = 1.0, density = 17.0, inclination = model.incl)

    Index1 = 3.0
    Index2 = 3.0
    r_br_g = 6.0
    Rin_ms = 1.0
    Rout_ms = 400.0
    m2 = XS_Relconv(Index1, Index2, r_br_g, 0.0, model.incl, Rin_ms, Rout_ms, 0.0)
    
    invokemodel!(our_output, our_domain, m1)
    invokemodel!(our_output, our_domain, m2)

    # return the result from the origin domain excluding the extended bins
    output .= our_output[length(our_low_bins)+1:length(our_low_bins)+length(output)]
end

energy = collect(range(1.0, 10.0, 150))
m = invokemodel(energy, OurTestModel())
plot(energy[1:end-1],m./diff(energy),xlim=(2.0, 8.0),xlabel="Energy (keV)")
