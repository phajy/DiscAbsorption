### A Pluto.jl notebook ###
# v0.20.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ a549e35f-3151-4043-ab57-7029119a13d5
begin
	# make sure we have the necessary packages in Pluto environment
	import Pkg
	Pkg.add("Plots")
	Pkg.add("SpectralFitting")
	Pkg.add("XSPECModels")
	Pkg.add("Relxill")
	Pkg.add("Warmabs")
	Pkg.add("PlutoUI")
end

# ╔═╡ 345c379e-5092-46e7-adf7-decbbfcce731
using Plots, SpectralFitting, XSPECModels, Relxill, Warmabs, PlutoUI

# ╔═╡ 19887774-4c34-11f0-27ad-27b1f250b3d5
# Pluto notebook to interactively explore disc absorption models

# ╔═╡ a6104fbc-09ea-4ef6-a915-25c10b9c7f48
function ISCO(a::Float64)
    Z_1 = 1.0+((1-a^2)^(1/3))*((1+a)^(1/3)+(1-a)^(1/3))
    Z_2 = (3*a^2+Z_1^2)^(1/2)
    if a >= 0 
        r = 3+Z_2-((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    else
        r = 3+Z_2+((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    end
end

# ╔═╡ 4ba2e218-81a5-4e14-bb63-8977142b5ee4
begin
	struct DiscAbsModel{T} <: AbstractSpectralModel{T,Additive}
	    # xillver
	    "Normalisation."
	    K::T
	    K_i::T
	    Γ::T
	    A_Fe::T
	    logξ::T
	    density::T
	    θ::T
	    # relconv
	    index::T
	    a::T
	    r_in::T
	    r_abs::T
	    r_out::T
	    limb::T
	    # outer disc
	    K_o::T
	    # warm absorber
	    column::T
	    abs_logξ::T
	    vturb::T
	    z::T
	    # power law
	    K_pl::T
	end
	
	function DiscAbsModel(;
	    K = FitParam(1.0, frozen = true, lower_limit = 0.0, upper_limit = 1.0),
	    # xillver (inner disc)
	    K_i = FitParam(6.4e-12, frozen = false, lower_limit = 0.0, upper_limit = 1.0e-10),
	    Γ = FitParam(2.58, frozen = false, lower_limit = 1.0, upper_limit = 5.0),
	    A_Fe = FitParam(3.5, frozen = true, lower_limit = 0.0, upper_limit = 10.0),
	    logξ = FitParam(1.0, frozen = false, lower_limit = 0.0, upper_limit = 4.0),
	    density = FitParam(17.0, frozen = true, lower_limit = 15.0, upper_limit = 19.0),
	    θ = FitParam(70.0, frozen = false, lower_limit = 4.0, upper_limit = 86.0),
	    # relconv (inner disc)
	    index = FitParam(3.0, frozen = true, lower_limit = 0.0, upper_limit = 10.0),
	    a = FitParam(0.998, frozen = true, lower_limit = 0.0, upper_limit = 0.998),
	    # note r_in and r_abs in units of r_ISCO
	    r_in = FitParam(1.0, frozen = true, lower_limit = 1.0, upper_limit = 10.0),
	    r_abs = FitParam(3.0, frozen = false, lower_limit = 1.0, upper_limit = 10.0),
	    r_out = FitParam(400.0, frozen = true, lower_limit = 200.0, upper_limit = 800.0),
	    limb = FitParam(0.0, frozen = true, lower_limit = 0.0, upper_limit = 1.0),
	    # xillver (outer disc)
	    K_o = FitParam(9.6e-12, frozen = false, lower_limit = 0.0, upper_limit = 1.0e-10),
	    # warm absorber
	    column = FitParam(0.0, frozen = false, lower_limit = -3.0, upper_limit = 2.0, error = 0.1),
	    abs_logξ = FitParam(3.0, frozen = false, lower_limit = -4.0, upper_limit = 5.0),
	    vturb = FitParam(200.0, frozen = true, lower_limit = 0.0, upper_limit = 1000.0),
	    z = FitParam(0.066, frozen = true, lower_limit = 0.0, upper_limit = 1.0),
	    # power law
	    K_pl = FitParam(1.3e-7, frozen = false, lower_limit = 0.0, upper_limit = 1.0e-5)
	    )
	    DiscAbsModel(K, K_i, Γ, A_Fe, logξ, density, θ, index, a, r_in, r_abs, r_out, limb, K_o, column, abs_logξ, vturb, z, K_pl)
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
	
	    # diagnositcs
	    @info "Invoking DiscAbsModel with parameters:", 
	          "K = ", model.K,
	          " K_i = ", model.K_i,
	          " Γ = ", model.Γ,
	          " A_Fe = ", model.A_Fe,
	          " logξ = ", model.logξ,
	          " density = ", model.density,
	          " θ = ", model.θ,
	          " index = ", model.index,
	          " a = ", model.a,
	          " r_in = ", model.r_in,
	          " r_abs = ", model.r_abs,
	          " r_out = ", model.r_out,
	          " limb = ", model.limb,
	          " K_o = ", model.K_o,
	          " column = ", model.column,
	          " abs_logξ = ", model.abs_logξ,
	          " vturb = ", model.vturb,
	          " z = ", model.z,
	          " K_pl = ", model.K_pl
	
	    # INNER DISC MODEL
	
	    # setup inner disc xillver model
	    m1 = XillverD5()
	    if (model.K_i < 0.0)
	        @warn "K_i is negative, setting to 0.0"
	        m1.K = 0.0
	    else
	        m1.K = model.K_i * 1.0e-12
	    end
	    m1.Γ = model.Γ
	    m1.A_Fe = model.A_Fe
	    m1.logXi = model.logξ
	    m1.density = model.density
	    m1.inclination = model.θ   
	    # setup inner disc absorption model
	    m2 = XS_WarmAbsorber()
	    m2.column = model.column
	    m2.rlogxi = model.abs_logξ
	    m2.Cabund =  1.0
	    m2.Nabund =  1.0
	    m2.Oabund =  1.0
	    m2.Fabund =  1.0
	    m2.Neabund = 1.0
	    m2.Naabund = 1.0
	    m2.Mgabund = 1.0
	    m2.Alabund = 1.0
	    m2.Siabund = 1.0
	    m2.Pabund =  1.0
	    m2.Sabund =  1.0
	    m2.Clabund = 1.0
	    m2.Arabund = 1.0
	    m2.Kabund =  1.0
	    m2.Caabund = 1.0
	    m2.Scabund = 1.0
	    m2.Tiabund = 1.0
	    m2.Vabund =  1.0
	    m2.Crabund = 1.0
	    m2.Mnabund = 1.0
	    m2.Feabund = model.A_Fe
	    m2.Coabund = 1.0
	    m2.Niabund = 1.0
	    m2.Cuabund = 1.0
	    m2.Znabund = 1.0
	    m2.write_outfile = 0.0
	    m2.outfile_idx = 0.0
	    m2.vturb = model.vturb
	    m2.redshift = model.z
	    # setup inner disc relconv model
	    m3_r_in = model.r_in * ISCO(model.a)
	    m3_r_out = model.r_abs * ISCO(model.a)
	    m3_r_break = 0.5*(m3_r_in + m3_r_out)
	    m3 = XS_Relconv()
	    m3.index1 = model.index
	    m3.index2 = model.index
	    m3.r_break = m3_r_break
	    m3.a = model.a
	    m3.θ_obs = model.θ
	    m3.inner_r = m3_r_in
	    m3.outer_r = m3_r_out
	    m3.limb = model.limb
	    # evaluate inner disc model which is relconv(xillver)
	    invokemodel!(our_output, our_domain, m1)
	    # absorption (multiplicative; not convolution but might as well not apply unless xillver model is non-zero)
	    xillver_output = copy(our_output)
	    invokemodel!(our_output, our_domain, m2)
	    our_output .= our_output .* xillver_output
	    # convolution only works if there is something to convolve with (it crashes without the following condition)
	    if maximum(our_output) > 0.0
	        # relconv
			# include the following line if commenting out the convolution
			# @info "WARNING: inner disc convolution turned off"
	        invokemodel!(our_output, our_domain, m3)
	    end
	    # save the inner disc output
	    inner_disc_output = copy(our_output)
	
	    # OUTER DISC MODEL
	
	    # setup outer disc xillver model
	    m4 = XillverD5()
	    if (model.K_o < 0.0)
	        @warn "K_o is negative, setting to 0.0"
	        m4.K = 0.0
	    else
	        m4.K = model.K_o * 1.0e-12
	    end
	    m4.Γ = model.Γ
	    m4.A_Fe = model.A_Fe
	    m4.logXi = model.logξ
	    m4.density = model.density
	    m4.inclination = model.θ
	    # setup outer disc relconv model
	    m5_r_in = m3_r_out
	    m5_r_out = model.r_out * ISCO(model.a)
	    m5_break = 0.5*(m5_r_in + m5_r_out)
	    m5 = XS_Relconv()
	    m5.index1 = model.index
	    m5.index2 = model.index
	    m5.r_break = m5_break
	    m5.a = model.a
	    m5.θ_obs = model.θ
	    m5.inner_r = m5_r_in
	    m5.outer_r = m5_r_out
	    m5.limb = model.limb
	    # evaluate outer disc model which is relconv(xillver)
	    invokemodel!(our_output, our_domain, m4)
	    if maximum(our_output) > 0.0
	        invokemodel!(our_output, our_domain, m5)
	    end
	    # save the outer disc output
	    outer_disc_output = copy(our_output)
	
	    # POWER LAW MODEL
	
	    # setup power law model
	    m6 = XS_PowerLaw()
	    if model.K_pl < 0.0
	        @warn "K_pl is negative, setting to 0.0"
	        m6.K = 0.0
	    else
	        m6.K = model.K_pl * 1.0e-7
	    end
	    m6.a = model.Γ
	    # evaluate power law model
	    invokemodel!(our_output, our_domain, m6)
	    # save the power law output
	    pl_output = copy(our_output)
	
	    # ADD UP THE MODELS
	    our_output = copy(inner_disc_output)
	    our_output .= our_output .+ outer_disc_output
	    our_output .= our_output .+ pl_output
	
	    # return the result from the origin domain excluding the extended bins
	    output .= our_output[length(our_low_bins)+1:length(our_low_bins)+length(output)]
	end
	
end

# ╔═╡ 24dbb752-404a-4e00-bb73-1a91f01d66d7
energy = collect(range(1.0, 12.0, 150))

# ╔═╡ af5a94d2-b2d0-48f3-aa5b-c45cef81daef
model = DiscAbsModel()

# ╔═╡ f005205a-577f-4f39-a976-7742fe64c340
md"""
| Parameter        | Value | Units |
|------------------|-------|-------|
| log Column | $(@bind mod_column Slider(-3:0.1:2, default = 1.0,show_value = true)) | |
| Absorption log ξ | $(@bind mod_abs_logξ Slider(-4:0.1:5, default = 3.0, show_value = true)) | |
| Transition radius | $(@bind mod_r_abs Slider(1.25:0.1:3.0, default = 2.5, show_value = true)) | |
| Inner normalisation | $(@bind mod_K_i Slider(0:0.01:2, default = 1.25, show_value = true))| $10^{-12}$ |
| Outer normalisation | $(@bind mod_K_o Slider(0:0.01:2, default = 0.5, show_value = true)) | $10^{-12}$ |
| Power law mormailisation | $(@bind mod_K_pl Slider(0:0.01:1, default = 0.5, show_value = true)) | $10^{-7}$ |
"""

# ╔═╡ e391d101-cc1b-472a-8242-1bc92c23ac39
begin
	model.r_abs = mod_r_abs
	model.column = mod_column
	model.abs_logξ = mod_abs_logξ
	model.K_i = mod_K_i
	model.K_o = mod_K_o
	model.K_pl = mod_K_pl
	model
	full_model = invokemodel(energy, model)
	# plot model multiplied by energy squared
	p = plot(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*full_model./diff(energy), xscale=:log10, yscale=:log10, xlim=(2.0, 12.0), xlabel="Energy (keV)", xticks = ([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]), label="Full model", legend = :bottomleft)
	# individual components
	# inner
	model.K_i = mod_K_i
	model.K_o = 0.0
	model.K_pl = 0.0
	inner_model = invokemodel(energy, model)
	plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*inner_model./diff(energy), label = "Inner disc")
	# outer
	model.K_i = 0.0
	model.K_o = mod_K_o
	model.K_pl = 0.0
	inner_model = invokemodel(energy, model)
	plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*inner_model./diff(energy), label = "Outer disc")
	# powerlaw
	model.K_i = 0.0
	model.K_o = 0.0
	model.K_pl = mod_K_pl
	inner_model = invokemodel(energy, model)
	plot!(energy[1:end-1], energy[1:end-1].*energy[1:end-1].*inner_model./diff(energy), label = "Power law")
	# reset parameters
	model.K_i = mod_K_i
	model.K_o = mod_K_o
	model.K_pl = mod_K_pl
	p
end

# ╔═╡ Cell order:
# ╠═19887774-4c34-11f0-27ad-27b1f250b3d5
# ╠═a549e35f-3151-4043-ab57-7029119a13d5
# ╠═345c379e-5092-46e7-adf7-decbbfcce731
# ╟─a6104fbc-09ea-4ef6-a915-25c10b9c7f48
# ╟─4ba2e218-81a5-4e14-bb63-8977142b5ee4
# ╠═24dbb752-404a-4e00-bb73-1a91f01d66d7
# ╠═af5a94d2-b2d0-48f3-aa5b-c45cef81daef
# ╟─f005205a-577f-4f39-a976-7742fe64c340
# ╟─e391d101-cc1b-472a-8242-1bc92c23ac39
