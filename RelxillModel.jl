using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots, LaTeXStrings #sets all used packages

# ISCO function for manual range setting
function ISCO(a::Float64)
    Z_1 = 1.0+((1-a^2)^(1/3))*((1+a)^(1/3)+(1-a)^(1/3))
    Z_2 = (3*a^2+Z_1^2)^(1/2)
    if a >= 0 
 r = 3+Z_2-((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    else
 r = 3+Z_2+((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    end
end

# ensures reclonv gets positive non zero fluxes by setting nay such values to system ϵ
#= function SpectralFitting._invoke_guard!(output, domain, model::XS_Relconv{<:Number})
    for i in eachindex(output)
        if output[i] <= 0
            output[i] = eps(Float64)
        end
        #throw("BAD")
    end
    SpectralFitting.invoke!(output, domain, model)
end =#

#load and initialises data
begin
    #paths for spectra
    DATADIR = "data"
    STATE = ["lowflux","midflux","highflux"]
    PATH = joinpath.(DATADIR,STATE)

    #just for the lof flux for now
    i = 1
    SPEC = joinpath(PATH[i], "joined_spec_grp.pha")
    BKGD = joinpath(PATH[i], "joined_spec.bak")
    RMF = joinpath(PATH[i], "joined_spec.rsp")
    ARF = joinpath(PATH[i], "joined_spec.arf")

    #create data and increase scale to help with normalisation fitting 
    data = OGIPDataset(SPEC,background=BKGD,response=RMF,ancillary=ARF)
    regroup!(data) ; normalize!(data) ; drop_bad_channels!(data) ; mask_energies!(data, 1.0, 10.0)

    data.spectrum.data .*= 10e8
    data.spectrum.errors .*= 10e8
end
#model = XS_Relconv()(XillverD5()+GaussianLine())

begin
Base.@kwdef struct Constant{T} <: AbstractSpectralModel{T,Multiplicative}
    K::T = FitParam(1.0, frozen = true)
end

function SpectralFitting.invoke!(output, input, model::Constant)
    output .= model.K
end
end

inner_disk = XS_Relxill()+Constant(K = FitParam(-1.0,frozen = true))*(XS_Relconv()(GaussianLine()))
outer_disk = XS_Relxill()
comp_model  = PhotoelectricAbsorption()*(inner_disk+outer_disk)

function patcher!(p)
    p.a1.a = clamp(p.a1.a, 0, 0.998)
    p.a3.a = clamp(p.a3.a, 0, 0.998)
    p.a1.inner_r = ISCO(p.a1.a)
    p.a1.outer_r = p.a1.inner_r > p.a1.outer_r ? p.a1.inner_r*1.1 : p.a1.outer_r
    p.a3.inner_r = p.a1.outer_r
    p.a1.r_break = (p.a1.inner_r+p.a1.outer_r)/2
    p.a3.r_break = (p.a3.inner_r+p.a3.outer_r)/2
    
end
patched_model = ParameterPatch(comp_model; patch = patcher!)

begin
    prob = FittingProblem(patched_model => data)
    append!(prob.data.extension[1].high, logrange(12.1, 50, 100))
    append!(prob.data.extension[1].low, logrange(0.1, 0.9, 100))
    details(prob)
end 

function FreezeAll!(model)
    for p in SpectralFitting.parameter_vector(model)
        p.frozen = true
    end
    details(prob)
end

function ThawAll!(model)
    for p in SpectralFitting.parameter_vector(model)
        p.frozen = false
    end
    details(prob)
end

function ApplyResult(model, result)
    all_params = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)
    for (p, r) in zip(SpectralFitting.parameter_vector(model), all_params)
        set_value!(p, r)
    end
    model
    details(prob)
end

FreezeAll!(patched_model)

begin
    bind!(prob, (1, :a1, :index1) => (1, :c1, :index1) => (1, :a3, :index1))
    bind!(prob, (1, :a1, :index2) => (1, :c1, :index2) => (1, :a3, :index2))
    bind!(prob, (1, :a1, :r_break) => (1, :c1, :r_break))
    bind!(prob, (1, :a1, :a) => (1, :c1, :a) => (1, :a3, :a))
    bind!(prob, (1, :a1, :θ_obs) => (1, :c1, :θ_obs) => (1, :a3, :θ_obs))
    bind!(prob, (1, :a1, :inner_r) => (1, :c1, :inner_r))
    bind!(prob, (1, :a1, :outer_r) => (1, :c1, :outer_r))
    bind!(prob, (1, :a1, :outer_r) => (1, :a3, :inner_r))
    bind!(prob, (1, :a1, :z) => (1, :a3, :z))
    bind!(prob, (1, :a1, :Gamma) => (1, :a3, :Gamma))
    bind!(prob, (1, :a1, :logxi) => (1, :a3, :logxi))    
    bind!(prob, (1, :a1, :Afe) => (1, :a3, :Afe))
    bind!(prob, (1, :a1, :Ecut) => (1, :a3, :Ecut))
    bind!(prob, (1, :a1, :refl_frac) => (1, :a3, :refl_frac))
    
    patched_model.a1.inner_r = ISCO(0.998)
    patched_model.a1.θ_obs = 70
    patched_model.a1.z = 0.0658
    
    patched_model.a2.μ = 6.8
    patched_model.a2.σ = 0.001

    patched_model.a1.outer_r = 6
end

patched_model.a1.K.frozen = false
patched_model.a3.K.frozen = false

details(prob)
##
result = fit(prob, LevenbergMarquadt(), verbose = true)#, max_iter = 10)

ApplyResult(patched_model,result)

patched_model.a1.a.frozen = false
patched_model.a1.θ_obs.frozen = false
patched_model.a1.Gamma.frozen = false
patched_model.a1.logxi.frozen = false
patched_model.a1.Afe.frozen = false
patched_model.a1.refl_frac.frozen = false


##

details(prob)
result = fit(prob, LevenbergMarquadt(), verbose = true, max_iter = 100)

ApplyResult(patched_model, result)

begin
i=1
COLORS_point = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_bars = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_model = ["#00a676","#0052cd","#a619c5"]
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10,color=COLORS_point[i],markerstrokecolor=COLORS_bars[i])
plot!(result, xlims=(1.0, 10.0),yscale = :log10, xscale = :log10,color=COLORS_model[i])
end 