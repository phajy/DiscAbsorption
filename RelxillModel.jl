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
function SpectralFitting._invoke_guard!(output, domain, model::XS_Relconv{<:Number})
    for i in eachindex(output)
        if output[i] <= 0
            output[i] = eps(Float64)
        end
        #throw("BAD")
    end
    SpectralFitting.invoke!(output, domain, model)
end

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

model = XS_Relxill()+XS_Relconv()(GaussianLine())

function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.r_break = (p.c1.inner_r+p.c1.outer_r)/2    
end

patched_model = ParameterPatch(model; patch = patcher!)

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

patched_model.a1.K.frozen = false
details(prob)

result = fit(prob, LevenbergMarquadt(), verbose = true)#, max_iter = 10)

ApplyResult(patched_model,result)

patched_model.a1.θ_obs = 60
patched_model.a1.inner_r = ISCO(0.998)
patched_model.a1.z = 0.0658

bind!(prob, (1, :a1, :index1) => (1, :c1, :index1))
bind!(prob, (1, :a1, :index2) => (1, :c1, :index2))
bind!(prob, (1, :a1, :r_break) => (1, :c1, :r_break))
bind!(prob, (1, :a1, :a) => (1, :c1, :a))
bind!(prob, (1, :a1, :θ_obs) => (1, :c1, :θ_obs))
bind!(prob, (1, :a1, :inner_r) => (1, :c1, :inner_r))
bind!(prob, (1, :a1, :outer_r) => (1, :c1, :outer_r))

patched_model.a2.K = -1
patched_model.a2.K.upper_limit = 0
patched_model.a2.K.lower_limit = -Inf64
patched_model.a2.μ = 6.8

patched_model.a1.a.frozen = false
patched_model.a1.θ_obs.frozen = false
patched_model.a1.Gamma.frozen = false
patched_model.a1.logxi.frozen = false
patched_model.a1.Afe.frozen = false

patched_model.a2.K.frozen = false
patched_model.a2.σ = 1e-4

details(prob)
#begin

result = fit(prob, LevenbergMarquadt(), verbose = true)#, max_iter = 10)

ApplyResult(patched_model, result)




begin
i=1
COLORS_point = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_bars = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_model = ["#00a676","#0052cd","#a619c5"]
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10,color=COLORS_point[i],markerstrokecolor=COLORS_bars[i])
plot!(result, xlims=(1.0, 10.0),yscale = :log10, xscale = :log10,color=COLORS_model[i])
end 

function calc_residuals(result)
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result[1]
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    @. (obj - y) / sqrt(var)
end

domain = SpectralFitting.plotting_domain(data)

rp = hline([0], linestyle = :dash, legend = false)
plot!(rp,domain, calc_residuals(result), seriestype = :stepmid)

details(prob)


ApplyResult(patched_model, result)
patched_model.a1.Gamma = 1
energy= collect(range(1,12,1000))
fullmodel = invokemodel(energy, patched_model)

plot!(energy[1:end-1],fullmodel,yscale=:log10,xscale=:log10)

patched_model.a1.K = 0

Gaussian =  invokemodel(energy, patched_model)

plot!(energy[1:end-1],Gaussian,yscale=:log10,xscale=:log10)