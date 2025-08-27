using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots, LaTeXStrings #sets all used packages

# ISCO function for manual range setting
function ISCO(a::Float64)
    Z_1 = 1.0 + ((1 - a^2)^(1 / 3)) * ((1 + a)^(1 / 3) + (1 - a)^(1 / 3))
    Z_2 = (3 * a^2 + Z_1^2)^(1 / 2)
    if a >= 0
        r = 3 + Z_2 - ((3 - Z_1) * (3 + Z_1 + 2 * Z_2))^(1 / 2)
    else
        r = 3 + Z_2 + ((3 - Z_1) * (3 + Z_1 + 2 * Z_2))^(1 / 2)
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
    STATE = ["lowflux", "midflux", "highflux"]
    PATH = joinpath.(DATADIR, STATE)

    #just for the lof flux for now
    i = 1
    SPEC = joinpath(PATH[i], "joined_spec_grp.pha")
    BKGD = joinpath(PATH[i], "joined_spec.bak")
    RMF = joinpath(PATH[i], "joined_spec.rsp")
    ARF = joinpath(PATH[i], "joined_spec.arf")

    #create data and increase scale to help with normalisation fitting 
    data = OGIPDataset(SPEC, background=BKGD, response=RMF, ancillary=ARF)
    regroup!(data)
    normalize!(data)
    drop_bad_channels!(data)
    mask_energies!(data, 1.0, 10.0)

    data.spectrum.data .*= 10e8
    data.spectrum.errors .*= 10e8
end
#model = XS_Relconv()(XillverD5()+GaussianLine())

model = XS_Relxill() + Constant() * XS_Relconv()(GaussianLine())

function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.r_break = (p.c1.inner_r + p.c1.outer_r) / 2
end

patched_model = ParameterPatch(model; patch=patcher!)

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

patched_model.m1.value.frozen = true
patched_model.a1.K = 1.0
patched_model.a1.inner_r = ISCO(0.998)
patched_model.a1.θ_obs = 60.0
patched_model.a1.θ_obs.frozen = false
patched_model.a1.Gamma.frozen = false
patched_model.a1.logxi.frozen = false
patched_model.a1.refl_frac.frozen = false
patched_model.a1.K.frozen = false
patched_model.a1.z = 0.0658
patched_model.a2.K = 1.0e7
patched_model.a2.K.frozen = true
patched_model.a2.μ = 6.9
patched_model.a2.μ.lower_limit = 6.7
patched_model.a2.μ.upper_limit = 6.9
patched_model.a2.μ.frozen = true
patched_model.a2.σ = 1.0e-3
patched_model.c1.outer_r = 10.0
patched_model.c1.outer_r.frozen = false
# details(prob)

# result = fit(prob, LevenbergMarquadt(), verbose=true)#, max_iter = 10)

# ApplyResult(patched_model, result)

bind!(prob, (1, :a1, :index1) => (1, :c1, :index1))
bind!(prob, (1, :a1, :index2) => (1, :c1, :index2))
bind!(prob, (1, :a1, :r_break) => (1, :c1, :r_break))
bind!(prob, (1, :a1, :a) => (1, :c1, :a))
bind!(prob, (1, :a1, :θ_obs) => (1, :c1, :θ_obs))
# bind!(prob, (1, :a1, :inner_r) => (1, :c1, :inner_r))
# bind!(prob, (1, :a1, :outer_r) => (1, :c1, :outer_r))

patched_model.m1.value.frozen = false
patched_model.m1.value = -1.0e-5
patched_model.m1.value.upper_limit = 0.0
patched_model.m1.value.lower_limit = -1.0e-3

patched_model.c1.inner_r = 4.0
patched_model.c1.inner_r.frozen = true
patched_model.c1.outer_r = 5.0
patched_model.c1.outer_r.frozen = true

details(prob)
#begin

result = fit(prob, LevenbergMarquadt(), verbose=true)#, max_iter = 10)

ApplyResult(patched_model, result)




begin
    i = 1
    COLORS_point = ["#3b8a00", "#5317d4", "#cd1e69"]
    COLORS_bars = ["#3b8a00", "#5317d4", "#cd1e69"]
    COLORS_model = ["#00a676", "#0052cd", "#a619c5"]
    plot(data, xlims=(1.0, 10.0), yscale=:log10, xscale=:log10, color=COLORS_point[i], markerstrokecolor=COLORS_bars[i])
    plot!(result, xlims=(1.0, 10.0), yscale=:log10, xscale=:log10, color=COLORS_model[i])
end

function calc_residuals(result)
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result[1]
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    @. (obj - y) / sqrt(var)
end

domain = SpectralFitting.plotting_domain(data)

rp = hline([0], linestyle=:dash, legend=false)
plot!(rp, domain, calc_residuals(result), seriestype=:stepmid)

details(prob)


ApplyResult(patched_model, result)
patched_model.a1.Gamma = 1
energy = collect(range(1, 12, 1000))
fullmodel = invokemodel(energy, patched_model)

plot(energy[1:end-1], fullmodel, yscale=:log10, xscale=:log10)

patched_model.a1.K = 0

gamodel = invokemodel(energy, patched_model)

plot!(energy[1:end-1], gamodel, yscale=:log10, xscale=:log10)

# use MCMC to fit two power laws to data for illustrative purposes
using StatsPlots
using Turing

model = PowerLaw() + PowerLaw()
model.a1.a = 3.5
model.a1.K = 100.0
model.a2.a = 1.8
model.a2.K = 300.0
model

@model function mcmc_model(objective, stddev, f)
    K1 ~ truncated(Normal(100.0, 10.0); lower = 0.0)
    a1 ~ Normal(3.5, 0.5)
    K2 ~ truncated(Normal(300.0, 10.0); lower = 0.0)
    a2 ~ Normal(1.8, 0.5)
    pred = f(K1, a1, K2, a2)
    return objective ~ MvNormal(pred, stddev)
end

config = FittingConfig(FittingProblem(model => data))

mm = mcmc_model(
    get_objective_single(config),
    sqrt.(get_objective_variance_single(config)),
    get_invoke_wrapper_single(config),
)

chain = sample(mm, NUTS(), 5_000)

plot(chain)

import PairPlots, Makie, CairoMakie

table = (; # named tuple syntax
    K1 = vec(chain["K1"]),
    a1 = vec(chain["a1"]),
    K2 = vec(chain["K2"]),
    a2 = vec(chain["a2"])
)

PairPlots.pairplot(table)

# use MCMC to fit data (does not work)
# using StatsPlots
# using Turing

# model = XS_Relxill()
# model.K.upper_limit = 10.0
# model.index1.frozen = false
# model.index1.upper_limit = 10.0
# model.z = 0.0658
# model.logxi.frozen = false
# model.Afe.frozen = false
# model.Afe.upper_limit = 10.0
# model.refl_frac.upper_limit = 10.0
# model

# @model function mcmc_model(objective, stddev, f)
#     K ~ Normal(2.0, 1.0)
#     index1 ~ Normal(3.0, 2.0)
#     θ_obs ~ truncated(Normal(60.0, 15.0); lower = 5, upper = 85)
#     Gamma ~ Normal(2.0, 0.5)
#     logxi ~ Normal(2.0, 1.0)
#     AFe ~ Normal(1.0, 2.0)
#     refl_frac ~ Normal(2.0, 1.0)
#     pred = f(K, index1, θ_obs, Gamma, logxi, AFe, refl_frac)
#     return objective ~ MvNormal(pred, stddev)
# end

# config = FittingConfig(FittingProblem(model => data))

# mm = mcmc_model(
#     get_objective_single(config),
#     sqrt.(get_objective_variance_single(config)),
#     get_invoke_wrapper_single(config),
# )

# chain = sample(mm, NUTS(), 5_000, autodiff=:finite)
