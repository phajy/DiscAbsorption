using SpectralFitting, XSPECModels, Relxill, Plots

include("gradus-lamp-post.jl")

function FreezeAll(model)
    for p in SpectralFitting.parameter_tuple(model)
        p.frozen = true
    end
end

function calc_residuals(result)
    r = result
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    @. (obj - y) / sqrt(var)
end

function plot_res(dataA,dataB,result)
    domainA = SpectralFitting.plotting_domain(dataA)
    domainB = SpectralFitting.plotting_domain(dataB)
    dataplot = plot(dataA, yscale=:log10, xscale=:log10, label="$(dataA.user_data.observation_id)"*"FPMA", color=:black, msc=:black, alpha=0.2, xlabel=false)
    plot!(dataplot, dataB, label="$(dataA.user_data.observation_id)"*"FPMB", color=:red, msc=:red, alpha=0.2)
    plot!(dataplot, result[1], color=:black)
    plot!(dataplot, result[2], color=:red)
    resplot = hline([0], linestyle = :dash, color=:blue, xlabel="Energy (keV)", ylabel="Residuals", label=false)
    plot!(resplot, domainA, calc_residuals(result[1]),seriestype=:stepmid, color=:black, alpha=0.7,label = "FPMA χ^2=$(round(sum(result[1].stats)))",)
    plot!(resplot, domainB, calc_residuals(result[2]),seriestype=:stepmid, color=:red, alpha=0.7,label = "FPMB χ^2=$(round(sum(result[2].stats)))",)
    plot(dataplot, resplot, layout = (2,1), link=:x, xscale=:log10, xlims=(3,79), xticks=([3,4,5,6,7,8,9,10,20,30,40,50,60,70,80], ["3", "4", "5", "6", "7", "8","9","10","20","30","40","50","60","70","80"]))
end

function ApplyResult(result,modelA,modelB)
        all_paramsA = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)[1:16]
        all_paramsB = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)[17:end]

    for (p, r) in zip(SpectralFitting.parameter_vector(modelA), all_paramsA)
        set_value!(p, r)
    end

    for (p, r) in zip(SpectralFitting.parameter_vector(modelB), all_paramsB)
        set_value!(p, r)
    end
    details(prob)
end

convmodel = LampPost(
    h = FitParam(5.,lower_limit = 5., upper_limit = 100., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(0.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(Inf, lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85, frozen = false),
    a = FitParam(0.7,lower_limit=0.0,upper_limit=0.998, frozen = false)
)
    
specmodel = XillverD5(
    K = FitParam(0.0,frozen = true),
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 2., frozen = false),
    A_Fe = FitParam(1.0,lower_limit = 0.1, upper_limit = 100., frozen = true),
    logXi = FitParam(3.,lower_limit= 3., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false), 
    inclination = FitParam(30.,lower_limit=27.,upper_limit=33., frozen = false)
)
        
PL = PowerLaw(
    a = FitParam(2.0,lower_limit=1.0,upper_limit=3.0, frozen = false)
)

Abs = PhotoelectricAbsorption(
    ηH = FitParam(0.0,lower_limit=0.0,upper_limit=3.0, frozen = true)
)
convolution_model = AsConvolution(convmodel)
modelA = Constant(value = FitParam(1.0, frozen=true))*Abs*(PL+convolution_model(specmodel))
modelB = Constant(value = FitParam(1.0, frozen=false))*Abs*(PL+convolution_model(specmodel))

PATH = "/Users/er19801/DiscAbsorption/data/NuSTAR/"

SPECA = joinpath(PATH, "nu80402315002A01_sr_grp.pha")
SPECB = joinpath(PATH, "nu80402315002B01_sr_grp.pha")

dataA = OGIPDataset(SPECA)
dataB = OGIPDataset(SPECB)

regroup!(dataA) ; normalize!(dataA) ; drop_bad_channels!(dataA) ; mask_energies!(dataA, 3.0, 79.0)
regroup!(dataB) ; normalize!(dataB) ; drop_bad_channels!(dataB) ; mask_energies!(dataB, 3.0, 79.0)

prob = FittingProblem(modelA => dataA, modelB => dataB)
details(prob)
begin
    bind!(prob, (1, :m2, :ηH) => (2, :m2, :ηH))
    bind!(prob, (1, :a1, :K) => (2, :a1, :K))
    bind!(prob, (1, :a1, :a) => (1, :a2, :Γ) => (2, :a1, :a) => (2, :a2, :Γ))
    #bind!(prob, (1, :c1, :K) => (2, :c1, :K))
    bind!(prob, (1, :c1, :h) => (2, :c1, :h))
    bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination) => (2, :c1, :θ) => (2, :a2, :inclination))
    bind!(prob, (1, :c1, :a) => (2, :c1, :a))
    bind!(prob, (1, :a2, :K) => (2, :a2, :K))
    bind!(prob, (1, :a2, :A_Fe) => (2, :a2, :A_Fe))
    bind!(prob, (1, :a2, :logXi) => (2, :a2, :logXi))
    bind!(prob, (1, :a2, :density) => (2, :a2, :density))
    details(prob)
end

result = fit(prob, LevenbergMarquadt(), autodiff = :finite, verbose = true)

ApplyResult(result,modelA,modelB)
plot_res(dataA,dataB,result)


result_PL_fit = deepcopy(result)

begin
modelA.a2.K = 1 
modelA.a2.K.frozen = false 
modelB.a2.K = 1 
modelB.a2.K.frozen = false 
bind!(prob, (1, :a2, :K) => (2, :a2, :K))
details(prob)
end

result = fit(prob, LevenbergMarquadt(), autodiff = :finite, verbose = true)

result_LP_fit = deepcopy(result)

plot_res(dataA,dataB,result)

ApplyResult(result,modelA,modelB)
details(prob)
modelA.m2.ηH = 0
modelA.m2.ηH.frozen = false
modelA.a2.A_Fe = 1
modelA.a2.A_Fe.frozen = true
modelA.c1.model.θ = 27