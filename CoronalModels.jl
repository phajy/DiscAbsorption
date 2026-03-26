using SpectralFitting, XSPECModels, Relxill, Plots

include("gradus-lamp-post.jl")
function FreezeAll(model)
    for p in SpectralFitting.parameter_tuple(model)
        p.frozen = true
    end
end

function calc_residuals(result)
    # select which result we want (only have one, but for generalisation to multi-model fits)
    r = result[1]
    y = calculate_objective!(r, r.u)
    obj, var = get_objective(r), get_objective_variance(r)
    @. (obj - y) / sqrt(var)
end

function plot_res(data,result)
    domain = SpectralFitting.plotting_domain(data)
    dp = plot(data, yscale=:log10, xscale=:log10)
    plot!(dp, result, label = "χ^2=$(round(sum(result.stats)))")
    rp = hline([0], linestyle = :dash, legend = false, color=:black)
    plot!(rp, domain, calc_residuals(result),seriestype=:stepmid)
    plot(dp, rp, layout = (2,1), link=:x,xscale=:log10,xlims=(3,79),xticks = ([3,4,5,6,7,8,9,10,20,30,40,50,60,70,80], ["3", "4", "5", "6", "7", "8","9","10","20","30","40","50","60","70","80"]))
end

function ApplyResult(result,model)
        all_params = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)
    for (p, r) in zip(SpectralFitting.parameter_vector(model), all_params)
        set_value!(p, r)
    end
    model
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
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 2., frozen = false),
    A_Fe = FitParam(1.0,lower_limit = 0.1, upper_limit = 100., frozen = false),
    logXi = FitParam(3.,lower_limit= 3., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false), 
    inclination = FitParam(30.,lower_limit=27.,upper_limit=33., frozen = false)
)
        
PL = PowerLaw(
    a = FitParam(2.0,lower_limit=1.0,upper_limit=3.0, frozen = false)
)

Abs = PhotoelectricAbsorption(
    ηH = FitParam(0.86,lower_limit=0.0,upper_limit=3.0, frozen = true)
)
convolution_model = AsConvolution(convmodel)
model = Abs*(PL+convolution_model(specmodel))

PATH = "/Users/er19801/DiscAbsorption/data/NuSTAR/"

SPECA = joinpath(PATH, "nu80402315002A01_sr_grp.pha")
#SPECB = joinpath(PATH, "nu80402315002B01_sr_grp.pha")

dataA = OGIPDataset(SPECA)
#dataB = OGIPDataset(SPECB)

regroup!(dataA) ; normalize!(dataA) ; drop_bad_channels!(dataA) ; mask_energies!(dataA, 3.0, 79.0)

prob = FittingProblem(model => dataA)
bind!(prob, (1, :a1, :a) => (1, :a2, :Γ))
bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination))
details(prob)

model.a2.K = 0
model.a2.K.frozen = true
result = fit(prob, LevenbergMarquadt(), autodiff = :finite)

plot_res(dataA,result)

result_PL_fit = deepcopy(result)

ApplyResult(result,model)

model.a2.K = 1
model.a2.K.frozen = false
model.a1.K.frozen = true
model.a1.a.frozen = true
details(prob)
result = fit(prob, LevenbergMarquadt(), autodiff = :finite)

ApplyResult(result,model)


