using SpectralFitting, XSPECModels, Relxill, Plots
include("gradus-lamp-post.jl")

function FreezeAll(model)
    for p in SpectralFitting.parameter_vector(model)
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
    plot!(resplot, domainA, calc_residuals(result[1]),seriestype=:stepmid, color=:black, alpha=0.7,label = "FPMA",)
    plot!(resplot, domainB, calc_residuals(result[2]),seriestype=:stepmid, color=:red, alpha=0.7,label = "FPMB",)
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

function notice!(data,lowE,highE)
    if lowE < 1.6 
        @warn "Minumum energy provided is lower than domain. Setting minumum to 1.6keV"
        lowE = 1.6
    end
    if highE > 165.4
        @warn "Maximum energy provided is higher than domain. Setting maximum to 165.4keV"
        highE = 165.4
    end
    mask = getproperty(data, :data_mask)
    lowPI = Int((lowE-1.6)/0.04)+1
    highPI = Int((highE-1.6)/0.04)+1
    for i in eachindex(mask)
        if lowPI <= i <= highPI
            mask[i] = 1
        end
    end
    setproperty!(data, :data_mask, mask)
    println("noticed channels $(lowPI-1)-$(highPI-1)")
end

function ignore!(data,lowE,highE)
    if lowE < 1.6 
        @warn "Minumum energy provided is lower than domain. Setting minumum to 1.6keV"
        lowE = 1.6
    end
    if highE > 165.4
        @warn "Maximum energy provided is higher than domain. Setting maximum to 165.4keV"
        highE = 165.4
    end
    mask = getproperty(data, :data_mask)
    lowPI = Int((lowE-1.6)/0.04)+1
    highPI = Int((highE-1.6)/0.04)+1
    for i in eachindex(mask)
        if lowPI <= i <= highPI
            mask[i] = 0
        end
    end
    setproperty!(data, :data_mask, mask)
    println("ignored channels $(lowPI-1)-$(highPI-1)")
end

convmodel = LampPost(
    h = FitParam(5.,lower_limit = 1.5, upper_limit = 100., frozen = false),
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
        
PL = CutoffPL(
    Γ = FitParam(2.0,lower_limit=1.0,upper_limit=3.0, frozen = false),
    β = FitParam(100.0,lower_limit=10.0,upper_limit=600.0, frozen = false)
)

Abs = PhotoelectricAbsorption(
    ηH = FitParam(0.0,lower_limit=0.0,upper_limit=3.0, frozen = true)
)
convolution_model = AsConvolution(convmodel)
modelA = Constant(value = FitParam(1.0, frozen=true))*Abs*(PL+convolution_model(specmodel))
modelB = Constant(value = FitParam(1.0, frozen=false))*Abs*(PL+convolution_model(specmodel))

#PATH = "/Users/er19801/DiscAbsorption/data/NuSTAR/"
PATH = "/data/typhon2/DariusM/NuStar_Data/MAXI_J1348-630/80402315002/products"

SPECA = joinpath(PATH, "nu80402315002A01_sr_1000.pha")
SPECB = joinpath(PATH, "nu80402315002B01_sr_1000.pha")

dataA = OGIPDataset(SPECA)
dataB = OGIPDataset(SPECB)

regroup!(dataA) ; normalize!(dataA) ; drop_bad_channels!(dataA); mask_energies!(dataA,3.0,79.0) 
regroup!(dataB) ; normalize!(dataB) ; drop_bad_channels!(dataB); mask_energies!(dataB,3.0,79.0) 

begin
prob = FittingProblem(modelA => dataA, modelB => dataB)
    bind!(prob, (1, :m2, :ηH) => (2, :m2, :ηH))
    bind!(prob, (1, :a1, :K) => (2, :a1, :K))
    bind!(prob, (1, :a1, :Γ) => (1, :a2, :Γ) => (2, :a1, :Γ) => (2, :a2, :Γ))
    bind!(prob, (1, :a1, :β) => (2, :a1, :β))
    #bind!(prob, (1, :c1, :K) => (2, :c1, :K))
    bind!(prob, (1, :c1, :h) => (2, :c1, :h))
    bind!(prob, (1, :c1, :θ) => (1, :a2, :inclination) => (2, :c1, :θ) => (2, :a2, :inclination))
    bind!(prob, (1, :c1, :a) => (2, :c1, :a))
    bind!(prob, (1, :a2, :K) => (2, :a2, :K))
    bind!(prob, (1, :a2, :A_Fe) => (2, :a2, :A_Fe))
    bind!(prob, (1, :a2, :logXi) => (2, :a2, :logXi))
    bind!(prob, (1, :a2, :density) => (2, :a2, :density))
    FreezeAll(modelA)
    FreezeAll(modelB)
    modelA.a1.K.frozen = false
    modelA.a1.Γ.frozen = false
    modelA.a1.β.frozen = false
    modelB.m1.value.frozen = false
    details(prob)
end

#= @time begin
result = fit(prob, LevenbergMarquadt(), autodiff = :finite, verbose = true)
end
result_PL_fit = deepcopy(result)
ApplyResult(result,modelA,modelB)
#plot_res(dataA,dataB,result) =#

begin
modelA.c1.h.frozen = false
modelA.c1.θ.frozen = false
modelA.c1.a.frozen = false
modelA.a2.K.frozen = false 
modelA.a2.K = 5 
modelB.a2.K.frozen = false 
modelB.a2.K = 1 
bind!(prob, (1, :a2, :K) => (2, :a2, :K))
modelA.a2.A_Fe.frozen = false
modelA.a2.logXi.frozen = false
modelA.a2.density.frozen = false
details(prob)
end


@time begin
result = fit(prob, LevenbergMarquadt(), autodiff = :finite, verbose = true, maxIter=10000)
end


result_LP_fit = deepcopy(result)
ApplyResult(result,modelA,modelB)
details(prob)
#plot_res(dataA,dataB,result)