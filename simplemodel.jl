using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots, LaTeXStrings #sets all used packages

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


    data.spectrum.data .*= 10e6
    data.spectrum.errors .*= 10e6
end

#define composite model
#comp_model  = PhotoelectricAbsorption()*(XS_Relconv()(AutoCache(XS_WarmAbsorber(),abstol=1e-9)*XillverD5())+XS_Relconv()(XillverD5())+PowerLaw())
#inner_disk = XS_Relconv()(AutoCache(XS_WarmAbsorber(),abstol=1e-9)*XillverD5())
inner_disk = XS_Relconv()(XillverD5()+GaussianLine())
outer_disk = XS_Relconv()(XillverD5())
comp_model  = PhotoelectricAbsorption()*(inner_disk+outer_disk+PowerLaw())

mo = XS_Relxill()

#patcher function which sets apropriate radii ranges
function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c2.a = clamp(p.c2.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.outer_r = p.c1.inner_r > p.c1.outer_r ? p.c1.inner_r*1.1 : p.c1.outer_r
    p.c2.inner_r = p.c1.outer_r
    p.c1.r_break = (p.c1.inner_r+p.c1.outer_r)/2
    p.c2.r_break = (p.c2.inner_r+p.c2.outer_r)/2
    
end

#define the fitting problem on the patched model 
patched_comp_model = ParameterPatch(comp_model; patch = patcher!)

# extend the data range for convolution purposes 
begin
    prob = FittingProblem(patched_comp_model => data)
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

begin
    patched_comp_model.m1.ηH = 0.168

    patched_comp_model.c1.index1 = 12
    #patched_comp_model.c1.index2 = 3
    #patched_comp_model.c1.a = 0.998
    patched_comp_model.c1.θ_obs = 60
    patched_comp_model.c1.inner_r = ISCO(0.998)
    patched_comp_model.c1.outer_r = 15
    #patched_comp_model.c1.limb = 0 

    #patched_comp_model.a1.K = 1
    #patched_comp_model.a1.Γ = 2 
    patched_comp_model.a1.A_Fe = 2
    #patched_comp_model.a1.density = 17
    patched_comp_model.a1.inclination = 60
    
    patched_comp_model.a2.K = -1
    patched_comp_model.a2.K.upper_limit = 0
    patched_comp_model.a2.K.lower_limit = -Inf64
    patched_comp_model.a2.μ = 6.8
    patched_comp_model.a2.σ = 1
    
    #patched_comp_model.c2.index1 = 3
    #patched_comp_model.c2.index2 = 3
    #patched_comp_model.c2.a = 0.998
    patched_comp_model.c2.θ_obs = 60
    patched_comp_model.c2.inner_r = 3
    #patched_comp_model.c2.outer_r = 400
    #patched_comp_model.c1.limb = 0

    #patched_comp_model.a3.K = 1
    #patched_comp_model.a3.Γ = 2 
    patched_comp_model.a3.A_Fe = 2
    #patched_comp_model.a3.density = 17
    patched_comp_model.a3.inclination = 60

    #patched_comp_model.a4.K = 1
    #patched_comp_model.a4.a =2

    bind!(prob, (1, :c1, :index1) => (1, :c1, :index2) => (1, :c2, :index1) => (1, :c2, :index2))
    bind!(prob, (1, :a1, :Γ) => (1, :a3, :Γ) => (1, :a4, :a))
    bind!(prob, (1, :c1, :a) => (1, :c2, :a))
    bind!(prob, (1, :c1, :θ_obs,) => (1, :a1, :inclination) => (1, :c2, :θ_obs) => (1, :a3, :inclination))
    bind!(prob, (1, :c1, :outer_r) => (1, :c2, :inner_r))
    bind!(prob, (1, :a1, :A_Fe) =>(1, :a3, :A_Fe))
    bind!(prob, (1, :a1, :density) => (1, :a3, :density))
    bind!(prob, (1, :a1, :logXi) => (1, :a3, :logXi))
end

FreezeAll!(patched_comp_model)

patched_comp_model.a1.K.frozen = false
patched_comp_model.a3.K.frozen = false
patched_comp_model.a4.K.frozen = false


#Fit the model to the data
result = fit(prob, LevenbergMarquadt(), verbose = true)#, max_iter = 10)



#begin
function ApplyResult(model, result)
    all_params = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)
    for (p, r) in zip(SpectralFitting.parameter_vector(model), all_params)
        set_value!(p, r)
    end
    model
    details(prob)
end


ApplyResult(patched_comp_model,result)

details(prob)

begin
    patched_comp_model.c1.a.frozen = false
    patched_comp_model.c1.θ_obs.frozen = false
    patched_comp_model.c1.outer_r.frozen = false    

    patched_comp_model.a1.K.frozen = false
    patched_comp_model.a1.Γ.frozen = false
    patched_comp_model.a1.A_Fe.frozen = false
    patched_comp_model.a1.logXi.frozen = false
    patched_comp_model.a1.density.frozen = false

    patched_comp_model.a2.K.frozen = false
    patched_comp_model.a2.σ.frozen = false
    
    patched_comp_model.a3.K.frozen = false
    patched_comp_model.a4.K.frozen = false

    details(prob)
end


begin 
    result = fit(prob, LevenbergMarquadt(), verbose = true, max_iter = 100)
    ApplyResult(patched_comp_model,result)
#plot the results 
end

begin
i=1
COLORS_point = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_bars = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_model = ["#00a676","#0052cd","#a619c5"]
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10,color=COLORS_point[i],markerstrokecolor=COLORS_bars[i])
plot!(result, xlims=(1.0, 10.0),yscale = :log10, xscale = :log10,color=COLORS_model[i])
end