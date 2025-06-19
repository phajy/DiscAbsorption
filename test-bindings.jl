using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots #sets all used packages

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

#paths for spectra
DATADIR = "data"
STATE = ["lowflux","midflux","highflux"]
PATH = joinpath.(DATADIR,STATE)

#for i in eachindex(PATH)

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


#define composite model
comp_model  = PhotoelectricAbsorption()*(XS_Relconv()(AutoCache(XS_WarmAbsorber(),abstol=1e-9)*XillverD5())+XS_Relconv()(XillverD5())+PowerLaw())
#comp_model  = PhotoelectricAbsorption()*(XS_Relconv()(XS_WarmAbsorber()*XillverD5())+XS_Relconv()(XillverD5())+PowerLaw())

#patcher function which sets apropriate radii ranges
function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c2.a = clamp(p.c2.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.outer_r = p.c1.inner_r > p.c1.outer_r ? p.c1.inner_r*1.1 : p.c1.outer_r
    p.c2.inner_r = p.c1.outer_r
    p.c1.r_break = (p.c1.inner_r+p.c1.outer_r)/2
    p.c2.r_break = (p.c2.inner_r+p.c2.outer_r)/2
    
    @show p.a3.K
    @show p.m2.column, p.m2.rlogxi
end

#define the fitting problem on the patched model 
patched_comp_model = ParameterPatch(comp_model; patch = patcher!)
prob = FittingProblem(patched_comp_model => data)
details(prob)

#bind, thaw and set parameters function for later calling
function setup_prob_and_model(prob, patched_comp_model)
    #freeze nH and set to galctic comumn 
    patched_comp_model.m1.ηH = 0.168
    patched_comp_model.m1.ηH.frozen = true

    #bind index across relconv
    patched_comp_model.c1.index1.frozen = false
    patched_comp_model.c1.index2.frozen = false
    patched_comp_model.c2.index1.frozen = false
    patched_comp_model.c2.index2.frozen = false
    bind!(prob, )
    bind!(prob, (1, :c1, :index1) => (1, :c1, :index2) => (1, :c2, :index1) => (1, :c2, :index2) => (1, :a1, :Γ) => (1, :a3, :a) => (1, :a2, :Γ))


    #bind spin a
    patched_comp_model.c1.a.frozen = false
    patched_comp_model.c2.a.frozen = false
    bind!(prob, (1, :c1, :a) => (1, :c2, :a))

    #bind inclination
    patched_comp_model.c1.θ_obs.frozen = false
    patched_comp_model.c2.θ_obs.frozen = false
    bind!(prob, (1, :c1, :θ_obs,) => (1, :a1, :inclination) => (1, :c2, :θ_obs) => (1, :a2, :inclination))
    patched_comp_model.c1.θ_obs = 60.0

    # bind radii ??
    patched_comp_model.c1.inner_r = ISCO(0.998)
    patched_comp_model.c1.outer_r = 3
    patched_comp_model.c1.outer_r.frozen = false
    patched_comp_model.c1.inner_r.frozen = false
    bind!(prob, (1, :c1, :outer_r) => (1, :c2, :inner_r))

    #bind Fe abundance
    bind!(prob, (1, :m2, :Feabund) => (1, :a1, :A_Fe) => (1, :a2, :A_Fe))

    #freeze all abundaces
    begin
    patched_comp_model.m2.Cabund.frozen = true
    patched_comp_model.m2.Nabund.frozen = true
    patched_comp_model.m2.Oabund.frozen = true
    patched_comp_model.m2.Fabund.frozen = true
    patched_comp_model.m2.Neabund.frozen = true
    patched_comp_model.m2.Naabund.frozen = true
    patched_comp_model.m2.Mgabund.frozen = true
    patched_comp_model.m2.Alabund.frozen = true
    patched_comp_model.m2.Siabund.frozen = true
    patched_comp_model.m2.Pabund.frozen = true
    patched_comp_model.m2.Sabund.frozen = true
    patched_comp_model.m2.Clabund.frozen = true
    patched_comp_model.m2.Arabund.frozen = true
    patched_comp_model.m2.Kabund.frozen = true
    patched_comp_model.m2.Caabund.frozen = true
    patched_comp_model.m2.Scabund.frozen = true
    patched_comp_model.m2.Tiabund.frozen = true
    patched_comp_model.m2.Vabund.frozen = true
    patched_comp_model.m2.Crabund.frozen = true
    patched_comp_model.m2.Mnabund.frozen = true
    patched_comp_model.m2.Coabund.frozen = true
    patched_comp_model.m2.Niabund.frozen = true
    patched_comp_model.m2.Cuabund.frozen = true
    patched_comp_model.m2.Znabund.frozen = true
    patched_comp_model.m2.write_outfile.frozen = true
    patched_comp_model.m2.outfile_idx.frozen = true
    patched_comp_model.m2.model.vturb = 200
    patched_comp_model.m2.vturb.frozen = true
    patched_comp_model.m2.column.frozen = false
    patched_comp_model.m2.rlogxi.frozen = false
    patched_comp_model.m2.Feabund.frozen = false
    patched_comp_model.m2.model.redshift = 0.0658
    patched_comp_model.m2.redshift.frozen = true
    patched_comp_model.m2.model.column = 1
    patched_comp_model.m2.model.rlogxi = 4
    end

    #bind disk density
    patched_comp_model.a1.density.frozen = false
    bind!(prob, (1, :a1, :density) => (1, :a2, :density))

    #thaw ionisation
    patched_comp_model.a1.logXi.frozen = false
    bind!(prob, (1, :a1, :logXi) => (1, :a2, :logXi))
end


setup_prob_and_model(prob, patched_comp_model)
details(prob)

for (name, p) in zip(SpectralFitting.parameter_names(patched_comp_model), SpectralFitting.parameter_vector(patched_comp_model))
    if name != :K
        p.frozen = true
    else
        p.frozen = false
        p.value = 1.0
    end
end
patched_comp_model.a3.K = 1e5
patched_comp_model.a3.K.upper_limit = 1e10
details(prob)

#Fit the model to the data

result = fit(prob, LevenbergMarquadt(), verbose = true, max_iter = 10)



begin
    all_params = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)
    for (p, r) in zip(SpectralFitting.parameter_vector(patched_comp_model), all_params)
        set_value!(p, r)
    end
    patched_comp_model
    details(prob)
end

setup_prob_and_model(prob, patched_comp_model)
details(prob)

result = fit(prob, LevenbergMarquadt(), verbose = true, max_iter = 10)


begin
    all_params = SpectralFitting.update_free_parameters!(result.config.parameter_cache, result.u)
    for (p, r) in zip(SpectralFitting.parameter_vector(patched_comp_model), all_params)
        set_value!(p, r)
    end
    patched_comp_model
    details(prob)
end

#plot the results 

i=1
COLORS_point = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_bars = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_model = ["#00a676","#0052cd","#a619c5"]
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10,color=COLORS_point[i],markerstrokecolor=COLORS_bars[i])
plot!(result, xlims=(1.0, 10.0),yscale = :log10, xscale = :log10,color=COLORS_model[i])
