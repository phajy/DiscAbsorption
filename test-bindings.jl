using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots

temp = []
vals = @. log10(abs(temp))

plot(sign.(temp))

function SpectralFitting._invoke_guard!(output, domain, model::XS_Relconv{<:Number})
    replace!(output, 0 => 1e-8)
    if any(<=(0), output)
        throw("BAD")
    end
    SpectralFitting.invoke!(output, domain, model)
end

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

i = 1

SPEC = joinpath(PATH[i], "joined_spec_grp.pha")
BKGD = joinpath(PATH[i], "joined_spec.bak")
RMF = joinpath(PATH[i], "joined_spec.rsp")
ARF = joinpath(PATH[i], "joined_spec.arf")

data = OGIPDataset(SPEC,background=BKGD,response=RMF,ancillary=ARF)
regroup!(data) ; normalize!(data) ; drop_bad_channels!(data) ; mask_energies!(data, 1.0, 12.0)

#define composite model
comp_model  = XS_NeutralHydrogenAbsorption()*(XS_Relconv()(XS_WarmAbsorber()*XillverD5())+XS_Relconv()(XillverD5())+PowerLaw())

function get_comp(model, args...)
    bk = getfield(model, :model)

    m = bk
    for a in args
        m = getfield(m, a)
    end
    m
end

function plot_model(prob)
    c1 = get_comp(prob.model.m[1],:right,:left,:left)
    c2 = get_comp(prob.model.m[1],:right,:left,:right)

    domain = SpectralFitting.make_model_domain(ContiguouslyBinned(), data)

    o1 = invokemodel(domain, c1)
    o2 = invokemodel(domain, c2)

    p = plot(domain[1:end-1], o1, yscale = :log10, xscale = :log10)
    plot!(p, domain[1:end-1], o2)
    p
end

function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c2.a = clamp(p.c2.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.r_break = (p.c1.inner_r+p.c1.outer_r)/2
    p.c2.r_break = (p.c2.inner_r+p.c2.outer_r)/2
    @show p.c1.r_break, p.c1.inner_r, p.c1.outer_r
    @show p.c2.r_break, p.c2.inner_r, p.c2.outer_r
    @show p.c1.inner_r
end

patched_comp_model = ParameterPatch(comp_model; patch = patcher!)

prob = FittingProblem(patched_comp_model => data)

function testprob(prob)
    conf = FittingConfig(prob)
    ps = SpectralFitting._get_parameters(conf.parameter_cache, 1)
    ps[conf.parameter_bindings[1]]
end

#= for p in SpectralFitting.parameter_vector(patched_comp_model)
    p.frozen = true
end =#
    
details(prob)
#bind thaw and set parameters

#thaw nH and set to galctic comumn 
patched_comp_model.m1.nH = 0.168

#bind index across relconv
patched_comp_model.c1.index1.frozen = false
patched_comp_model.c1.index2.frozen = false
patched_comp_model.c2.index1.frozen = false
patched_comp_model.c2.index2.frozen = false
bind!(prob, (1, :c1, :index1) => (1, :c1, :index2) => (1, :c2, :index1) => (1, :c2, :index2))
testprob(prob)


#bind spin a
patched_comp_model.c1.a.frozen = false
patched_comp_model.c2.a.frozen = false
bind!(prob, (1, :c1, :a) => (1, :c2, :a))
testprob(prob)
#bind inclination
patched_comp_model.c1.θ_obs.frozen = false
patched_comp_model.c2.θ_obs.frozen = false
bind!(prob, (1, :c1, :θ_obs,) => (1, :a1, :inclination) => (1, :c2, :θ_obs) => (1, :a2, :inclination))
testprob(prob)
# bind radii ??
patched_comp_model.c1.inner_r = ISCO(0.998)
patched_comp_model.c1.outer_r = 8
patched_comp_model.c1.outer_r.frozen = false
patched_comp_model.c1.inner_r.frozen = false
bind!(prob, (1, :c1, :outer_r) => (1, :c2, :inner_r))
testprob(prob)
#bind Fe abundance
bind!(prob, (1, :m2, :Feabund) => (1, :a1, :A_Fe) => (1, :a2, :A_Fe))
testprob(prob)
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
patched_comp_model.m2.vturb = 200
patched_comp_model.m2.vturb.frozen = true
patched_comp_model.m2.redshift = 0.0658
patched_comp_model.m2.redshift.frozen = true
end

#bind powerlaw
bind!(prob, (1, :a1, :Γ) => (1, :a3, :a) => (1, :a2, :Γ))
testprob(prob)
#bind disk density
bind!(prob, (1, :a1, :density) => (1, :a2, :density))

#thaw ionisation
bind!(prob, (1, :a1, :logXi) => (1, :a2, :logXi))


details(prob)


#Fit the model to the data

result = fit(prob, LevenbergMarquadt(), verbose = true)

update_model!(patched_comp_model, result)

#plot the results 
i=1
COLORS_point = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_bars = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_model = ["#00a676","#0052cd","#a619c5"]
plot(data,xlims=(1.0, 12.0),yscale=:log10,xscale=:log10,color=COLORS_point[i],markerstrokecolor=COLORS_bars[i])
plot!(result, xlims=(1.0, 12.0),yscale = :log10, xscale = :log10,color=COLORS_model[i])
 