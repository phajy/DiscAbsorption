using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots

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
DATADIR = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/"
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

function patcher!(p)
    p.c1.a = clamp(p.c1.a, 0, 0.998)
    p.c2.a = clamp(p.c2.a, 0, 0.998)
    p.c1.inner_r = ISCO(p.c1.a)
    p.c1.r_break = (p.c1.inner_r+p.c1.outer_r)/2
    p.c2.r_break = (p.c2.inner_r+p.c2.outer_r)/2
@show p.c1.inner_r
end



patched_comp_model = ParameterPatch(comp_model; patch = patcher!)

prob = FittingProblem(patched_comp_model => data)
for p in SpectralFitting.parameter_vector(patched_comp_model)
    p.frozen = true
end
    
details(prob)
#bind thaw and set parameters

#thaw nH and set to galctic comumn 
patched_comp_model.m1.nH.frozen = false
patched_comp_model.m1.nH = 0.168

#bind index across relconv
patched_comp_model.c1.index1.frozen = false
bind!(prob, (1, :c1, :index1) => (1, :c1, :index2) => (1, :c2, :index1) => (1, :c2, :index2))

#bind spin a
patched_comp_model.c1.a.frozen = false
patched_comp_model.c2.a.frozen = false
patched_comp_model.c1.a = 0.998
bind!(prob, (1, :c1, :a) => (1, :c2, :a))
##
#= #bind inclination
patched_comp_model.c1.θ_obs.frozen = false
bind!(prob, (1, :c1, :θ_obs,) => (1, :a1, :inclination) => (1, :c2, :θ_obs) => (1, :a2, :inclination))
##
# bind radii ??
patched_comp_model.c1.inner_r = ISCO(0.998)
patched_comp_model.c1.outer_r = 8
bind!(prob, (1, :c1, :outer_r) => (1, :c2, :inner_r))
patched_comp_model.c1.outer_r.frozen = false
patched_comp_model.c1.r_break = 4.0
patched_comp_model.c2.r_break = 200.0
##
#bind Fe abundance
patched_comp_model.m2.Feabund.frozen = false
bind!(prob, (1, :m2, :Feabund) => (1, :a1, :A_Fe) => (1, :a2, :A_Fe))
##
#bind powerlaw
patched_comp_model.a1.Γ.frozen = false
bind!(prob, (1, :a3, :a) => (1, :a2, :Γ) => (1, :a1, :Γ))
##
#bind disk density
patched_comp_model.a1.density.frozen = false
patched_comp_model.a2.density.frozen = false
bind!(prob, (1, :a1, :density) => (1, :a2, :density))
##
#thaw all norms
patched_comp_model.a1.K.frozen = false
patched_comp_model.a2.K.frozen = false
patched_comp_model.a3.K.frozen = false
##
#thaw ionisation
patched_comp_model.m2.column.frozen = false
patched_comp_model.m2.rlogxi.frozen = false
patched_comp_model.a1.logXi.frozen = false
patched_comp_model.a2.logXi.frozen = false
bind!(prob, (1, :a1, :logXi) => (1, :a2, :logXi))
##
#set redshift
patched_comp_model.m2.redshift = 0.0658
##
#set vturb
patched_comp_model.m2.vturb = 200

details(prob) =#
##
#Fit the model to the data

result = fit(prob, LevenbergMarquadt(), verbose = true)
#= update_model!(patched_comp_model, result)

#plot the results 
i=2
COLORS_point = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_bars = ["#3b8a00","#5317d4","#cd1e69"]
COLORS_model = ["#00a676","#0052cd","#a619c5"]
plot(data,xlims=(1.0, 12.0),yscale=:log10,xscale=:log10,color=COLORS_point[i],markerstrokecolor=COLORS_bars[i])
plot!(result, xlims=(1.0, 12.0),yscale = :log10, xscale = :log10,color=COLORS_model[i])
 =#