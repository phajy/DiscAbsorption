using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots
TABLEDIR = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/lowflux/xillverD-5.fits"

function ISCO(a::Float64)
    Z_1 = 1.0+((1-a^2)^(1/3))*((1+a)^(1/3)+(1-a)^(1/3))
    Z_2 = (3*a^2+Z_1^2)^(1/2)
    if a >= 0 
        r = 3+Z_2-((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    else
        r = 3+Z_2+((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    end
end

#= function full_Model(;
    q = FitParam(3.0, lower_limit = -10.0, upper_limit = 10.0,frozen = false),
    a = FitParam(0.998, lower_limit = -0.998, upper_limit = 0.998, frozen = false),
    θ = FitParam(30.0, lower_limit = 4.0, upper_limit = 86.0, frozen = false),
    R_in = FitParam(3.0, lower_limit = ISCO(a.value), upper_limit = 400.0, frozen = false),
    column1 = FitParam(0.0; lower_limit = -3.0, upper_limit = 2.0, frozen = false),
    logξ_1 = FitParam(0.0; lower_limit = -4.0, upper_limit = 5.0, frozen = false),
    A_Fe = FitParam(1.0, lower_limit = 0.5, upper_limit = 20.0, frozen = false),
    logξ_2 = FitParam(0.0; lower_limit = -4.0, upper_limit = 5.0, frozen = false),
    column2 = FitParam(0.0; lower_limit = -3.0, upper_limit = 2.0, frozen = false),
    Γ = FitParam(2.0,lower_limit = 1.2, upper_limit = 3.6, frozen = false),
    logξ_disk = FitParam(0.0; lower_limit = 0.0, upper_limit = 4.6989, frozen = false),
    dens = FitParam(17.0,lower_limit = 15.0, upper_limit = 19.0, frozen = false),
    
    V_turb = FitParam(200.0, frozen = true),
    R_out = FitParam(400.0, frozen = true),
    
    H = FitParam(0.168,frozen = true),
    z = FitParam(0.0658, frozen = true),
)
model = XS_NeutralHydrogenAbsorption(
    nH = H
)*(XS_Relconv(
    index1 = q,
    index2 = q,
    r_break = FitParam(((ISCO(a.value)+R_in.value)/2), frozen = true),
    a = a,
    θ_obs = θ,
    inner_r = FitParam(ISCO(a.value), frozen = true),
    outer_r = R_in,
    limb = FitParam(0.0, frozen = true),
)(XS_WarmAbsorber(
    column = column1,
    rlogxi = logξ_1,
    Cabund = FitParam(1.0; frozen = true),
    Nabund = FitParam(1.0; frozen = true),
    Oabund = FitParam(1.0; frozen = true),
    Fabund = FitParam(1.0; frozen = true),
    Neabund = FitParam(1.0; frozen = true),
    Naabund = FitParam(1.0; frozen = true),
    Mgabund = FitParam(1.0; frozen = true),
    Alabund = FitParam(1.0; frozen = true),
    Siabund = FitParam(1.0; frozen = true),
    Pabund = FitParam(1.0; frozen = true),
    Sabund = FitParam(1.0; frozen = true),
    Clabund = FitParam(1.0; frozen = true),
    Arabund = FitParam(1.0; frozen = true),
    Kabund = FitParam(1.0; frozen = true),
    Caabund = FitParam(1.0; frozen = true),
    Scabund = FitParam(1.0; frozen = true),
    Tiabund = FitParam(1.0; frozen = true),
    Vabund = FitParam(1.0; frozen = true),
    Crabund = FitParam(1.0; frozen = true),
    Mnabund = FitParam(1.0; frozen = true),
    Feabund = A_Fe,
    Coabund = FitParam(1.0; frozen = true),
    Niabund = FitParam(1.0; frozen = true),
    Cuabund = FitParam(1.0; frozen = true),
    Znabund = FitParam(1.0; frozen = true),
    write_outfile = FitParam(0.0; frozen = true),
    outfile_idx = FitParam(0.0; frozen = true),
    vturb = V_turb,
    Redshift = z,
)*XS_WarmAbsorber(
    column = column2,
    rlogxi = logξ_2,
    Cabund = FitParam(1.0; frozen = true),
    Nabund = FitParam(1.0; frozen = true),
    Oabund = FitParam(1.0; frozen = true),
    Fabund = FitParam(1.0; frozen = true),
    Neabund = FitParam(1.0; frozen = true),
    Naabund = FitParam(1.0; frozen = true),
    Mgabund = FitParam(1.0; frozen = true),
    Alabund = FitParam(1.0; frozen = true),
    Siabund = FitParam(1.0; frozen = true),
    Pabund = FitParam(1.0; frozen = true),
    Sabund = FitParam(1.0; frozen = true),
    Clabund = FitParam(1.0; frozen = true),
    Arabund = FitParam(1.0; frozen = true),
    Kabund = FitParam(1.0; frozen = true),
    Caabund = FitParam(1.0; frozen = true),
    Scabund = FitParam(1.0; frozen = true),
    Tiabund = FitParam(1.0; frozen = true),
    Vabund = FitParam(1.0; frozen = true),
    Crabund = FitParam(1.0; frozen = true),
    Mnabund = FitParam(1.0; frozen = true),
    Feabund = A_Fe,
    Coabund = FitParam(1.0; frozen = true),
    Niabund = FitParam(1.0; frozen = true),
    Cuabund = FitParam(1.0; frozen = true),
    Znabund = FitParam(1.0; frozen = true),
    write_outfile = FitParam(0.0; frozen = true),
    outfile_idx = FitParam(0.0; frozen = true),
    vturb = V_turb,
    Redshift = z,
)*XillverD5(
    K = FitParam(1.0),
    Γ = Γ,
    A_Fe = A_Fe,
    logXi = logξ_disk,
    density = dens,
    inclination = θ,
))+
XS_Relconv(
    index1 = q,
    index2 = q,
    r_break = FitParam(((R_in.value+R_out.value)/2), frozen = true),
    a = a,
    θ_obs = θ,
    inner_r = R_in,
    outer_r = R_out,
    limb = FitParam(0.0, frozen = true),
)(XillverD5(
    K = FitParam(1.0),
    Γ = Γ,
    A_Fe = A_Fe,
    logXi = logξ_disk,
    density = dens,
    inclination = θ,
))+
XS_PowerLaw(
    K = FitParam(1.0),
    a = Γ,
))
model
end =#

#Defineing the model

model  = XS_NeutralHydrogenAbsorption()*(XS_Relconv()(XS_WarmAbsorber()*XillverD5())+XS_Relconv()(XillverD5())+PowerLaw())

#freeze all paramaters

for p in SpectralFitting.parameter_tuple(model)
    p.frozen = true
end


#= energy = collect(range(1., 10, 1000))
flux  = invokemodel(energy,model)
plot(energy[1:end-1],flux,yscale = :log10, xscale = :log10) =#

#Load the data 
model

DATADIR = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/lowflux"

spectra = joinpath(DATADIR, "joined_spec_grp.pha")
BKG = joinpath(DATADIR, "joined_spec.bak")
RMF = joinpath(DATADIR, "joined_spec.rsp")
ARF = joinpath(DATADIR, "joined_spec.arf")

data = OGIPDataset(spectra,background=BKG,response=RMF,ancillary=ARF)
regroup!(data) ; normalize!(data) ; drop_bad_channels!(data) ; mask_energies!(data, 1.0, 10.0)

#define the fitting problem

prob = FittingProblem(model => data)

#Set paramaters and bindings
model
SpectralFitting.parameter_named_tuple(model).nH_1.value = 0.168

SpectralFitting.parameter_named_tuple(model).nH_1.value = 0.168

SpectralFitting.parameter_named_tuple(model).vturb_1.value = 200.0

SpectralFitting.parameter_named_tuple(model).redshift_1.value = 0.0658

SpectralFitting.parameter_named_tuple(model).K_1.frozen = false
SpectralFitting.parameter_named_tuple(model).K_2.frozen = false
SpectralFitting.parameter_named_tuple(model).K_3.frozen = false


SpectralFitting.parameter_named_tuple(model).r_break_1.frozen = false
SpectralFitting.parameter_named_tuple(model).r_break_1.value = 15

SpectralFitting.parameter_named_tuple(model).r_break_2.frozen = false
SpectralFitting.parameter_named_tuple(model).r_break_2.value = 2

#=,
z = FitParam(0.0658, frozen = true), =#

SpectralFitting.parameter_named_tuple(model).a_1.frozen = false
SpectralFitting.parameter_named_tuple(model).a_1.lower_limit = 1.2
SpectralFitting.parameter_named_tuple(model).a_1.upper_limit = 3.6

SpectralFitting.parameter_named_tuple(model).index1_1.frozen = false
SpectralFitting.parameter_named_tuple(model).index1_1.lower_limit = -10.0
SpectralFitting.parameter_named_tuple(model).index1_1.upper_limit = 10.0

SpectralFitting.parameter_named_tuple(model).a_2.frozen = false
SpectralFitting.parameter_named_tuple(model).a_2.lower_limit = -0.998
SpectralFitting.parameter_named_tuple(model).a_2.upper_limit = 0.998

SpectralFitting.parameter_named_tuple(model).inclination_1.frozen = false
SpectralFitting.parameter_named_tuple(model).inclination_1.lower_limit = 4.0
SpectralFitting.parameter_named_tuple(model).inclination_1.upper_limit = 86.0

SpectralFitting.parameter_named_tuple(model).outer_r_2.frozen = false
SpectralFitting.parameter_named_tuple(model).outer_r_2.value = 3
SpectralFitting.parameter_named_tuple(model).outer_r_2.lower_limit = 2.0
SpectralFitting.parameter_named_tuple(model).outer_r_2.upper_limit = 10.0

SpectralFitting.parameter_named_tuple(model).column_1.frozen = false
SpectralFitting.parameter_named_tuple(model).column_1.lower_limit = -3.0
SpectralFitting.parameter_named_tuple(model).column_1.upper_limit = 2.0

SpectralFitting.parameter_named_tuple(model).rlogxi_1.frozen = false
SpectralFitting.parameter_named_tuple(model).rlogxi_1.lower_limit = -4.0
SpectralFitting.parameter_named_tuple(model).rlogxi_1.upper_limit = 5.0

SpectralFitting.parameter_named_tuple(model).A_Fe_1.frozen = false
SpectralFitting.parameter_named_tuple(model).A_Fe_1.lower_limit = 0.5
SpectralFitting.parameter_named_tuple(model).A_Fe_1.upper_limit = 20.0

SpectralFitting.parameter_named_tuple(model).logXi_1.frozen = false
SpectralFitting.parameter_named_tuple(model).logXi_1.lower_limit = 0.0
SpectralFitting.parameter_named_tuple(model).logXi_1.upper_limit = 4.6989

SpectralFitting.parameter_named_tuple(model).density_1.frozen = false
SpectralFitting.parameter_named_tuple(model).density_1.lower_limit = 17.0
SpectralFitting.parameter_named_tuple(model).density_1.upper_limit = 19.0

bind!(prob, 1 => :a_1, 1 => :Γ_1, 1 => :Γ_2) 
bind!(prob, 1 => :index1_1, 1 => :index1_2, 1 => :index2_1, 1 => :index2_2) 
bind!(prob, 1 => :a_2, 1 => :a_3) 
bind!(prob, 1 => :inclination_1, 1 => :θ_obs_1, 1 => :inclination_2, 1 => :θ_obs_2) 
bind!(prob, 1 => :outer_r_2, 1 => :inner_r_1,1 => :r_break_2) 
bind!(prob, 1 => :r_break_1, 1 => :outer_r_1) 
bind!(prob, 1 => :A_Fe_1, 1 => :Feabund_1, 1 => :A_Fe_2) 
bind!(prob, 1 => :logXi_1, 1 => :logXi_2) 
bind!(prob, 1 => :density_1, 1 => :density_2) 

details(prob)
prob
#Fit the model to the data

result = fit(prob, LevenbergMarquadt())
update_model!(model, result)

#plot the results 
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10)
plotresult(data, [result], xlims=(1.0, 10.0),yscale = :log10, xscale = :log10)

