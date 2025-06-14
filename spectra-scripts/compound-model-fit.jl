# example fit of compound model to a spectrum

using SpectralFitting, XSPECModels, Plots, Relxill, Warmabs

include("compound-model.jl")

model = DiscAbsModel()
model.A_Fe = 3.5
model.logξ.frozen = false
model.θ = 70.0
model.θ.frozen = false
model.index.frozen = false
model.r_abs.frozen = false
model.column = 1.9
model.abs_logξ = 4.0
model

#Load the data 

# DATADIR = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/lowflux"
DATADIR = "/Users/phajy/Documents/GitHub/DiscAbsorption/data/lowflux"

spectra = joinpath(DATADIR, "joined_spec_grp.pha")
BKG = joinpath(DATADIR, "joined_spec.bak")
RMF = joinpath(DATADIR, "joined_spec.rsp")
ARF = joinpath(DATADIR, "joined_spec.arf")

data = OGIPDataset(spectra,background=BKG,response=RMF,ancillary=ARF)
regroup!(data) ; normalize!(data) ; drop_bad_channels!(data) ; mask_energies!(data, 1.0, 10.0)

plot(data, xscale=:log10, yscale=:log10)

#define the fitting problem

prob = FittingProblem(model => data)

#Fit the model to the data
result = fit(prob, LevenbergMarquadt(); autodiff = :finite)
update_model!(model, result)

#plot the results 
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10)
# plotresult(data, [result], xlims=(1.0, 10.0),yscale = :log10, xscale = :log10)
plot!(result)
