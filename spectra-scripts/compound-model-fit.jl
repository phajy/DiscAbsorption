# example fit of compound model to a spectrum

using SpectralFitting, XSPECModels, Plots, Relxill

include("compound-model.jl")

model = DiscAbsModel()
model.logξ.frozen = true
model.r_abs.frozen = true
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

#define the fitting problem

prob = FittingProblem(model => data)

#Fit the model to the data
result = fit(prob, LevenbergMarquadt())
update_model!(model, result)

#plot the results 
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10)
# plotresult(data, [result], xlims=(1.0, 10.0),yscale = :log10, xscale = :log10)
plot!(result)
