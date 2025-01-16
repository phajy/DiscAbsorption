using SpectralFitting, Plots, XSPECModels, OptimizationOptimJL
#cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/highflux")

DATADIR = ["/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/highflux","/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/midflux","/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/lowflux"]
i=1
spectra = joinpath(DATADIR[i], "joined_spec_grp.pha")
background = joinpath(DATADIR[i], "joined_spec.bak")
RMF = joinpath(DATADIR[i], "joined_spec.rsp")
ARF = joinpath(DATADIR[i], "joined_spec.arf")

data = OGIPDataset(spectra,background=background,response=RMF,ancillary=ARF)
regroup!(data)
normalize!(data)
drop_bad_channels!(data)
mask_energies!(data, 1.0, 10.0)
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10)

##

model = PowerLaw() + BlackBody()
prob = FittingProblem(model => data)
details(prob)

result = fit(prob, NelderMead())
update_model!(model, result)

plotresult(data, [result], xlims=(1.0, 10.0),yscale = :log10, xscale = :log10)