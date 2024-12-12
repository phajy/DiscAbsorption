using SpectralFitting, Plots, XSPECModels, OptimizationOptimJL

cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/")
files = readdir("FlareFilt")
for n = eachindex(files)
    global EvtsFiles = filter(x -> occursin("clean", x), files)
end

flux = ["lowflux","midflux","highflux"]

n=1
k=7

DATADIR = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectraGroup/"*flux[n]
prefix = chop(EvtsFiles[k], tail=10)

spectra = joinpath(DATADIR, prefix*"src_spec.fits")
background = joinpath(DATADIR, prefix*"bkg_spec.fits")
RMF = joinpath(DATADIR, prefix*"src_rmf.fits")
ARF = joinpath(DATADIR, prefix*"src_arf.fits")
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/JoinedSpectra/")
spec = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/JoinedSpectra/highflux_src_spec.fits"
resp = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra/JoinedSpectra/highflux_src_rsp.fits"
#data = OGIPDataset(spectra, background=background, response=RMF, ancillary=ARF)
data = OGIPDataset(spec,response=resp)
regroup!(data)
normalize!(data)
drop_bad_channels!(data)
mask_energies!(data, 1.0, 10.0)
data
plot(data,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10)

model = PowerLaw() + BlackBody()
prob = FittingProblem(model => data)
details(prob)

result = fit(prob, NelderMead())
update_model!(model, result)

plotresult(data, [result], xlims=(1.0, 10.0),yscale = :log10, xscale = :log10)