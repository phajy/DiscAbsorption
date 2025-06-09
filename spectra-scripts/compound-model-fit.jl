# example fit of compound model to a spectrum

using SpectralFitting, XSPECModels, Plots, Relxill

include("compound-model.jl")

function ISCO(a::Float64)
    Z_1 = 1.0+((1-a^2)^(1/3))*((1+a)^(1/3)+(1-a)^(1/3))
    Z_2 = (3*a^2+Z_1^2)^(1/2)
    if a >= 0 
        r = 3+Z_2-((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    else
        r = 3+Z_2+((3-Z_1)*(3+Z_1+2*Z_2))^(1/2)
    end
end

model = XS_PowerLaw() + OurTestModel()

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
