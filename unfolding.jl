using Interpolations, SpectralFitting

begin
    #paths for spectra
    DATADIR = "data"
    STATE = ["lowflux","midflux","highflux"]
    PATH = joinpath.(DATADIR,STATE)

    #just for the lof flux for now
    i = 3
    SPEC = joinpath(PATH[i], "joined_spec_grp.pha")
    BKGD = joinpath(PATH[i], "joined_spec.bak")
    RMF = joinpath(PATH[i], "joined_spec.rsp")
    ARF = joinpath(PATH[i], "joined_spec.arf")

    #create data and increase scale to help with normalisation fitting 
    data = OGIPDataset(SPEC,background=BKGD,response=RMF,ancillary=ARF)
    regroup!(data) ; normalize!(data) ; drop_bad_channels!(data) ; mask_energies!(data, 1.0, 10.0)
end

energy = SpectralFitting.spectrum_energy(data)
enerror = (data.energy_high-data.energy_low)[findall(x -> x == 1, data.data_mask)]./2

spectrum = filter(!iszero,(data.spectrum.data.*data.data_mask))
specerror = filter(!iszero,(data.spectrum.errors.*data.data_mask))

background = data.background.data[findall(x -> x == 1, data.data_mask)]
backerror = data.background.errors[findall(x -> x == 1, data.data_mask)]

EffectiveArea = linear_interpolation((data.ancillary.bins_low.+data.ancillary.bins_low)./2,data.ancillary.effective_area)

y = (spectrum .- background).*energy.^2

yer = (((specerror .* energy.^2).^2) .+ ((backerror .* energy.^2).^2) .+ ((enerror.*2y./energy).^2)).^0.5

scatter(energy,y./EffectiveArea(energy),legend=false,xerror=enerror,yerror = yer./EffectiveArea(energy),markershape = :x, markerstrokewidth = 0.5,markersize= 5.0,xlims=(1.0, 10.0),yscale=:log10,xscale=:log10)

#plot((data.ancillary.bins_low.+data.ancillary.bins_low)./2,data.ancillary.effective_area,xscale = :log10,yscale=:log10)