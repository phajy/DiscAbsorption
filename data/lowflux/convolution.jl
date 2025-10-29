# Check of convolution code

using SpectralFitting, XSPECModels, Relxill, Plots

energy = collect(range(1, 12, 1000))

model = XS_Relconv()(XS_Gaussian())
model.c1.a = 0.0
model.a1.σ = 1.0e-3
model_inv = invokemodel(energy, model)
plot(energy[1:end-1], model_inv)

model.a1.σ = 1.0e-4
model_inv = invokemodel(energy, model)
plot!(energy[1:end-1], model_inv)
