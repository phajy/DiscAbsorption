using SpectralFitting, XSPECModels, Relxill, Warmabs, Plots, ProgressMeter

println("Please type a model to start")


rows = 1000
columns = 1200

minE = 1.0
maxE = 12.0
energy = collect(logrange(minE, maxE, columns))

model = XS_WarmAbsorber()
GaussianLine
model.column = 2.0
Par = :rlogxi
typeof(Par)
minPar = 1.0
maxPar = 5.0 

parameter = collect(range(minPar,maxPar,rows))
model = XS_WarmAbsorber()
#= begin
    model.column = 0
    model.rlogxi = 0
    model.Cabund = 1
    model.Nabund = 1
    model.Oabund = 1
    model.Fabund = 1
    model.Neabund = 1
    model.Naabund = 1
    model.Mgabund = 1
    model.Alabund = 1
    model.Siabund = 1
    model.Pabund = 1
    model.Sabund = 1
    model.Clabund = 1
    model.Arabund = 1
    model.Kabund = 1
    model.Caabund = 1
    model.Scabund = 1
    model.Tiabund = 1
    model.Vabund = 1
    model.Crabund = 1
    model.Mnabund = 1
    model.Feabund = 1
    model.Coabund = 1
    model.Niabund = 1
    model.Cuabund = 1
    model.Znabund = 1
    model.write_outfile = 0
    model.outfile_idx = 0
    model.vturb = 200
    model.redshift = 0.0
end =#

T = zeros(rows, columns-1)

p = Progress(rows, desc="Working")

for i in eachindex(parameter)
    setproperty!(model, Par, parameter[i])
    #model.rlogxi = parameter[i]
    T[i,:] = invokemodel(energy,model)
    next!(p)
end 
model

field(model, Par) = 1.0

tick_labels = round.(10 .^ range(log10(minE), log10(maxE), length(collect(minE:maxE))),digits=2)
tick_positions = [(log10(E/minE) / log10(maxE/minE))*(columns-2) + 1 for E in tick_labels]

heatmap(T,colorbar_title="Transmitted fraction", xrotation=90,
xticks=(tick_positions,string.(tick_labels)),
yticks=(Integer.(ceil.(collect(range(1,rows,(length(collect(minPar:maxPar))-1))))),string.(collect(minPar:maxPar))),
xlabel = "Energy (keV)", 
ylabel = string(Par))