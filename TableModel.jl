using SpectralFitting, XSPECModels, Relxill, Warmabs, CFITSIO, Plots, Base.Threads, Suppressor
Threads.nthreads() = 25
include("gradus-lamp-post.jl")
# as an example I'm initiating a model and setting some parameters. The upper and lower limits set will be the high and low values for which spectra will be created 
convmodel = LampPost(    
    h = FitParam(1.5,lower_limit = 1.5, upper_limit = 15., frozen = false),
    E = FitParam(1.0,lower_limit = 1., upper_limit = 10., frozen = true),
    R_in = FitParam(-1.,lower_limit= -Inf,frozen = true),
    R_out = FitParam(400., lower_limit=-Inf, frozen = true), 
    θ = FitParam(30.,lower_limit=7,upper_limit=85, frozen = true),
    a = FitParam(0.998,lower_limit=-0.998,upper_limit=0.998, frozen = true))

specmodel = XillverD5(
    Γ = FitParam(2.3,lower_limit = 1, upper_limit = 3., frozen = true),
    A_Fe = FitParam(1.0,lower_limit = 1., upper_limit = 100., frozen = true),
    logXi = FitParam(3.,lower_limit= 0., upper_limit = 4.,frozen = false),
    density = FitParam(17., lower_limit=15., upper_limit=19., frozen = false), 
    inclination = FitParam(30.,lower_limit=7,upper_limit=85, frozen = true))


convolution_model = AsConvolution(convmodel)
model = convolution_model(specmodel)

full_model_vals , model_names = SpectralFitting._all_parameters_with_names(model)

if typeof(model).name.name == :CompositeModel
    full_model_symbols = [Symbol.(k) for k in split.(model_names,".")]
else
    full_model_symbols = [[Symbol.(k)] for k in stack(split.(model_names,"."))[2,:]]
end

deleteat!(full_model_vals, findall(x->:K in x,full_model_symbols))
deleteat!(model_names, findall(x->:K in x,full_model_symbols))
deleteat!(full_model_symbols, findall(x->:K in x,full_model_symbols))

SymbolToValues = Dict(full_model_symbols.=>full_model_vals)
NamestoValues = Dict(model_names.=>full_model_vals)

free_param_symbols = filter(x -> SpectralFitting.isfree(SymbolToValues[x]), full_model_symbols)
free_param_names = filter(x -> SpectralFitting.isfree(NamestoValues[x]), model_names)
free_param_values = filter(x -> SpectralFitting.isfree(x), full_model_vals)
frozen_param_symbols = filter(x -> !SpectralFitting.isfree(SymbolToValues[x]), full_model_symbols)
frozen_param_names = filter(x -> !SpectralFitting.isfree(NamestoValues[x]), model_names)
frozen_param_values = filter(x -> !SpectralFitting.isfree(x), full_model_vals)

#function MakeTable(model,outName)
    SPECTRA_Units = "photons/cm^2/s"
    Out_path = "LampPostTest4.fits"
    REDSHIFT = "F"
    ESCALE = "F"
    logged = [0, 0, 0, 0, 0]
    NumbVals = [5, 5, 5, 5, 5]
    ENERGIES_Nbins = 100
    E_Min = 0.1
    E_Max = 20.0
    
    Model_Name = splitpath(Out_path)[end]
    file_name = splitpath(Out_path)[end]
    if length(file_name) <5 
        file_name *= ".fits"
        Out_path *= ".fits"
    elseif file_name[end-4:end] !== ".fits"
        file_name = file_name*".fits"
        Out_path *= ".fits"
    elseif file_name[end-4:end] == ".fits"
        Model_Name = Model_Name[1:end-5]
    end

if modelkind(model) == Additive()
    AddModel = "T"
else
    AddModel = "F"
end

    greek_map = Dict(
        'α' => "alpha",  'β' => "beta",  'γ' => "gamma",  'δ' => "delta",
        'ϵ' => "epsilon",  'ζ' => "zeta",  'η' => "eta",  'θ' => "theta",
        'ι' => "iota",  'κ' => "kappa",  'λ' => "lambda",  'μ' => "mu",
        'ν' => "nu",  'ξ' => "xi",  'ο' => "omicron",  'π' => "pi",
        'ρ' => "rho",  'σ' => "sigma",  'τ' => "tau",  'υ' => "upsilon",
        'ϕ' => "phi", 'χ' => "chi", 'ψ' => "psi", 'ω' => "omega",
        'Α' => "Alpha",  'Β' => "Beta",  'Γ' => "Gamma",  'Δ' => "Delta",
        'Ε' => "Epsilon",  'Ζ' => "Zeta",  'Η' => "Eta",  'Θ' => "Theta",
        'Ι' => "Iota",  'Κ' => "Kappa",  'Λ' => "Lambda",  'Μ' => "Mu",
        'Ν' => "Nu",  'Ξ' => "Xi",  'Ο' => "Omega",  'Π' => "Pi",
        'Ρ' => "Rho",  'Σ' => "Sigma",  'Τ' => "Tau",  'Υ' => "Upsilon",
        'Φ' => "Phi", 'Χ' => "Chi", 'Ψ' => "Psi", 'Ω' => "Omega")

# Primary Header 

f = fits_clobber_file(Out_path)
fits_create_empty_img(f)
fits_write_key(f,"MODLNAME", Model_Name, "the name of the model")
fits_write_key(f,"MODLUNIT", SPECTRA_Units, "the units for the model")
fits_write_key(f,"REDSHIFT", REDSHIFT, "whether the model contains redshift as a parameter")
fits_write_key(f,"ESCALE", ESCALE, " whether escale is to be a parameter")
fits_write_key(f,"ADDMODEL", AddModel, "whether the model is additive or not")
fits_write_key(f,"LOELIMIT", 0, "the model value for energies below those tabulated")
fits_write_key(f,"HIELIMIT", 0, "the model value for energies above those tabulated")

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUVERS", "1.1.0", "format version")

# Parameters 

function ReplaceGreek(TxtVec) 
    for i in eachindex(TxtVec)
        for (g, latin) in greek_map
            TxtVec[i] = replace(TxtVec[i], g => latin)
        end
    end
end

ReplaceGreek(free_param_names)
ReplaceGreek(frozen_param_names)
free_param_names = String.(stack(split.(free_param_names,"."))[2,:])
frozen_param_names = String.(stack(split.(frozen_param_names,"."))[2,:])

params = []
for i in eachindex(free_param_values)
    if logged[i] == 1 
        push!(params,collect(logrange(SpectralFitting.get_lowerlimit(free_param_values[i]),SpectralFitting.get_upperlimit(free_param_values[i]),NumbVals[i])))
    elseif logged[i] == 0
        push!(params,collect(range(SpectralFitting.get_lowerlimit(free_param_values[i]),SpectralFitting.get_upperlimit(free_param_values[i]),NumbVals[i])))
    end
end

PARAMETERS_colsdef = [("NAME", "12A", ""),("METHOD", "J", ""),("INITIAL", "E", ""),("DELTA", "E", ""),("MINIMUM", "E", ""),("BOTTOM", "E", ""),("TOP", "E", ""),("MAXIMUM", "E", ""),("NUMBVALS", "J", ""),("VALUE", string(maximum(NumbVals))*"E", "")]
fits_create_binary_tbl(f, 0, PARAMETERS_colsdef, "PARAMETERS")
fits_write_col(f, 1, 1, 1, free_param_names)
fits_write_col(f, 2, 1, 1, logged) 
fits_write_col(f, 3, 1, 1, SpectralFitting.get_value.(free_param_values))
fits_write_col(f, 4, 1, 1, SpectralFitting.get_error.(free_param_values))
fits_write_col(f, 5, 1, 1, SpectralFitting.get_lowerlimit.(free_param_values))
fits_write_col(f, 6, 1, 1, SpectralFitting.get_lowerlimit.(free_param_values))
fits_write_col(f, 7, 1, 1, SpectralFitting.get_upperlimit.(free_param_values))
fits_write_col(f, 8, 1, 1, SpectralFitting.get_upperlimit.(free_param_values))
fits_write_col(f, 9, 1, 1, NumbVals)
for i in eachindex(params)
    fits_write_col(f, 10, i, 1,params[i])
end

for i in eachindex(frozen_param_values)
    fits_write_key(f,frozen_param_names[i],SpectralFitting.get_value(frozen_param_values[i]),"physical parameter held constant") 
end
fits_write_key(f,"NINTPARM", length(free_param_names), "the number of interpolated parameters")
fits_write_key(f,"NADDPARM", 0, "the number of additional parameters")
fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "PARAMETERS", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")

# Energies

ENERGIES_colsdef = [("ENERG_LO", "E", "keV"),("ENERG_HI", "E", "keV")]
fits_create_binary_tbl(f, ENERGIES_Nbins, ENERGIES_colsdef, "ENERGIES")
Energies = collect(logrange(E_Min,E_Max,ENERGIES_Nbins+1))
E_low = Energies[1:end-1]
E_high = Energies[2:end]
fits_write_col(f, 1, 1, 1, E_low)
fits_write_col(f, 2, 1, 1, E_high)

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "ENERGIES", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")

# Spectra
function addparams(A,B)
    out = []
    for a in A
        for b in B
            push!(out,[a; b])
        end
    end
    return out
end

function multiplyparams(V)
    out = V[1]
    for i in 2:length(V)
        out = addparams(out,V[i])
    end
    return out
end

iter_params = multiplyparams(params)

SPECTRA_colsdef = [("PARAMVAL", string(length(free_param_values))*"E", ""),("INTPSPEC", string(ENERGIES_Nbins)*"E", SPECTRA_Units)]
fits_create_binary_tbl(f, prod(length.(params)), SPECTRA_colsdef, "SPECTRA")
fits_write_col(f, 1, 1, 1, vec(stack(iter_params)))

a = eachindex(iter_params)
chunks = Iterators.partition(a, cld(length(a), Threads.nthreads()))

tasks = map(chunks) do chunk
    Threads.@spawn for j in chunk
        #println(j,"/",length(iter_params),"<",Threads.threadid(),">")
        ps = iter_params[j]
        for i in eachindex(ps)
            if length(free_param_symbols[i]) == 1
                setproperty!(model,free_param_symbols[i][1],ps[i])
            else
                n,m = free_param_symbols[i]
                try
                    setproperty!(getproperty(model,n), m, ps[i])
                catch
                    setproperty!(getproperty(getproperty(model,n),:model), m, ps[i])
                end
             end
            fits_write_col(f, 2, j, 1, invokemodel(Energies,model).parent[:,1]) #CHECK THIS <<<<------------<<<<
        end
    end
end

fetch.(tasks)

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "MODEL SPECTRA", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")
close(f)