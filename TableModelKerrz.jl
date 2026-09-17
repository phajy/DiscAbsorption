using XSPECModels, Relxill, CFITSIO, Plots, Base.Threads, Statistics

t1 = time()
include("KerrRingLinetest.jl")
Threads.nthreads() = 8
min_grp_size = 100

NumbVals = [3, 3, 3, 3, 3, 3, 3, 3]
logged = [0, 1, 0, 0, 0, 1, 0, 0]
deltas = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,]
plineprof=(
    #lineprof parameters
    theta = collect(range(5.0, 85.0, NumbVals[5])),
    #reflection paramters
    A_Fe = collect(logrange(0.1, 10.0, NumbVals[6])),
    logXi = collect(range(1.0, 4.0, NumbVals[7])),
    density = collect(range(15.0, 19.0, NumbVals[8])),
)

pemissivity=(
    #emissivity parameters
    r = collect(range(1.5, 10.0, NumbVals[1])),
    h = collect(logrange(1.5, 50.0, NumbVals[2])),
    a = collect(range(0.1, 0.998, NumbVals[3])),
    PhotonIndex = collect(range(1.0, 3.0, NumbVals[4]))
)
pfree = merge(plineprof,pemissivity)
pfrz=(
R_in = 0.0,
R_out = 400.0
)

SPECTRA_Units = "photons/cm^2/s"
Out_path = "Ring_Kerrz-test-thread.fits"
REDSHIFT = "F"
ESCALE = "F"
AddModel = "T"

ENERGIES_Nbins = 1000
E_Min = 0.1
E_Max = 80.0
    
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

free_param_names = collect(String.(keys(pfree)))
frozen_param_names = collect(String.(keys(pfrz)))

free_param_values = collect(values(pfree))
frozen_param_values = collect(values(pfrz))

PARAMETERS_colsdef = [("NAME", "12A", ""),("METHOD", "J", ""),("INITIAL", "E", ""),("DELTA", "E", ""),("MINIMUM", "E", ""),("BOTTOM", "E", ""),("TOP", "E", ""),("MAXIMUM", "E", ""),("NUMBVALS", "J", ""),("VALUE", string(maximum(NumbVals))*"E", "")]
fits_create_binary_tbl(f, 0, PARAMETERS_colsdef, "PARAMETERS")
fits_write_col(f, 1, 1, 1, free_param_names)
fits_write_col(f, 2, 1, 1, logged) 
fits_write_col(f, 3, 1, 1, median.(free_param_values))
fits_write_col(f, 4, 1, 1, deltas)
fits_write_col(f, 5, 1, 1, [ps[1] for ps in free_param_values])
fits_write_col(f, 6, 1, 1, [ps[1] for ps in free_param_values])
fits_write_col(f, 7, 1, 1, [ps[end] for ps in free_param_values])
fits_write_col(f, 8, 1, 1, [ps[end] for ps in free_param_values])
fits_write_col(f, 9, 1, 1, NumbVals)

for i in eachindex(free_param_values)
    fits_write_col(f, 10, i, 1,free_param_values[i])
end

for i in eachindex(frozen_param_values)
    fits_write_key(f,frozen_param_names[i],frozen_param_names[i],"physical parameter held constant") 
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

iter_params = multiplyparams(free_param_values)

SPECTRA_colsdef = [("PARAMVAL", string(length(free_param_values))*"E", ""),("INTPSPEC", string(ENERGIES_Nbins)*"E", SPECTRA_Units)]
fits_create_binary_tbl(f, length(iter_params), SPECTRA_colsdef, "SPECTRA")
fits_write_col(f, 1, 1, 1, vec(stack(iter_params)))

em_param_values = collect(values(pemissivity))
emissivity_iter = multiplyparams(em_param_values)


# Process in chunks: compute in parallel, then write serially
# Chunk size controls memory usage (number of spectra held in memory at once)
chunk_size_em = min(min_grp_size, cld(length(emissivity_iter), Threads.nthreads()))
em_chunks = Iterators.partition(eachindex(emissivity_iter), chunk_size_em)

#precalculate all the emissivity profiles
for chunk in em_chunks
    chunk_indices = collect(chunk)
    n_in_chunk = length(chunk_indices)
    
    # Pre-allocate buffer for this chunk's results
     # Compute spectra in parallel within the chunk
    Threads.@threads for local_idx in 1:n_in_chunk
        j = chunk_indices[local_idx]
        ps = emissivity_iter[j]
        ensure_emissivity_ring(ps[3],ps[2],ps[1],photon_index=ps[4])
    end
end

chunk_size_free = min(min_grp_size, length(emissivity_iter), cld(length(iter_params), Threads.nthreads()))
chunks_free = Iterators.partition(eachindex(iter_params), chunk_size_free)

for chunk in chunks_free
    chunk_indices = collect(chunk)
    n_in_chunk = length(chunk_indices)
    
    # Pre-allocate buffer for this chunk's results
    chunk_results = Vector{Vector{Float64}}(undef, n_in_chunk)
    
    # Compute spectra in parallel within the chunk
    Threads.@threads for local_idx in 1:n_in_chunk
        j = chunk_indices[local_idx]
        ps = iter_params[j]        
        local_model = FullModelRingKerrz_mod(
                r = FitParam(ps[5]),
                h = FitParam(ps[6]),
                Γ = FitParam(ps[8]),
                R_in = FitParam(frozen_param_values[1]),
                R_out = FitParam(frozen_param_values[2]), 
                θ = FitParam(ps[1]),
                a = FitParam(ps[7]),
                A_Fe = FitParam(ps[2]),
                logXi = FitParam(ps[3]),
                density = FitParam(ps[4]),)
        chunk_results[local_idx] = invokemodel(Energies, local_model).parent[:, 1]
    end
    
    # Write results serially (thread-safe)
    for local_idx in 1:n_in_chunk
         j = chunk_indices[local_idx]
        println("Writing $j")
        fits_write_col(f, 2, j, 1, chunk_results[local_idx])
    end
end

fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
fits_write_key(f,"HDUCLAS2", "MODEL SPECTRA", "")
fits_write_key(f,"HDUVERS", "1.0.0", "format version")
close(f)
#54145.161617 seconds

