using SpectralFitting, XSPECModels, Relxill, Warmabs, CFITSIO, Plots, Base.Threads, Suppressor
# Spectra
include("TableModel.jl")
Threads.nthreads() = 5
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
a = eachindex(iter_params)
chunks = Iterators.partition(a, cld(length(a), Threads.nthreads()))

SPECTRA_colsdef = [("PARAMVAL", string(length(free_param_values))*"E", ""),("INTPSPEC", string(ENERGIES_Nbins)*"E", SPECTRA_Units)]

tasks = map((enumerate(chunks))) do (N, chunk)
    Threads.@spawn begin
        println("<",Threads.threadid(),">","{$N}")
        f = fits_clobber_file("$N.fits")
        fits_create_binary_tbl(f, chunks.n, SPECTRA_colsdef, "SPECTRA")
        for (k,j) in enumerate(chunk)
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
            end
            fits_write_col(f, 2, k, 1, invokemodel(Energies,model).parent[:,1]) #CHECK THIS <<<<------------<<<<
        end
        fits_write_key(f,"HDUCLASS", "OGIP", "format conforms to OGIP standard")
        fits_write_key(f,"HDUCLAS1", "XSPEC TABLE MODEL", "")
        fits_write_key(f,"HDUCLAS2", "MODEL SPECTRA", "")
        fits_write_key(f,"HDUVERS", "1.0.0", "format version")
        close(f)
        println("finished chunk","{$i}")
    end
end
fetch.(tasks)