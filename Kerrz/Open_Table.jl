using Gradus

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

rm("Kerrz/Table", recursive=true)

kerrz = "/data/typhon2/DariusM/kerrz/zig-out/bin/kerrz"

mkdir("Kerrz/Table")

rs = collect(range(1.5,10.0,10))
hs = collect(logrange(1.5,50.0,10))
Γs = collect(range(1.0,3.0,10))
params = [rs, hs, Γs]

listparam = multiplyparams(params)

for p in listparam
    emisivity_out_file = "Kerrz/Table/emsvty_g$(p[3])_h$(p[2])_r$(p[1])"
    run(`$kerrz emissivity --velocity corotate --photon-index $(p[3]) --nthreads $(Threads.nthreads()) --ring-like h:$(p[2]),x:$(p[1]) --output $emisivity_out_file.dat`)
end