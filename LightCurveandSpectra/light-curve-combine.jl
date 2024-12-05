using FITSIO
using Plots
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FlareFilt") #rm IRASTEST
files = readdir()
for n = eachindex(files)
    global LtCrvs = filter(x -> occursin("bkgsub", x), files)
end

LtCrvs

#bunching light curves
k = FITS(LtCrvs[1])
rate = read(k[2],"RATE")
error = read(k[2],"ERROR")
time = read(k[2],"TIME")
time = time .- time[1]
time_gap = read(k[2],"TIME")
close(k)
observations = []
l = length(LtCrvs)
for i = 1:(l-1)
    fits1 = FITS(LtCrvs[i])
    fits2 = FITS(LtCrvs[i+1])

    next_rate = read(fits2[2],"RATE")
    next_error = read(fits2[2],"ERROR")
    next_time = read(fits2[2], "TIME")

    time_diff = next_time[1] - time[length(time)] + 100 
    nogap_time = next_time .- time_diff
    append!(observations,nogap_time[1])
    append!(rate,next_rate)
    append!(error,next_error)
    append!(time,nogap_time)
    append!(time_gap,next_time)
    close(fits1)
    close(fits2)
end

f = FITS("LightCurve_nogap.fits","w")
write(f, Dict("RATE"=>rate, "ERROR"=>error, "TIME"=>time))
close(f)

f = FITS("LightCurve_gap.fits","w")
write(f, Dict("RATE"=>rate, "ERROR"=>error, "TIME"=>time_gap))
close(f)

function objective(cutoff, rate, target; f = Base.:>)
    counts = sum(filter(x -> !isnan(x) && f(x, first(cutoff)), rate))
    (counts - target)^2
end

target = sum(filter(!isnan, rate)) / 3

objective(3, rate, target; f = Base.:>)

using Optim 

result = optimize(p -> objective(p, rate, target), [2.0], NelderMead())
high, _ = first(Optim.minimizer(result)), Optim.minimum(result)

result = optimize(p -> objective(p, rate, target; f = Base.:<), [2.0], NelderMead())
low, _ = first(Optim.minimizer(result)), Optim.minimum(result)

mat = [time;;rate;;error]

total = sum(filter(!isnan,mat[:,2]))*100
target = total/3

low_filt = mat[mat[:,2] .< low, :]
low_sum = sum(low_filt[:,2])*100

high_filt = mat[mat[:,2] .> high, :]
high_sum = sum(high_filt[:,2])*100

mid_filt = mat[low .< mat[:,2] .< high, :]
mid_sum = sum(mid_filt[:,2])*100

plot(bar(low_filt[:,1],low_filt[:,2],linecolor=:green,fillcolor = :green, fillalpha=1,label="Low flux",xlabel = "Time (s)", ylabel = "Count rate (s^-1)"))
plot(bar!(mid_filt[:,1],mid_filt[:,2],linecolor=:orange,fillcolor = :orange, fillalpha=1,label="Medium flux"))
plot(bar!(high_filt[:,1],high_filt[:,2],linecolor=:pink,fillcolor = :pink, fillalpha=1,label="High flux"))
plot(hline!([low, high],color = :grey,label="cutoff"))
plot(vline!(observations,color = :grey,label="exposure"))

println(high)
println(low)

run(`tabgtigen table=LightCurve_gap.fits expression='RATE<$low' gtiset=lowfluxgti.fits`)
run(`tabgtigen table=LightCurve_gap.fits expression='RATE>$high' gtiset=highfluxgti.fits`)
run(`tabgtigen table=LightCurve_gap.fits expression='RATE<$high && RATE>$low' gtiset=midfluxgti.fits`)