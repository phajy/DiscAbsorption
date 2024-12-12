
##########  REMOVE THE FILES INDICATED BEFORE RUNNING  ##########

using FITSIO
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra")
mkdir("JoinedSpectra")
cd("JoinedSpectra")
highfiles = filter(x -> occursin("time", x), readdir("../highflux", join=true))
midfiles = filter(x -> occursin("time", x), readdir("../midflux", join=true))
lowfiles = filter(x -> occursin("time", x), readdir("../lowflux", join=true))
flux = ["high","mid","low"]
region = ["src_";"bkg_"]
file = ["spec";;"rmf";;"arf"]

matrix = region .* file

for i in axes(matrix,1), j in axes(matrix,2)
    w = open("highflux_"*matrix[i,j]*".txt", "w")
    for k in eachindex(highfiles)
        highpref = chop(highfiles[k], tail=14)
        write(w, highpref*matrix[i,j]*".fits")
        write(w,"\n")
    end
    close(w)

    w = open("midflux_"*matrix[i,j]*".txt", "w")
    for k in eachindex(midfiles)
        midpref = chop(midfiles[k], tail=14)
        write(w, midpref*matrix[i,j]*".fits")
        write(w,"\n")
    end
    close(w)

    w = open("lowflux_"*matrix[i,j]*".txt", "w")
    for k in eachindex(lowfiles)
        lowpref = chop(lowfiles[k], tail=14)
        write(w, lowpref*matrix[i,j]*".fits")
        write(w,"\n")
    end
    close(w)
end
outlist = readdir()
highbkgarf = outlist[1]
highbkgrmf = outlist[2]
highbkgspec = outlist[3]
highsrcarf = outlist[4]
highsrcrmf = outlist[5]
highsrcspec = outlist[6]
lowbkgarf = outlist[7]
lowbkgrmf = outlist[8]
lowbkgspec = outlist[9]
lowsrcarf = outlist[10]
lowsrcrmf = outlist[11]
lowsrcspec = outlist[12]
midbkgarf = outlist[13]
midbkgrmf = outlist[14]
midbkgspec = outlist[15]
midsrcarf = outlist[16]
midsrcrmf = outlist[17]
midsrcspec = outlist[18]

directory = "/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra"
highflux = readlines(highsrcspec)
#for i in eachindex(highflux)
i = 1 
file = highflux[i]
bkg = directory*chop(readlines(highbkgspec)[i], head=2)
arf = directory*chop(readlines(highsrcarf)[i], head=2)
rmf = directory*chop(readlines(highsrcrmf)[i], head=2)
run(`grppha $file "chkey RESPFILE $rmf" "chkey ANCRFILE $arf chkey BACKFILE=$bkg"`)
#end
# 
#run(`addspec infil=$highsrcspec outfil="highflux_src_spec.fits" qaddrmf="yes" qsubback="yes"`)

#run(`epicspeccombine pha=$highsrcspec bkg=$highbkgspec rmf=$highsrcrmf arf=$highsrcarf filepha=highflux_src_spec.fits filersp=highflux_src_rsp.fits`)
#run(`epicspeccombine pha=$midsrcspec bkg=$midbkgspec rmf=$midsrcrmf arf=$midsrcarf filepha=midflux_src_spec.fits filersp=midflux_src_rsp.fits`)
#run(`epicspeccombine pha=$lowsrcspec bkg=$lowbkgspec rmf=$lowsrcrmf arf=$lowsrcarf filepha=lowflux_src_spec.fits filersp=lowflux_src_rsp.fits`)