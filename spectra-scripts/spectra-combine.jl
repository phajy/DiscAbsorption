
##########  REMOVE THE FILES INDICATED BEFORE RUNNING  ##########
using FITSIO
using Dates
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FluxSplitSpectra")
dirs = readdir()
for dir in dirs
#dir = dirs[1]
cd(dir)
srcfiles = filter(x -> occursin("src_spec.fits", x), readdir()) 
w = open("spec_list.txt", "w")
arflist = []
arfweight = []
for file in srcfiles
    #file = srcfiles[1]
    f = FITS(file)
    header = read_header(f[1])
    exposure = (Dates.seconds((Dates.DateTime(header["DATE-END"]) - Dates.DateTime(header["DATE-OBS"]))))
    backfile = chop(file,tail=13)*"bkg_spec.fits"
    ancrfile = chop(file,tail=13)*"src_arf.fits"
    respfile = chop(file,tail=13)*"src_rmf.fits"
    run(`fthedit $file keyword=BACKFILE operation=add value=$backfile comment='Name of background file'`)
    run(`fthedit $file keyword=ANCRFILE operation=add value=$ancrfile comment='Name of ARF file'`)
    run(`fthedit $file keyword=RESPFILE operation=add value=$respfile comment='Name of responce file'`)
    write(w, file)
    write(w,"\n")
    push!(arflist,ancrfile)
    push!(arfweight,exposure)
end
close(w)
arfweight = string.(arfweight/sum(arfweight))
w = open("arf_list.txt", "w")
for i in eachindex(arflist)
    write(w, arflist[i])
    write(w, " ")
    write(w, arfweight[i])
    write(w,"\n")
end
close(w)
run(`addspec infil="spec_list.txt" outfil="joined_spec" qaddrmf="yes" qsubback="yes"`)
run(`addarf @arf_list.txt out_ARF=joined_spec.arf`)
run(`grppha joined_spec.pha joined_spec_grp.pha comm="group min 20 & exit"`)
cd("..")
end