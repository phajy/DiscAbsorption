
##########  REMOVE THE *U014* FILES BEFORE RUNNING THIS ONE  ########## 


using FITSIO
using Suppressor
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FlareFilt") #rm IRASTEST
files = readdir(FlareFilt)
for n = eachindex(files)
    global EvtsFiles = filter(x -> occursin("clean", x), files)
end

EvtsFiles

for n in eachindex(EvtsFiles);

EPN = EvtsFiles[n]
EPNSplit = split(EPN,"_")
pref = EPNSplit[1]*"_"*EPNSplit[2]*"_"*EPNSplit[3]
SAS_CCF = "../"*EPNSplit[2]*"/ODF/ccf.cif"

f = FITS(EPN)
header = read_header(f[2])
RAWX = header["SRCPOSX"]
RAWY = header["SRCPOSY"] 
CCDNR = header["CCDSRC"]
close(f)
bkgRAWY = RAWY - 50
withenv("SAS_CCF"=>SAS_CCF,"SAS_CCFPATH"=>"/opt/local/XMM/ccf/")do
    output_src = @capture_out run(`ecoordconv imageset=$EPN coordtype=raw x=$RAWX y=$RAWY ccdno=$CCDNR`)
    output_bkg = @capture_out run(`ecoordconv imageset=$EPN coordtype=raw x=$RAWX y=$bkgRAWY ccdno=$CCDNR`)
    open(pref*"_CoordTemp_src.txt", "w") do io
        write(io, output_src)
    end
    open(pref*"_CoordTemp_bkg.txt", "w") do io
        write(io, output_bkg)
    end
end 
srccoords = split(readlines(pref*"_CoordTemp_src.txt")[5]," ")
srcX = srccoords[4]
srcY = srccoords[5]
bkgcoords = split(readlines(pref*"_CoordTemp_bkg.txt")[5]," ")
bkgX = bkgcoords[4]
bkgY = bkgcoords[5]

srcexpression = "(#XMMEA_EP&&(PATTERN<=4)&&((X,Y) IN circle($srcX, $srcY,600))&&(PI in [300:10000]))"
bkgexpression = "(#XMMEA_EP&&(PATTERN<=4)&&((X,Y) IN circle($bkgX, $bkgY,1200))&&(PI in [300:10000]))"
rm(pref*"_CoordTemp_src.txt")
rm(pref*"_CoordTemp_bkg.txt")

srclightcurveName = pref*"_raw_src_lightcurve.fits"
bkglightcurveName = pref*"_raw_bkg_lightcurve.fits"
srcbkgsubName = pref*"_src_bkgsub_lightcurve.fits"
withenv("SAS_CCF"=>SAS_CCF,"SAS_CCFPATH"=>"/opt/local/XMM/ccf/")do
    run(`evselect table=$EPN energycolumn=PI expression=$srcexpression withrateset=yes rateset=$srclightcurveName timebinsize=1000 maketimecolumn=yes makeratecolumn=yes`)
    run(`evselect table=$EPN energycolumn=PI expression=$bkgexpression withrateset=yes rateset=$bkglightcurveName timebinsize=1000 maketimecolumn=yes makeratecolumn=yes`)
    run(`epiclccorr srctslist=$srclightcurveName eventlist=$EPN outset=$srcbkgsubName bkgtslist=$bkglightcurveName withbkgset=yes applyabsolutecorrections=yes`)
end

end