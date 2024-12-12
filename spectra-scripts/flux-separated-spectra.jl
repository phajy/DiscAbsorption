

##########  RUN THE GTI FUNCTION FOR HIGH MID AND LOW FLUX BEFORE THIS ONE  ########## 

using FITSIO
using Suppressor
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FlareFilt/") #rm IRASTEST
EvtsFiles = filter(x -> occursin("clean", x), readdir())
cd("..")
mkdir("FluxSplitSpectra")
cd("FluxSplitSpectra")
flux = ["highflux","midflux","lowflux"]
mkdir("highflux")
mkdir("midflux")
mkdir("lowflux")

for n in eachindex(EvtsFiles)
    EPN=EvtsFiles[n]
    EPNSplit = split(EPN,"_")
    pref = EPNSplit[1]*"_"*EPNSplit[2]*"_"*EPNSplit[3]
    SAS_CCF = "../"*EPNSplit[2]*"/ODF/ccf.cif"
    EPNClean = "../FlareFilt/"*EPN

    f = FITS(EPNClean)
    header = read_header(f[2])
    RAWX = header["SRCPOSX"]
    RAWY = header["SRCPOSY"] 
    CCDNR = header["CCDSRC"]
    close(f)
    bkgRAWY = RAWY - 50
    withenv("SAS_CCF"=>SAS_CCF,"SAS_CCFPATH"=>"/opt/local/XMM/ccf/")do
        output_src = @capture_out run(`ecoordconv imageset=$EPNClean coordtype=raw x=$RAWX y=$RAWY ccdno=$CCDNR`)
        output_bkg = @capture_out run(`ecoordconv imageset=$EPNClean coordtype=raw x=$RAWX y=$bkgRAWY ccdno=$CCDNR`)
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
    src_expression = "(#XMMEA_EP&&(PATTERN<=4)&&(FLAG == 0)&&(PI in [200:15000])&&((X,Y) IN circle($srcX, $srcY,600)))"
    bkg_expression = "(#XMMEA_EP&&(PATTERN<=4)&&(FLAG == 0)&&(PI in [200:15000])&&((X,Y) IN circle($bkgX, $bkgY,1200)))"
    rm(pref*"_CoordTemp_src.txt")
    rm(pref*"_CoordTemp_bkg.txt")
    SAS_CCF = "../"*SAS_CCF
    
    for k in eachindex(flux)
        cd(flux[k])
        EPNtimeFiltName=pref*"_time_filt.fits"
        EPNsrcFiltName=pref*"_src_filt.fits"
        EPNbkgFiltName=pref*"_bkg_filt.fits"
        EPNsrcSpecName=pref*"_src_spec.fits"
        EPNbkgSpecName=pref*"_bkg_spec.fits"
        EPNPileUpName=pref*"_epat.ps"
        EPNsrcRMFName=pref*"_src_rmf.fits"
        EPNbkgRMFName=pref*"_bkg_rmf.fits"
        EPNsrcARFName=pref*"_src_arf.fits"
        EPNbkgARFName=pref*"_bkg_arf.fits"
        gti = "../../FlareFilt/"*flux[k]*"gti.fits"
        time_expression = "GTI($gti,TIME)"
        EPNCleanPath = joinpath("..",EPNClean)
        withenv("SAS_CCF"=>SAS_CCF,"SAS_CCFPATH"=>"/opt/local/XMM/ccf/")do
            run(`evselect table=$EPNCleanPath filtertype=expression filteredset=$EPNtimeFiltName expression=$time_expression`)
            run(`evselect table=$EPNtimeFiltName energycolumn='PI' filteredset=$EPNsrcFiltName filtertype='expression' expression=$src_expression spectrumset=$EPNsrcSpecName spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479`)
            run(`evselect table=$EPNtimeFiltName energycolumn='PI' filteredset=$EPNbkgFiltName filtertype='expression' expression=$bkg_expression spectrumset=$EPNbkgSpecName spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479`)
            run(`epatplot set=$EPNsrcFiltName plotfile=$EPNPileUpName useplotfile=yes withbackgroundset=yes backgroundset=$EPNbkgFiltName`)
            run(`rmfgen rmfset=$EPNsrcRMFName spectrumset=$EPNsrcSpecName`)
            run(`rmfgen rmfset=$EPNbkgRMFName spectrumset=$EPNbkgSpecName `)
            run(`arfgen arfset=$EPNsrcARFName spectrumset=$EPNsrcSpecName withrmfset=yes rmfset=$EPNsrcRMFName withbadpixcorr=yes badpixlocation=$EPNtimeFiltName setbackscale=yes`)
            run(`arfgen arfset=$EPNbkgARFName spectrumset=$EPNbkgSpecName withrmfset=yes rmfset=$EPNbkgRMFName withbadpixcorr=yes badpixlocation=$EPNtimeFiltName setbackscale=yes`)
            #run(`specgroup spectrumset=$EPNsrcSpecName mincounts=20`)
        end
        cd("..")
    end
end


