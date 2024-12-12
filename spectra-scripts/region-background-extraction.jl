using FITSIO
using Suppressor

cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/FlareFilt") #rm IRASTEST
dirs = readdir()
clean = filter(x -> occursin("clean", x), readdir())

for n in eachindex(clean)
    EPN=clean[n]
    EPNSplit = split(EPN,"_")
    pref = EPNSplit[1]*"_"*EPNSplit[2]*"_"*EPNSplit[3]
    EPNsrcName = pref*"_src.fits"
    EPNbkgName = pref*"_bkg.fits"
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
    srcexpression = "(#XMMEA_EP&&(PATTERN<=4)&&((X,Y) IN circle($srcX, $srcY,600)))"
    bkgexpression = "(#XMMEA_EP&&(PATTERN<=4)&&((X,Y) IN circle($bkgX, $bkgY,1200)))"
    rm(pref*"_CoordTemp_src.txt")
    rm(pref*"_CoordTemp_bkg.txt")
    run(`evselect table=$EPN filteredset=$EPNsrcName filtertype='expression' expression=$srcexpression`)
    run(`evselect table=$EPN filteredset=$EPNbkgName filtertype='expression' expression=$bkgexpression`)

end

#evselect table=EPN_0780560101_U002_clean.fits filteredset='EPN_S002_src.fits' filtertype='expression' expression='((RAWX,RAWY) IN circle(37.5909359813425, 190.094882273403,7.5)&&CCDNR==4)'
#evselect table=EPN_0780560101_U002_clean.fits filteredset='EPN_S002_bkg.fits' filtertype='expression' expression='((RAWX,RAWY) IN circle(37.5909359813425, 140.094882273403,15)&&CCDNR==4)'
