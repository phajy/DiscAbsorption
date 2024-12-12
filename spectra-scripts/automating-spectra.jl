#cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/0792180601/PileUp")
cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/0780560101/PROC")
#cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/0780561301/PROC")

cwd = pwd()
gtiin = "gti.txt"
gtiout = "gtiset.fits"
sourcein = "EPN_src_ex.reg"
backgroundin = "EPN_bkg.reg"

indir = readdir()
EPNEvts = filter(x -> occursin("EPN", x) && occursin(".fits", x) &&! occursin("LtCrv", x), indir)
##MOSEvts = filter(x -> occursin("EMOS", x) && occursin(".fits", x) &&! occursin("LtCrv", x), indir)

#filtering flaring events
#filter the field on the good time interval
##run(`gtibuild file=$gtiin table=gtiset.fits`)
println("gtibuild file=$gtiin table=$gtiout")
time_expression = "'GTI($gtiout,TIME)'"
EPNEvts_time_filt = []
for n in eachindex(EPNEvts)
    EPN=EPNEvts[n]
    EPNSplit = split(EPN,"_")
    EPNFiltName=EPNSplit[1]*"_"*EPNSplit[2]*"_time_filt.fits"
    push!(EPNEvts_time_filt,EPNFiltName)
    println("evselect table=$EPN filtertype=expression filteredset=$EPNFiltName expression=$time_expression")
    ##run(`evselect table=$EPN filtertype=expression filteredset=$EPNFiltName expression='GTI(gtiset.fits,TIME)'`)
end
#extracting source and background region filters and spectra
src_reg = readlines(sourcein)
src_expression = "";
for  n in eachindex(src_reg);
    lines = src_reg[n] 
    if lines[1] == '-';
        lines = replace(lines,"-circle"=>"!((X,Y) in CIRCLE")
        global src_expression = src_expression*lines*")&&"
    elseif first(lines) != '-';
        lines = replace(lines,"circle"=>"((X,Y) in CIRCLE")
        global src_expression = src_expression*lines*")&&"
    end
end
src_expression = chop(src_expression, tail=2)


bkg_reg = readlines(backgroundin)
bkg_expression = "";
for  n in eachindex(bkg_reg);
    lines = bkg_reg[n] 
    if lines[1] == '-';
        lines = replace(lines,"-circle"=>"!((X,Y) in CIRCLE")
        global bkg_expression = bkg_expression*lines*")&&"
    elseif first(lines) != '-';
        lines = replace(lines,"circle"=>"((X,Y) in CIRCLE")
        global bkg_expression = bkg_expression*lines*")&&"
    end
end
bkg_expression = chop(bkg_expression, tail=2)

EPNEvts_src_filt = []
EPNEvts_bkg_filt = []
for n in eachindex(EPNEvts_time_filt)
    EPN=EPNEvts_time_filt[n]
    EPNSplit = split(EPN,"_")
    EPNsrcFiltName=EPNSplit[1]*"_"*EPNSplit[2]*"_src_filt.fits"
    EPNbkgFiltName=EPNSplit[1]*"_"*EPNSplit[2]*"_bkg_filt.fits"
    EPNsrcSpecName=EPNSplit[1]*"_"*EPNSplit[2]*"_src_spec.fits"
    EPNbkgSpecName=EPNSplit[1]*"_"*EPNSplit[2]*"_bkg_spec.fits"
    push!(EPNEvts_src_filt,EPNsrcFiltName)
    push!(EPNEvts_bkg_filt,EPNbkgFiltName)
    println("evselect table=$EPN energycolumn='PI' filteredset=$EPNsrcFiltName filtertype='expression' expression='$src_expression' spectrumset=$EPNsrcSpecName spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479")
    println("evselect table=$EPN energycolumn='PI' filteredset=$EPNbkgFiltName filtertype='expression' expression='$bkg_expression' spectrumset=$EPNbkgSpecName spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479")
end

EPNEvts_src_filt
EPNEvts_bkg_filt

#checking for pile PileUp
for n in eachindex(EPNEvts_src_filt)
    EPNsrc=EPNEvts_src_filt[n]
    EPNbkg=EPNEvts_bkg_filt[n]
    EPNSplit = split(EPNsrc,"_")
    EPNPileUpName=EPNSplit[1]*"_"*EPNSplit[2]*"_epat.ps"
    ##run(`epatplot set=$EPNsrc plotfile=$EPNPileUpName useplotfile=yes withbackgroundset=yes backgroundset=$EPNbkg`)
    println("epatplot set=$EPNsrc plotfile=$EPNPileUpName useplotfile=yes withbackgroundset=yes backgroundset=$EPNbkg")
end
#create rmf and arf for both source and background in order to set backscale
##run(`rmfgen rmfset=mos1_src_rmf.fits spectrumset=mos1_src_pi.fits `)
##run(`rmfgen rmfset=mos1_bkg_rmf.fits spectrumset=mos1_bkg_pi.fits `)
##run(`arfgen arfset=mos1_src_arf.fits spectrumset=mos1_src_pi.fits withrmfset=yes rmfset=mos1_src_rmf.fits withbadpixcorr=yes badpixlocation=mos1_filt_time.fits setbackscale=yes`)
##run(`arfgen arfset=mos1_bkg_arf.fits spectrumset=mos1_bkg_pi.fits withrmfset=yes rmfset=mos1_bkg_rmf.fits withbadpixcorr=yes badpixlocation=mos1_filt_time.fits setbackscale=yes`)
