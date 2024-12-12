    #on my data takes about 6min to run
    #defining variables 
"""
    function run_spectra(;
    PROCDir = pwd(),
    gtiin = "gti.txt",
    gtiout = "gtiset.fits",
    sourcein = "EPN_src.reg",
    backgroundin = "EPN_bkg.reg",
    SAS_CCF = "../../ODF/ccf.cif",
    SAS_CCFPATH = "/opt/local/XMM/ccf/",
    OutDirName = "prod",
    OutDir = ".",
    mincounts = 20)

The above are chanable arguments as well as the default values they take. 
heasoft and ESAS must also be initialised before running this program

""" 
function run_spectra(;
    PROCDir = pwd(),
    gtiin = "gti.txt",
    gtiout = "gtiset.fits",
    sourcein = "EPN_src.reg",
    backgroundin = "EPN_bkg.reg",
    SAS_CCF = "../../ODF/ccf.cif",
    SAS_CCFPATH = "/opt/local/XMM/ccf/",
    OutDirName = "prod",
    OutDir = ".",
    mincounts = 20)

    #setting up the correct working directory and files
    cd(PROCDir)
    indir = readdir()
    mkdir(joinpath(OutDir,OutDirName))
    cd(joinpath(OutDir,OutDirName))
    println(pwd())
    EPNEvts = filter(x -> occursin("EPN", x) && occursin(".fits", x) &&! occursin("LtCrv", x), indir)
    #mabey introduce in future to work on mos. will be a pain though 
    ##MOSEvts = filter(x -> occursin("EMOS", x) && occursin(".fits", x) &&! occursin("LtCrv", x), indir)

    #filtering flaring events from gti.txt file in directory
    gtiin = joinpath(PROCDir,gtiin)
    run(`gtibuild file=$gtiin table=gtiset.fits`)
    time_expression = "GTI($gtiout,TIME)" #setting the time filtering expression. may not be nescasary. try moving to command later

    # reading the region file and making an expression for filtering. It must be in ciao format and only contain circle regions but should be able to handle anuli
    src_reg = readlines(joinpath(PROCDir,sourcein))
    global src_expression = ""
    for  n in eachindex(src_reg);
        lines = src_reg[n] 
        if lines[1] == '-';
            lines = replace(lines,"-circle"=>"!((X,Y) in CIRCLE")
            global (src_expression = src_expression*lines*")&&")
        elseif first(lines) != '-';
            lines = replace(lines,"circle"=>"((X,Y) in CIRCLE")
            global (src_expression = src_expression*lines*")&&")
        end
    end
    src_expression = chop(src_expression, tail=2)

    # same for reading the background region file
    bkg_reg = readlines(joinpath(PROCDir,backgroundin))
    global bkg_expression = ""
    for  n in eachindex(bkg_reg);
        lines = bkg_reg[n] 
        if lines[1] == '-';
            lines = replace(lines,"-circle"=>"!((X,Y) in CIRCLE")
            global (bkg_expression = bkg_expression*lines*")&&")
        elseif first(lines) != '-';
            lines = replace(lines,"circle"=>"((X,Y) in CIRCLE")
            global (bkg_expression = bkg_expression*lines*")&&")
        end
    end
    bkg_expression = chop(bkg_expression, tail=2)
    println(pwd())

    #creating the gti filtered file, source and background region filters, source and background spectra, the test for pile up and the source and background rmf and arf.
    for n in eachindex(EPNEvts)
        EPN=EPNEvts[n]
        EPNSplit = split(EPN,"_")
        EPNPrefix = EPNSplit[1]*"_"*EPNSplit[2]
        EPNtimeFiltName=EPNPrefix*"_time_filt.fits"
        EPNsrcFiltName=EPNPrefix*"_src_filt.fits"
        EPNbkgFiltName=EPNPrefix*"_bkg_filt.fits"
        EPNsrcSpecName=EPNPrefix*"_src_spec.fits"
        EPNbkgSpecName=EPNPrefix*"_bkg_spec.fits"
        EPNPileUpName=EPNPrefix*"_epat.ps"
        EPNsrcRMFName=EPNPrefix*"_src_rmf.fits"
        EPNbkgRMFName=EPNPrefix*"_bkg_rmf.fits"
        EPNsrcARFName=EPNPrefix*"_src_arf.fits"
        EPNbkgARFName=EPNPrefix*"_bkg_arf.fits"
        EPN = joinpath(PROCDir,EPN)

        withenv("SAS_CCF"=>SAS_CCF,"SAS_CCFPATH"=>SAS_CCFPATH)do
            run(`evselect table=$EPN filtertype=expression filteredset=$EPNtimeFiltName expression=$time_expression`)
            run(`evselect table=$EPNtimeFiltName energycolumn='PI' filteredset=$EPNsrcFiltName filtertype='expression' expression=$src_expression spectrumset=$EPNsrcSpecName spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479`)
            run(`evselect table=$EPNtimeFiltName energycolumn='PI' filteredset=$EPNbkgFiltName filtertype='expression' expression=$bkg_expression spectrumset=$EPNbkgSpecName spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479`)
            run(`epatplot set=$EPNsrcFiltName plotfile=$EPNPileUpName useplotfile=yes withbackgroundset=yes backgroundset=$EPNbkgFiltName`)
            run(`rmfgen rmfset=$EPNsrcRMFName spectrumset=$EPNsrcSpecName`)
            run(`rmfgen rmfset=$EPNbkgRMFName spectrumset=$EPNbkgSpecName `)
            run(`arfgen arfset=$EPNsrcARFName spectrumset=$EPNsrcSpecName withrmfset=yes rmfset=$EPNsrcRMFName withbadpixcorr=yes badpixlocation=$EPNtimeFiltName setbackscale=yes`)
            run(`arfgen arfset=$EPNbkgARFName spectrumset=$EPNbkgSpecName withrmfset=yes rmfset=$EPNbkgRMFName withbadpixcorr=yes badpixlocation=$EPNtimeFiltName setbackscale=yes`)
            run(`specgroup spectrumset=$EPNsrcSpecName mincounts=$mincounts`)
        end
    end
end