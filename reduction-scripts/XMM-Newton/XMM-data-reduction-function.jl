     


#This file defines the function from which the data reduction is run. see the file XMM-data-reduction.jl for information on default values and running the function
"""
    function run_reduction(;
        CCF_Path = "/opt/local/XMM/ccf/", path to calibration files
        ODF_Path = pwd(),                 path to ODF directory
        OutDir = "..",                    location he output directory is created
        OutDirName = "PROC",              name of output directory
        EM = true,                        true or false as to if emproc command is run to process th files for the MOS instruments
        EP = true,                        true or false as to if epproc command is run to process th files for the PN instrument
        Filt = true,                      true or false as to if the standard filters are aplied to the output events file
        LtCrv = true,                     true or false as to if a high energy light curve for the whole image is created from the output events file
        timebinsize = 100                 bin size for the output light curve if created
        )
The above are chanable arguments as well as the default values they take. 
heasoft and ESAS must also be initialised before running this program

"""
function run_reduction(;
    CCF_Path = "/opt/local/XMM/ccf/",
    ODF_Path = pwd(),
    OutDir = "..",
    OutDirName = "PROC",
    EM = true,
    EP = true,
    Filt = true,
    LtCrv = true,
    timebinsize = 100
    )
    
    #place operation inside the selected ODF directory
    cd(ODF_Path)
    
    #Run Cif builder with nescasary enviroment variables pointing to cwd and calibration files 
    withenv("SAS_CCFPATH"=>CCF_Path,"SAS_ODF"=>ODF_Path)do
        run(`cifbuild`)
    end

    #Run odf ingest nescasary enviroment variables
    withenv("SAS_CCFPATH"=>CCF_Path,"SAS_ODF"=>ODF_Path,"SAS_CCF"=>"$ODF_Path/ccf.cif")do
        run(`odfingest`)
    end
    
    #Finds the summary file created by odfingest
    SumFile = filter(x -> occursin("SUM.SAS", x), readdir())

    #Creates and moves to a new output directory named as given when function is called
    cd(OutDir)
    mkdir(OutDirName)
    cd(OutDirName)

    #Runs the emproc and eproc data processing if instructed to when function is called
    if EM == true
        withenv("SAS_CCFPATH"=>CCF_Path,"SAS_ODF"=>"$ODF_Path/"*SumFile[1],"SAS_CCF"=>"$ODF_Path/ccf.cif",)do
            run(`emproc`)
        end
    end

    if EP == true
        withenv("SAS_CCFPATH"=>CCF_Path,"SAS_ODF"=>"$ODF_Path/"*SumFile[1],"SAS_CCF"=>"$ODF_Path/ccf.cif",)do
            run(`epproc`)
        end
    end

    #Reads in the event files created and separates them into the posible multiple
    #exposures for each detector 
    EvFiles = filter(x -> occursin("Evts", x), readdir())
    EMOSEvts = filter(x -> occursin("EMOS", x), EvFiles)
    EPNEvts = filter(x -> occursin("EPN", x), EvFiles)

    #printing some information about the observation here may be usefull. Eg observation ID,
        #Loaction, Name?, number of exposures for each instrument and exposure time (future develpment)

    #The following will create and name apopriatly the filterd events files acording to the standard fiters and the light curve if instructed to do so
    if Filt == true &&  LtCrv == true 
        for n in eachindex(EMOSEvts);
            EMOS=EMOSEvts[n]
            EMOSSplit = split(EMOS,"_")
            EMOSFiltName=EMOSSplit[3]*"_"*EMOSSplit[4]*"_StdFilt.fits"
            EMOSLtcrvName=EMOSSplit[3]*"_"*EMOSSplit[4]*"_LtCrv.fits"
            run(`evselect table=$EMOS filtertype=expression filteredset=$EMOSFiltName expression='(PATTERN <= 12) && (PI in [200:12000]) && #XMMEA_EM'`)
            run(`evselect table=$EMOSFiltName withrateset=Y rateset=$EMOSLtcrvName maketimecolumn=Y timebinsize=$timebinsize makeratecolumn=yes expression='#XMMEA_EM && (PI>10000) && (PATTERN==0)'`)
        end
        for n in eachindex(EPNEvts);
            EPN=EPNEvts[n]
            EPNSplit = split(EPN,"_")
            EPNFiltName=EPNSplit[3]*"_"*EPNSplit[4]*"_StdFilt.fits"
            EPNLtcrvName=EPNSplit[3]*"_"*EPNSplit[4]*"_LtCrv.fits"
            run(`evselect table=$EPN filtertype=expression filteredset=$EPNFiltName expression='(PATTERN <= 4)&&(PI in [200:15000])&&#XMMEA_EP&&(FLAG == 0)'`)
            run(`evselect table=$EPNFiltName withrateset=Y rateset=$EPNLtcrvName maketimecolumn=Y timebinsize=$timebinsize makeratecolumn=yes expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'`)
        end
    elseif Filt == true &&  LtCrv == false
        for n in eachindex(EMOSEvts);
            EMOS=EMOSEvts[n]
            EMOSSplit = split(EMOS,"_")
            EMOSFiltName=EMOSSplit[3]*"_"*EMOSSplit[4]*"_StdFilt.fits"
            run(`evselect table=$EMOS filtertype=expression filteredset=$EMOSFiltName expression='(PATTERN <= 12) && (PI in [200:12000]) && #XMMEA_EM'`)
        end
        for n in eachindex(EPNEvts);
            EPN=EPNEvts[n]
            EPNSplit = split(EPN,"_")
            EPNFiltName=EPNSplit[3]*"_"*EPNSplit[4]*"_StdFilt.fits"
            run(`evselect table=$EPN filtertype=expression filteredset=$EPNFiltName expression='(PATTERN <= 4)&&(PI in [200:15000])&&#XMMEA_EP&&(FLAG == 0)'`)
        end
    elseif Filt == false &&  LtCrv == true
        for n in eachindex(EMOSEvts);
            EMOS=EMOSEvts[n]
            EMOSSplit = split(EMOS,"_")
            EMOSLtcrvName=EMOSSplit[3]*"_"*EMOSSplit[4]*"_LtCrv.fits"
            EMOSName=EMOSSplit[3]*"_"*EMOSSplit[4]*".fits"
            cp(EMOS,EMOSName)
            run(`evselect table=$EMOS withrateset=Y rateset=$EMOSLtcrvName maketimecolumn=Y timebinsize=$timebinsize makeratecolumn=yes expression='#XMMEA_EM && (PI>10000) && (PATTERN==0)'`)
        end
        for n in eachindex(EPNEvts);
            EPN=EPNEvts[n]
            EPNSplit = split(EPN,"_")
            EPNLtcrvName=EPNSplit[3]*"_"*EPNSplit[4]*"_LtCrv.fits"
            EPNName=EPNSplit[3]*"_"*EPNSplit[4]*".fits"
            cp(EPN,EPNName)
            run(`evselect table=$EPN withrateset=Y rateset=$EPNLtcrvName maketimecolumn=Y timebinsize=$timebinsize makeratecolumn=yes expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'`)
        end
    elseif Filt == false &&  LtCrv == false
        for n in eachindex(EMOSEvts);
            EMOS=EMOSEvts[n]
            EMOSSplit = split(EMOS,"_")
            EMOSName=EMOSSplit[3]*"_"*EMOSSplit[4]*".fits"
            cp(EMOS,EMOSName)
        end
        for n in eachindex(EPNEvts);
            EPN=EPNEvts[n]
            EPNSplit = split(EPN,"_")
            EPNName=EPNSplit[3]*"_"*EPNSplit[4]*".fits"
            cp(EPN,EPNName)
        end 
    end
end
    #If i want to flesh this out eventualy I could look at suggesting a rate or time interval on which to filter
#this is as far as a general script for data reduction can go. as the next steps require 
    #decisions to be made about what a suitable good time interval and extration region are. 