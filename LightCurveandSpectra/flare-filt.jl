cd("/data/typhon2/DariusM/XMM_Data/IRAS13224-3809/") #rm IRASTEST
#reads observation directories and finds the EPN StdFilt event files for each exposure including the path from teh IRAS directory
dirs = readdir()
EvtsFiles = String[]
for n = 1:12 #12
    path = dirs[n]*"/PROC"
    a = filter(x -> occursin("StdFilt", x) && occursin("EPN", x), readdir(path, join=true))
    append!(EvtsFiles,a)
end

#create flare filtered directory and enter it 
mkdir("FlareFilt")
cd("FlareFilt")

Threads.@threads for n in eachindex(EvtsFiles);
    EPN=EvtsFiles[n]
    EPNpathSplit = split(EPN,"/")
    EPNSplit = split(EPNpathSplit[3],"_")
    pref = EPNSplit[1]*"_"*EPNpathSplit[1]*"_"*EPNSplit[2]
    EPNRateName = pref*"_rate.fits"
    EPNPath = joinpath("..",EPN)
    EPNgtiName = pref*"_gti.fits"
    EPNcleanNAME = pref*"_clean.fits"
    gtiexpression = "#XMMEA_EP&&gti($EPNgtiName,TIME)&&(PI>150)"
    SAS_CCF = "../"*EPNpathSplit[1]*"/ODF/ccf.cif"
    #running commands to filter the flaring events
    run(`evselect table=$EPNPath withrateset=Y rateset=$EPNRateName maketimecolumn=Y timebinsize=100 makeratecolumn=Y expression='#XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'`)
    run(`tabgtigen table=$EPNRateName expression='RATE<=0.3' gtiset=$EPNgtiName`)
    run(`evselect table=$EPNPath withfilteredset=Y filteredset=$EPNcleanNAME destruct=Y keepfilteroutput=T expression=$gtiexpression`)

end