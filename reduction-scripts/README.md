# reduction-scripts
## XXM-data-reduction
A simple script created for automating the data reduction of XMM observations on typhon. The script must be able to access the uncompressed observations files in the Observation Data File (`ODF`). The scrip is completely customisable as to where the output files are placed and where it looks for the `ODF`. By default, it will look in the current working directory for observation files and create an output directory above the current called `PROC`. All of this can be changed by calling the arguments `ODF_Path`, `OutDir` and `OutDirName` respectively (see implementation for full arguments list). `CCF_Path` defines the path to the calibration files required to run the SAS processes. `EM` and `EP` are Booleans which define if the processing is run for the MOS-CCD's and/or the pn-CCD respectively. `Filt` and `LtCrv` direct to create a the standard filtered events file and light curve for each exposure. `timebinsize` defines the size of time bins used for creating the light curve file.   

Please be aware, longer observations may take several tens of minutes (10-30min). Once the script is called, interrupting may require the created flies to be deleted before running again.

## Requirements
In order to be able pass arguments from the command line when running code you need to have ArgParse in your julia environment. This can be added with
```
Pkg.add("ArgParse")
```
either added globally to your julia module or in your environment.
Before the script is called `SAS` and therefore also `HESOFT` must be installed and initialised. 

## Implimentation
The reducition code is run as follows:
```
julia XMM-data-reduction.jl --CCF_Path --ODF_Path --OutDir --OutDirName --EM --EP --Filt --LtCrv --timebinsize
```


| Argument    | Default    |
| --- | --- |
|  CCF_Path   |  "/opt/local/XMM/ccf/"   |
|  ODF_Path   |  pwd()   |
|  OutDirName   |  "PROC"   |
|  OutDir   |  ".."   |
|  EM   |  true   |
|  EP   |  true   |
|  Filt   |  true   |
|  LtCrv   |  true   |
|  timebinsize   |  100  [time in seconds] |
