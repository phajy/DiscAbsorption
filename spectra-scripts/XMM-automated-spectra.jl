using ArgParse
include("automating-spectra-function.jl") 

function parse_arguments()
    s = ArgParseSettings()

    #Add arguments with their default values
    @add_arg_table s begin
        "--PROCDir"
            help="Path to directory containg events files"
            default=pwd()
        "--gtiin"
            help="Good time interval text file path"
            default="gti.txt"
        "--gtiout"
            help="Name for the processed goof time interval file"
            default="gtiset.fits"
        "--sourcein"
            help="Path to source region file"
            default="EPN_src.reg"
        "--backgroundin"
            help="Path to background region file"
            default="EPN_bkg.reg"
        "--SAS_CCF"
            help="Path to the summary calibration files"
            default="../../ODF/ccf.cif"
        "--SAS_CCFPATH"
            help="Path to the general calibration files"
            default="/opt/local/XMM/ccf/"
        "--OutDirName"
            help="Name of the created output directory"
            default="prod"
        "--OutDir"
            help="Path to the created output directory"
            default="."
        "--mincounts"
            help="minimum counts required for each time bin"
            arg_type=Int
            default=20
    end

    return parse_args(ARGS, s)
end

parsed_args = parse_arguments()

#Run reduction function
run_spectra(
    PROCDir = parsed_args["PROCDir"],
    gtiin = parsed_args["gtiin"],
    gtiout = parsed_args["gtiout"],
    sourcein = parsed_args["sourcein"],
    backgroundin = parsed_args["backgroundin"],
    SAS_CCF = parsed_args["SAS_CCF"],
    SAS_CCFPATH = parsed_args["SAS_CCFPATH"],
    OutDirName = parsed_args["OutDirName"],
    OutDir = parsed_args["OutDir"],
    mincounts = parsed_args["mincounts"])

