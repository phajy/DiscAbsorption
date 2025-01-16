using ArgParse
include("XMM-data-reduction-function.jl") 


function parse_arguments()
    s = ArgParseSettings()

    #Add arguments with their default values
    @add_arg_table s begin
        "--CCF_Path"
            help="Path to calibration files"
            default="/opt/local/XMM/ccf/"
        "--ODF_Path"
            help="Path to ODF directory"
            default=pwd()
        "--OutDir"
            help="Location of the output directory"
            default=".."
        "--OutDirName"
            help="Name of output directory"
            default="PROC"
        "--EM"
            help="Run emproc command for MOS instruments"
            arg_type=Bool
            default=true
        "--EP"
            help="Run epproc command for PN instrument"
            arg_type=Bool
            default=true
        "--Filt"
            help="Apply standard filters to events file"
            arg_type=Bool
            default=true
        "--LtCrv"
            help="Create a high energy light curve from the whole events file"
            arg_type=Bool
            default=true
        "--timebinsize"
            help="Bin size for the output light curve"
            arg_type=Int
            default=100
    end

    return parse_args(ARGS, s)
end

parsed_args = parse_arguments()

#Run reduction function
run_reduction(
    CCF_Path = parsed_args["CCF_Path"],
    ODF_Path = parsed_args["ODF_Path"],
    OutDir = parsed_args["OutDir"],
    OutDirName = parsed_args["OutDirName"],
    EM = parsed_args["EM"],
    EP = parsed_args["EP"],
    Filt = parsed_args["Filt"],
    LtCrv = parsed_args["LtCrv"],
    timebinsize = parsed_args["timebinsize"]
)