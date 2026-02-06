using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--L"
            help = "Number of sites in the chain."
            arg_type = Int
            default = 5
        "--model"
            arg_type = String
            default = "EHM_ITensor"
            help = "Name of the model or exec_file"
        "--code"
            help = "Name of the executable. i.e. Phase_Diagram"
        "--results"
            arg_type = String
            help = "Path for the results folder"
        "--J"
            help = "Hopping term."
            default = 1.0
        "--Uf"
            help = "Final on-site interaction strength."
            arg_type = Float64
        "--Vf"
            help = "Final NN interaction strength."
            arg_type = Float64
        "--U0"
            help = "Initial on-site interaction strength."
            default = 0.0
            arg_type = Float64
        "--V0"
            help = "Initial NN interaction strength."
            default = 0.0
            arg_type = Float64
        "--tau"
            help = "Total time for a quench protocol."
            required = false
            arg_type = Float64
        "--Npoints"
            help = "Number of points in a mesh grid for a phase diagram with U x V.
                    It makes a list with V_values = range(V0, stop=Vf, length=Npoints)"
            arg_type = Int
            default = 30
        "--nsweeps"
            help = "(DMRG Parameters) Number of sweeps."
            default = 6
            arg_type = Int
        "--m"
            help = "(DMRG Parameters) bond dimension."
            default = 10
            arg_type = Int
        "--pairs"
            help = "Only run pairs of U and V"
            default = 0
            arg_type = Int
        "--region"
            help = "Specific region to run"
    end
    return parse_args(s)
end
