using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-L"
            help = "Number of sites in the chain."
            arg_type = Int
        "--model"
            help = "Name of the model or exec_file"
        "--results"
            help = "Path for the results folder"
        "-J"
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
        "--nsweeps"
            help = "a positional argument"
            default = 6
            arg_type = Int
        "--Npoints"
            help = "Number of points in a mesh grid for a phase diagram with U x V. 
                    It makes a list with V_values = range(V0, stop=Vf, length=Npoints)"
            arg_type = Int
    end
    return parse_args(s)
end
