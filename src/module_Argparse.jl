using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--L"
            help = "Number of sites in the chain."
            arg_type = Int
            default = 4
        "--model"
            arg_type = String
            default = "EHM_ITensor"
            help = "Name of the model or exec_file"
        "--output_name"
            arg_type = String
            default = "results"
            help = "Name of the output."
        "--result_file_name"
            arg_type = String
            help = "Result file name."
            default = "result.jld2"
        "--previous_random_mps"
            arg_type = String
            help = "Optimizing the simulation storing previous ps0 instead of builing it every simulation."
            default = "previous_MPS.jld2"
        "--hamiltonian_mpos"
            arg_type = String
            help = "Optimizing the simulation storing previous ps0 instead of builing it every simulation."
            default = "hamiltonian_MPOS.jld2"
        "--site_inds_path"
            arg_type = String
            help = "Optimizing the simulation storing siteinds in first simulation and loading it in next."
            default = "siteinds.jld2"
        "--code"
            help = "Name of the executable. i.e. Phase_Diagram"
        "--results"
            arg_type = String
            help = "Path for the results folder"
            default = "../results/"
        "--preprocessing_path"
            arg_type = String
            help = "Path for the results folder"
            default = "../preprocessing/"
        "--J"
            help = "Hopping term."
            default = 1.0
        "--Uf"
            help = "Final on-site interaction strength."
            arg_type = Float64
            default = 6.0
        "--Vf"
            help = "Final NN interaction strength."
            arg_type = Float64
            default = 6.0
        "--U0"
            help = "Initial on-site interaction strength."
            default = -6.0
            arg_type = Float64
        "--V0"
            help = "Initial NN interaction strength."
            default = -6.0
            arg_type = Float64
        "--tau"
            help = "Total time for a quench protocol."
            required = false
            arg_type = Float64
        "--Npoints"
            help = "Number of points in a mesh grid for a phase diagram with U x V.
                    It makes a list with V_values = range(V0, stop=Vf, length=Npoints)"
            arg_type = Int
            default = 100
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
