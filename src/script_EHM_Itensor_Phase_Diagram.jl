#=
    Creates a DataFrame object and stores it in a .csv file. 
=#
using DelimitedFiles
using DataFrames
using CSV
using Printf

include("module_Argparse.jl")
include("module_Organize_Results.jl")
let
    parser = parse_commandline()

    # Model parameters
    L = parser["L"]
    model = parser["model"]
    code = parser["code"]

    results = parser["results"]
    J = parser["J"]

    U0 = parser["U0"]
    Uf = parser["Uf"]
    V0 = parser["V0"]
    Vf = parser["Vf"]

    Npoints = parser["Npoints"]

    U_values = round.(range(U0, stop=Uf, length=Npoints), digits=2)
    V_values = round.(range(V0, stop=Vf, length=Npoints), digits=2)

    pairs = parser["pairs"]
    
    if(pairs == 1)
        store_to_CSV_files(model, results, code, L, Npoints, U_values, V_values, is_only_pairs = true)
    else
        store_to_CSV_files(model, results, code, L, Npoints, U_values, V_values, is_only_pairs = false)
    end
end
