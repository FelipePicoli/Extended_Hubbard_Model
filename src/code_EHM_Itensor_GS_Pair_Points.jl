#= 
    DMRG for various U and V for the exteded Hubbard model.
=#
using ITensors
using ITensorMPS

using Random

using ArgParse

using DelimitedFiles
using DataFrames
using CSV
using Printf

using PrettyTables


let
    include("operators.jl")
    include("measures.jl")
    include("parser.jl")

    parser = parse_commandline() 

    # model parameters
    L = parser["L"]
    model = parser["model"]
    results = parser["results"]
    J = parser["J"]
    U0 = parser["U0"]
    Uf = parser["Uf"]
    V0 = parser["V0"]
    Vf = parser["Vf"]

    Npoints = parser["Npoints"]

    U_values = range(U0, stop=Uf, length=Npoints)
    V_values = range(V0, stop=Vf, length=Npoints)

    # dmrg parameters 
    nsweeps = parser["nsweeps"]
    m = parser["m"]
    maxdim = [50, 100, 200, 400, 800, 800, 1000, 1200, 1400]
    cutoff = [1E-14]

    Npart = floor(Int, L/2) 
    Nup = Npart + L % 2 
    Ndn = L - Nup 
    
    sites = siteinds("Electron", L; conserve_qns=true)

    for (U, V) in zip(U_values, V_values)
        compute_measures_single_point(sites, Nup, Ndn, U, V, nsweeps, maxdim, cutoff)
    end
end
