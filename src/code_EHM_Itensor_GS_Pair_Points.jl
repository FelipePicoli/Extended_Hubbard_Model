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

let
    include("module_Operators.jl")
    include("module_Measures.jl")
    include("module_Argparse.jl")
    include("module_DMRG_States_GS.jl")
    include("module_Organize_Results.jl")

    parser = parse_commandline() 

    # model parameters
    L = parser["L"]
    model = parser["model"]
    results = parser["results"]
    J = parser["J"]
    U = parser["U0"]
    V = parser["V0"]
    
    @show (U, V) 

    Npoints = parser["Npoints"]

    # dmrg parameters 
    nsweeps = parser["nsweeps"]
    m = parser["m"]
    maxdim = [50, 100, 200, 400, 800, 800, 1000, 1200]
    cutoff = [1E-14]
    Npart = floor(Int, L/2) 
    Nup = Npart + L % 2 
    Ndn = L - Nup 
    
    sites = siteinds("Electron", L; conserve_qns=true)
    
    println("Building H")
    H = H_EHM(L, J, U, V, sites)

    state = state_ehm_diagram(L, Nup, Ndn, U, V)
    
    println("Obtaining random MPS")
    psi0 = random_mps(sites, state; linkdims=m)
    
    println("Runnning DMRG")
    # Start DMRG calculation:
    energy, psi = dmrg(H, psi0; nsweeps, maxdim=maxdim, cutoff=cutoff)

    dict_results = compute_GS_measures(L, sites, psi)
    dict_results["energy"]  = energy

    upd, dnd, updn = density_operators(L, psi, sites)
    
    @show upd 
    @show dnd 
    @show updn

    # store_EHM_GS_measures_results(results, model, L, U, V, Npoints, dict_results)

    H = nothing
    psi0 = nothing
    psi = nothing
    GC.gc()
end
