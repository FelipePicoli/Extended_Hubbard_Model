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
    include("module_Fermionic_Operators.jl")
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
    maxdim = [50, 100, 200, 400, 800, 800, 1000, 1200, 1400, 1600, 1800]
    cutoff = [1E-8]

    Npart = floor(Int, L/2) 
    Nup = Npart + L % 2 
    Ndn = L - Nup 

    # Correction to BOW phase which i'm not considering.
    # Runing small values for now...
    sites = siteinds("Electron", L; conserve_qns=true)

    for (i, U) in enumerate(U_values)
        for (j, V) in enumerate(V_values)

            @show U, V 
            H = H_EHM(L, J, U, V, sites)
            #= 
                The best state for variational step of the DMRG algorithm.
            =#
            state = state_ehm_diagram(L, Nup, Ndn, U, V)
            psi0 = random_mps(sites, state; linkdims=m)

            # Start DMRG calculation:
            energy, psi = dmrg(H, psi0; nsweeps, maxdim=maxdim, cutoff=cutoff)

            dict_results = compute_GS_measures(L, sites, psi)
            dict_results["energy"]  = energy

            store_EHM_GS_measures_results(results, model, L, U, V, Npoints, dict_results)

            # Free up memory 
            H = nothing
            psi0 = nothing
            psi = nothing
            GC.gc()
        end
    end
end
