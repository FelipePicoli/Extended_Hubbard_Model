#=
    DMRG for various U and V for the exteded Hubbard model.
=#
using ITensors
using ITensorMPS
using Random
using ArgParse
using DelimitedFiles
using DataFrames
using Printf
using LinearAlgebra
using JLD2
let
    include("module_Fermionic_Operators.jl")
    include("module_Quantum_Info_Tools.jl")
    include("module_Preproccessing.jl")
    include("module_Argparse.jl")
    include("module_Data.jl")

    parser = parse_commandline()

    # model parameters
    L = parser["L"]
    J = parser["J"]
    U = parser["U0"]
    V = parser["V0"]

    model = parser["model"]
    results_path = parser["results"]
    result_file_name = parser["result_file_name"]

    paths = (
            sites            = joinpath(results_path, parser["site_inds_path"]),
            hamiltonian_mpos = joinpath(results_path, parser["hamiltonian_mpos"]),
            random_mps       = joinpath(results_path, parser["previous_random_mps"]),
            results          = results_path
        )

    # Setting the results data structure.
    results = data_structure(L)

    nsweeps = parser["nsweeps"]
    m = parser["m"]
    maxdim = [10, 20, 50, 100, 200, 400, 800]
    cutoff = [1E-14]

    Npart = floor(Int, L/2)
    Nup = Npart + L % 2
    Ndn = L - Nup

    H_hamiltonian, psi0, sites = preprocessing_simulation(L, J, U, V, Nup, Ndn, m, paths)

    println("Runnning DMRG")
    energy, psi = dmrg(H_hamiltonian, psi0; nsweeps, maxdim=maxdim, cutoff=cutoff)

    results["scalar"]["energy_GS"] = energy
    @show energy

    upd, dnd, updn = density_operators(L, psi, sites)

    charge_density = upd .+ dnd
    magnetization = (upd .- dnd) / 2

    for i in 1:L
        results["vector"]["single_site_entanglement_GS"][i] = single_site_entanglement(L, i, upd, dnd, updn)
    end

    compute_particle_rdm_quantities(L, sites, psi, results, upd, dnd, updn; rdm="1rdm")
    compute_particle_rdm_quantities(L, sites, psi, results, upd, dnd, updn; rdm="2rdm")

    results["vector"]["magnetization_GS"] = magnetization
    results["vector"]["charge_density_GS"] = charge_density
    results["vector"]["doublons_GS"] = updn

    store_results(results, result_file_name, paths.results)

    # Store the random_MPS to .jld2 file for next run.
    jldsave(paths.random_mps; psi)
    H = nothing
    psi0 = nothing
    GC.gc()
end

