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
using LinearAlgebra

let
    include("module_Fermionic_Operators.jl")
    include("module_Quantum_Info_Tools.jl")
    include("module_States_Ansatze.jl")
    include("module_Argparse.jl")
    include("module_Data.jl")

    parser = parse_commandline()

    # model parameters
    L = 4 # parser["L"]
    model = parser["model"]
    results = parser["results"]
    J = parser["J"]
    U = 0.0 # parser["U0"]
    V = 0.0 # parser["V0"]

    @show (L, J, U, V)

    # Setting the results data structure.
    results = Dict(
        "scalar" =>         Dict(
                                "energy_GS" => 0.0,
                                "von_neumann_1rdm_GS" => 0.0,
                                "quantum_coherence_1rdm_GS" => 0.0,
                                "von_neumann_2rdm_GS" => 0.0,
                                "quantum_coherence_2rdm_GS" => 0.0,
                            ),
        "vector" =>         Dict(
                                "single_site_entanglement_GS" => zeros(L),
                                "charge_density_GS" => zeros(L),
                                "magnetization_GS" => zeros(L),
                                "doublons_GS" => zeros(L),
                            ),
        "spectrum" =>       Dict(
                                "entanglement_spectrum_1rdm_GS" => zeros(1, 2*L),
                                "entanglement_spectrum_2rdm_GS" => zeros(1, L*(2*L-1)),
                            )
    )
    nsweeps = parser["nsweeps"]
    m = parser["m"]
    maxdim = [10, 20, 50, 100, 200, 400, 800]
    cutoff = [1E-14]
    Npart = floor(Int, L/2)
    Nup = Npart + L % 2
    Ndn = L - Nup

    sites = siteinds("Electron", L; conserve_qns=true)

    # println("Building H")
    H = H_EHM(L, J, U, V, sites)

    state = state_ehm_diagram(L, Nup, Ndn, U, V)

    # println("Obtaining random MPS")
    psi0 = random_mps(sites, state; linkdims=m)

    # println("Runnning DMRG")
    # Start DMRG calculation:
    energy, psi = dmrg(H, psi0; nsweeps, maxdim=maxdim, cutoff=cutoff)

    results["scalar"]["energy_GS"] = energy

    upd, dnd, updn = density_operators(L, psi, sites)

    charge_density = upd .+ dnd
    magnetization = (upd .- dnd) / 2

    for i in 1:L
        results["vector"]["single_site_entanglement_GS"][i] = single_site_entanglement(L, i, upd, dnd, updn)
    end

    @show results["vector"]["single_site_entanglement_GS"]

    rho2 = get_2_particle_RDM(psi, sites)

    compute_particle_rdm_quantities(L, sites, psi, results, upd, dnd, updn; rdm="1rdm")
    # compute_particle_rdm_quantities(L, sites, psi, results, upd, dnd, updn; rdm="2rdm")

    results["vector"]["magnetization_GS"] = magnetization
    results["vector"]["charge_density_GS"] = charge_density
    results["vector"]["doublons_GS"] = updn

    # store_EHM_GS_measures_results(results, model, L, U, V, dict_results)
    H = nothing
    psi0 = nothing
    GC.gc()
end
