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
    include("module_States_Ansatze.jl")
    include("module_Argparse.jl")
    include("module_Data.jl")

    parser = parse_commandline()

    # model parameters
    L = parser["L"]
    J = parser["J"]
    U0 = parser["U0"]
    V0 = parser["V0"]
    Uf = parser["Uf"]
    Vf = parser["Vf"]

    Npoints = parser["Npoints"]
    output_name = parser["output_name"]

    U_range = range(U0, Uf, length=Npoints)
    V_range = range(V0, Vf, length=Npoints)

    # Force the U, V values to have two decimal places.
    model = parser["model"]
    results_path = parser["results"]
    preprocessing_path = parser["preprocessing_path"]
    random_mps_path = joinpath(preprocessing_path, parser["previous_random_mps"])

    nsweeps = parser["nsweeps"]
    m = parser["m"]
    maxdim = [10, 20, 50, 100, 200, 400, 800]
    cutoff = [1E-14]

    Npart = floor(Int, L/2)
    Nup = Npart + L % 2
    Ndn = L - Nup

    sites = siteinds("Electron", L; conserve_qns=true)
    H_J, H_U, H_V = get_EHM_Hamiltonian_components(L, sites)
    psi0 = nothing
    for U in U_range, V in V_range

        results = data_structure(L)

        U_str, V_str = @sprintf("%.2f", U), @sprintf("%.2f", V)

        info_output = "$(output_name)_$(model)_L=$(L)_U=$(U_str)_V=$(V_str)_Npoints=$(Npoints).jld2"
        file_results = joinpath(results_path, info_output)

        if isfile(file_results)
            println("Simulation for U = $(U_str) and V = $(V_str) already run. Skipping.")
        else
            println("Running for U = $(U_str) and V = $(V_str).")

            println("Building Hamiltonian")
            H_hamiltonian = J*H_J + U*H_U + V*H_V
            truncate!(H_hamiltonian; cutoff=1e-15)

            # Use the resulting DMRG psi as next guess.
            if psi0 == nothing
                state = state_ehm_diagram(L, Nup, Ndn, U, V)
                psi0 = random_mps(sites, state; linkdims=m)
            end

            println("Runnning DMRG")
            energy, psi = dmrg(H_hamiltonian, psi0; nsweeps, maxdim=maxdim, cutoff=cutoff)

            results["scalar"]["energy_GS"] = energy
            # @show energy
            #=
                Must use the siteinds to obtain the correct ID.
            =#
            upd, dnd, updn = density_operators(L, psi, siteinds(psi))

            charge_density = upd .+ dnd
            magnetization = (upd .- dnd) / 2

            for i in 1:L
                results["vector"]["single_site_entanglement_GS"][i] = single_site_entanglement(L, i, upd, dnd, updn)
            end

            compute_particle_rdm_quantities(L, siteinds(psi), psi, results, upd, dnd, updn; rdm="1rdm")
            compute_particle_rdm_quantities(L, siteinds(psi), psi, results, upd, dnd, updn; rdm="2rdm")

            results["vector"]["magnetization_GS"] = magnetization
            results["vector"]["charge_density_GS"] = charge_density
            results["vector"]["doublons_GS"] = updn

            store_results(results, info_output, results_path)
            psi0 = psi
        end
        GC.gc()
    end
end

