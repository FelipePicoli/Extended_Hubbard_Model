#= 
    DMRG for various U and V for the exteded Hubbard model.
=#
#=
    Implement argparser inputs 
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

    # Correction to BOW phase which i'm not considering.
    epsilon = 10e-2 

    for (i, U) in enumerate(U_values)
        for (j, V) in enumerate(V_values)
            #=
                Add check if files with Uf and Vf already exists. 
            =#
            @show (i, U), (j, V)

            sites = siteinds("Electron", L; conserve_qns=true)
            H = H_EHM(L, J, U, V, sites)
            state = fill("Emp", L)

            if (abs(V + 1.5 * epsilon) >= 0 && abs(U + 1.5*epsilon) >= 0 || 
                (V + 1.5 * epsilon < 0 && U + 1.5*epsilon < 0))
                println("Metallic") 
                state = random_metallic_state(L, Nup, Ndn)
            elseif ((abs(V) >= abs(2 * U + epsilon)) || 
                (V > 0 && U < 0)) 
                state = random_cdw_state(L, Nup, Ndn)
                println("CDW") 
            elseif (abs(V) < abs(2 * U + epsilon)) 
                state = random_sdw_state(L, Nup, Ndn)
                println("SDW")
            else
                state = random_ps_state(L, Nup, Ndn)
                println("PS")
            end
            psi0 = random_mps(sites, state; linkdims=m)

            # Start DMRG calculation:
            energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
            #=
                Can reconstruct the densities from these three 
                results.
            =#
            upd, dnd, updn = density_operators(L, psi, sites)

            charge_density = upd .+ dnd
            # removed the 1/2 factor
            magnetization = (upd .- dnd) / 2
            
            # compute order parameters
            op_m_cdw = abs(m_cdw(L, charge_density))
            op_m_sdw = abs(m_sdw(L, magnetization))               
            #= 
                Implementation of the single-site entanglement 
                to be compared with the form of S = 1 - (1/L) \sum_{i} Tr (rho_i^2)
            =#
            avg_ss_entanglement = average_single_site_entanglement(L, upd, dnd, updn)

            # Reduced density matrix computations

            rho_1, N_bulk= build_1_particle_rdm(psi, upd, dnd)

            S_bulk, S_bulk_bits = von_neumann_entropy(rho_1)
            E_p_bulk, E_p_bulk_bits = S_bulk - log(N_bulk), S_bulk_bits - log2(N_bulk) 

            coh_1rdm = quantum_coherence(rho_1)

            rho_2 = build_2_particle_rdm(psi, sites)

            Q_2, Q_2_bits  = von_neumann_entropy(rho_2) 
            Q_2, Q_2_bits  = Q_2 - log(L*(L-1)/2), Q_2_bits - log2(L*(L-1)/2)

            # Check number of particles
            # @show flux(psi0)
            # @show flux(psi)
            coh_2rdm = quantum_coherence(rho_2)

            info_input = @sprintf("XXXXX_%s_L=%d_U=%.2f_V=%.2f_NPoints=%d.txt", model, L, U, V, Npoints)

            filename = joinpath(results, replace(info_input, "XXXXX" => "charge_density"))
            writedlm(filename, charge_density)

            filename = joinpath(results, replace(info_input, "XXXXX" => "magnetization"))
            writedlm(filename, magnetization)

            filename = joinpath(results, replace(info_input, "XXXXX" => "doublons"))
            writedlm(filename, updn)

            filename = joinpath(results, replace(info_input, "XXXXX" => "E_GS"))
            writedlm(filename, energy)

            filename = joinpath(results, replace(info_input, "XXXXX" => "E_p"))
            writedlm(filename, E_p_bulk)

            filename = joinpath(results, replace(info_input, "XXXXX" => "E_p_bits"))
            writedlm(filename, E_p_bulk_bits)

            filename = joinpath(results, replace(info_input, "XXXXX" => "S"))
            writedlm(filename, avg_ss_entanglement)

            filename = joinpath(results, replace(info_input, "XXXXX" => "coh_1rdm"))
            writedlm(filename, coh_1rdm)

            filename = joinpath(results, replace(info_input, "XXXXX" => "coh_2rdm"))
            writedlm(filename, coh_2rdm)

            filename = joinpath(results, replace(info_input, "XXXXX" => "Q_2"))
            writedlm(filename, Q_2)

            filename = joinpath(results, replace(info_input, "XXXXX" => "Q_2_bits"))
            writedlm(filename, Q_2_bits)

            filename = joinpath(results, replace(info_input, "XXXXX" => "m_sdw"))
            writedlm(filename, op_m_sdw)

            filename = joinpath(results, replace(info_input, "XXXXX" => "m_cdw"))
            writedlm(filename, op_m_cdw)
            #=      
                Free up memory. 
            =# 
            H = nothing
            psi0 = nothing
            psi = nothing
            rho_1 = nothing
            rho_2 = nothing
            GC.gc()
            end
        end
    filename = joinpath(results, "U_vals_NPoints=$(Npoints).txt")
    writedlm(filename, U_values)
    filename = joinpath(results, "V_vals_NPoints=$(Npoints).txt")
    writedlm(filename, V_values)
end
