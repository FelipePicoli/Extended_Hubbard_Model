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
    incluce("module_DMRG_States_GS.jl")
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
    maxdim = [50, 100, 200, 400, 800, 800, 1000, 1200, 1400, 1600, 1800]
    cutoff = [1E-14]

    Npart = floor(Int, L/2) 
    Nup = Npart + L % 2 
    Ndn = L - Nup 
    
    sites = siteinds("Electron", L; conserve_qns=true)

    H = H_EHM(L, J, U, V, sites)

    state = state_ehm_diagram(L, Nup, Ndn, U, V)

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

    single_site_entanglement = average_single_site_entanglement(L, upd, dnd, updn)

    # Reduced density matrix computations
    rho_1 = build_1_particle_rdm(psi)

    S, S_bits = von_neumann_entropy(rho_1)
    E_p, E_p_bits = S - log(L), S_bits - log2(L) 

    coh_1rdm = quantum_coherence(rho_1)

    rho_2 = build_2_particle_rdm(psi, sites)

    Q_2, Q_2_bits  = von_neumann_entropy(rho_2) 
    Q_2, Q_2_bits  = Q_2 - log(L*(L-1)/2), Q_2_bits - log2(L*(L-1)/2)

    # Check number of particles
    # @show flux(psi0)
    # @show flux(psi)
    coh_2rdm = quantum_coherence(rho_2)

    store_EHM_GS_measures_results(results, model, L, U, V, charge_density, 
                                                        magnetization,
                                                        updn,
                                                        energy,
                                                        E_p, 
                                                        E_p_bits,
                                                        single_site_entanglement,
                                                        coh_1rdm,
                                                        coh_2rdm,
                                                        Q_2, 
                                                        Q_2_bits,
                                                        op_m_sdw,
                                                        op_m_cdw, 
                                                        NPoints)
    H = nothing
    psi0 = nothing
    psi = nothing
    GC.gc()
end
