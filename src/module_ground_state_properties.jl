function compute_GS_measures(L, sites, psi)
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
        Implementation of average single-site entanglement 
            S = 1 - (1/L) \sum_{i} Tr (rho_i^2)
    =#
    single_site_entanglement = average_single_site_entanglement(L, upd, dnd, updn)

    # Reduced density matrix computations
    rho_1 = build_1_particle_rdm(psi)

    S, S_bits = von_neumann_entropy(rho_1)
    E_p, E_p_bits = S - log(L), S_bits - log2(L) 

    coh_1rdm = quantum_coherence(rho_1)

    gaps_1rdm, xis_1rdm = entanglement_gap(rho_1) 

    rho_2 = build_2_particle_rdm(psi, sites)

    Q_2, Q_2_bits  = von_neumann_entropy(rho_2) 
    Q_2, Q_2_bits  = Q_2 - log(L*(L-1)/2), Q_2_bits - log2(L*(L-1)/2)

    # Check number of particles
    # @show flux(psi0)
    # @show flux(psi)
    coh_2rdm = quantum_coherence(rho_2)

    gaps_2rdm, xis_2rdm = entanglement_gap(rho_2) 
    
    dict = Dict(
        "charge_density" => charge_density, 
        "magnetization" => magnetization,
        "doublon" => updn,
        "energy" => energy,
        "E_p" => E_p, 
        "E_p_bits" => E_p_bits,
        "single_site_entanglement" => single_site_entanglement,
        "coh_1rdm" => coh_1rdm,
        "gaps_1rdm" => gaps_1rdm, 
        "xis_1rdm" => xis_1rdm,
        "coh_2rdm" => coh_2rdm,
        "Q_2" => Q_2, 
        "Q_2_bits" => Q_2_bits,
        "gaps_2rdm" => gaps_2rdm, 
        "xis_2rdm" => xis_2rdm,
        "op_m_sdw" => op_m_sdw,
        "op_m_cdw" => op_m_cdw, 
    ) 
    return dict
end

