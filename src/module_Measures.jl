include("module_Operators.jl")
using LinearAlgebra
# Von Neumann entropy in bits = configuration space.
function von_neumann_entropy(vals; atol=1e-12)
    return - sum(vals .* log.(vals)), - sum(vals .* log2.(vals))
end

function quantum_coherence(rho)
    dim = size(rho)[1]
    C = 0.0
    for i in 1:dim
        for j in (i+1):dim
            # maybe return here and force hermitian()
            C += abs(rho[i, j])
        end
    end
    return C 
end

#=
    Diagonalizes the matrix rho as 
    rho = \sum_i e^{-\xi_i} with x_{i+1} >= x_i 
    and returns the difference from the 
=#
function entanglement_gap(eigenvals; cutoff=1e-12)

    sorted_indices = sortperm(eigvals, rev=true)
    eigvals_sorted = eigvals[sorted_indices]

    eigvals_clipped = max.(eigvals_sorted, cutoff)
    entanglement_spectrum = -log.(eigvals_clipped)
    gaps = diff(entanglement_spectrum)

    return gaps, entanglement_spectrum
end

function entanglement_gap(eigvals; cutoff=1e-12)

    sorted_indices = sortperm(eigvals, rev=true)
    eigvals_sorted = eigvals[sorted_indices]
    eigvals_clipped = max.(eigvals_sorted, cutoff)

    xis = -log.(eigvals_clipped)

    return xis
end


using Statistics
function average_single_site_entanglement(L, up, dn, updn)
    single_site_entanglement = fill(0.0, L)
    for i in 1:L 
        w_2 = updn[i] 
        w_up = up[i] - w_2 
        w_dn = dn[i] - w_2 
        w_0 = 1 - w_up - w_dn - w_2 
        single_site_entanglement[i] = 1 - (w_2^2 + w_up^2 + w_dn^2 + w_0^2)
    end
    return Statistics.mean(single_site_entanglement)
end 

function m_sdw(L, Sj)
    m_sdw_val = 0.0 
    for j in 1:L 
        m_sdw_val += (-1)^(j-1) * Sj[j]
    end
    return m_sdw_val / L
end
function m_cdw(L, nj)
    m_cdw_val = 0 
    for j in 1:L 
        m_cdw_val += (-1)^(j-1) * (nj[j] - 1)
    end
    return m_cdw_val / L 
end 

function compute_GS_measures(L, sites, psi; cutoff=1e-12)
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
            \mathcal{L} = 1 - (1/L) \sum_{i} Tr (rho_i^2)
    =#
    single_site_entanglement = average_single_site_entanglement(L, upd, dnd, updn)

    # Reduced density matrix computations
    

    # one-particle RDM
    
    rho_1 = build_1_particle_rdm(psi)

    # diagonalize rho_1 
    vals_1 = eigen(Hermitian(rho_1)).values

    # Compute entanglement spectrum 
    sorted_indices = sortperm(vals_1, rev=true)
    eigvals_sorted = vals_1[sorted_indices]
    eigvals_clipped = max.(eigvals_sorted, cutoff)
    Omega_1rdm = -log.(eigvals_clipped)

    # Shifted Von Neumann entropy 
    pos_vals = vals_1[vals_1 .> cutoff]
    S_1rdm, S_1rdm_bits = - sum(vals_1 .* log.(vals_1)), - sum(vals_1 .* log2.(vals_1))
   
    E_p, E_p_bits = S_1rdm - log(L), S_1rdm_bits - log2(L)
    
    # quantum coherence
    coh_1rdm = sum(abs, rho_1) - sum(abs, diag(rho_1))

    # two-particle RDM
    rho_2 = build_2_particle_rdm(psi, sites)

    # diagonalize rho_2
    vals_2 = eigen(Hermitian(rho_2)).values

    # Compute entanglement spectrum 
    sorted_indices = sortperm(vals_2, rev=true)
    eigvals_sorted = vals_2[sorted_indices]
    eigvals_clipped = max.(eigvals_sorted, cutoff)
    Omega_2rdm = -log.(eigvals_clipped)
    
    # Shifted Von Neumann
    pos_vals = vals_2[vals_2 .> cutoff]
    S_2rdm, S_2rdm_bits = - sum(vals_2 .* log.(vals_2)), - sum(vals_2 .* log2.(vals_2))

    Q_2, Q_2_bits  = S_2rdm - log(L*(L-1)/2), S_2rdm_bits - log2(L*(L-1)/2)

    # quantum coherence
    coh_2rdm = sum(abs, rho_2) - sum(abs, diag(rho_2))

    dict = Dict(
        "charge_density" => charge_density, 
        "magnetization" => magnetization,
        "doublons" => updn,
        "E_p" => E_p, 
        "E_p_bits" => E_p_bits,
        "single_site_entanglement" => single_site_entanglement,
        "coh_1rdm" => coh_1rdm,
        "Omega_1rdm" => Omega_1rdm,
        "coh_2rdm" => coh_2rdm,
        "Q_2" => Q_2, 
        "Q_2_bits" => Q_2_bits,
        "Omega_2rdm" => Omega_2rdm,
        "op_m_sdw" => op_m_sdw,
        "op_m_cdw" => op_m_cdw, 
    ) 
    return dict
end
