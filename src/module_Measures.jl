include("module_Operators.jl")
using LinearAlgebra
# Von Neumann entropy in bits = configuration space.
function von_neumann_entropy(rho; atol=1e-12)
    vals = eigen(Hermitian(rho)).values
    vals = vals[vals .> atol]
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
function entanglement_gap(rho; cutoff=1e-12)
    rho = Hermitian(rho)
    eig = eigen(rho)
    eigvals = eig.values

    sorted_indices = sortperm(eigvals, rev=true)
    eigvals_sorted = eigvals[sorted_indices]

    eigvals_clipped = max.(eigvals_sorted, cutoff)
    entanglement_spectrum = -log.(eigvals_clipped)
    gaps = diff(entanglement_spectrum)

    return gaps, entanglement_spectrum
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
        "doublons" => updn,
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
