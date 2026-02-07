include("module_Fermionic_Operators.jl")
using LinearAlgebra

function single_site_entanglement(L, i, up, dn, updn)
    single_site_entanglement = 0.0
    w_2 = updn[i]
    w_up = up[i] - w_2
    w_dn = dn[i] - w_2
    w_0 = 1 - w_up - w_dn - w_2
    single_site_entanglement = 1 - (w_2^2 + w_up^2 + w_dn^2 + w_0^2)
    return single_site_entanglement
end

function compute_particle_rdm_quantities(L, sites, psi, dict_results, upd, dnd, updn; rdm="1rdm", cutoff=1e-12)

    # Reduced density matrix computations
    #
    shift = 0.0
    if occursin("1rdm", rdm)
        rho = get_1_particle_rdm(state, sites)
        shift = L
    elseif occursin("2rdm", rdm)
        rho = get_2_particle_rdm(psi, sites)
        shift = L * (L-1.0) / 2.0
    end

    # diagonalize rho
    eigenvalues = eigen(Hermitian(rho)).values
    eigenvalues = eigenvalues[eigenvalues .> cutoff]

    # Compute entanglement spectrum
    sorted_indices = sortperm(eigenvalues, rev=true)
    entanglement_spectrum = eigenvalues[sorted_indices]

    # Shifted Von Neumann entropy
    von_neumann, von_neumann_bits = - sum(eigenvalues .* log.(eigenvalues)), - sum(eigenvalues.* log2.(eigenvalues))
    von_neumann, von_neumann_bits = von_neumann - log(shift), von_neumann_bits - log2(shift)
    # quantum coherence
    quantum_coherence = sum(abs, rho) - sum(abs, diag(rho))


    # store to dict.
end
