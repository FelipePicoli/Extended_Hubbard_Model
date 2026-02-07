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
    shift = 0.0
    if occursin("1rdm", rdm)
        rho = get_1_particle_rdm(psi, sites)
        shift = L
    elseif occursin("2rdm", rdm)
        rho = get_2_particle_rdm(psi, sites)
        shift = L * (L-1.0) / 2.0
    end
    # diagonalize rho
    eigenvalues = eigen(Hermitian(rho)).values
    eigenvalues = max.(0.0, eigenvalues)

    eigenvals_sorted = sort(eigenvalues, rev=true)
    eigenvals_sorted = max.(cutoff, eigenvals_sorted)
    entanglement_spectrum = -log.(eigenvals_sorted)

    eigenvalues = eigenvalues[eigenvalues .> cutoff]
    if size(eigenvalues) == 0
        println("Error - no eigenvalues.")
        return 0.0
    end

    # Shifted Von Neumann entropy
    von_neumann = - sum(eigenvalues .* log.(eigenvalues))
    von_neumann_shifted = von_neumann - log(shift)

    # quantum coherence
    quantum_coherence = sum(abs, rho) - sum(abs, diag(rho))
    dict_results["scalar"]["von_neumann_$(rdm)_GS"] = von_neumann_shifted
    dict_results["scalar"]["quantum_coherence_$(rdm)_GS"] = quantum_coherence
    dict_results["vector"]["entanglement_spectrum_$(rdm)_GS"] = entanglement_spectrum
end
