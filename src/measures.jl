# Trying to use the MPS form to compute the reduced density matrix 
# properties. 
# For now comparing it to the single-site entanglement computed from correlators.
function average_entanglement_entropy(psi, L)
    # Average entropy taken over all sites. 
    S_avg = 0.0
    # Over all sites j
    for j = 1:L 

        # change orthogonality center to j
        psi = orthogonalize(psi, j)

        # tensor at site j
        A = psi[j]

        # prime the physical index 
        A_dag = dag(A)
        prime!(A_dag, "Site")

        rdm = A * A_dag
        # Diagonalize
        D, U = eigen(rdm)

        # Compute von Neumann entropy safely
        S = 0.0
        for n=1:dim(D, 1)
            p = D[n,n] 
            p_real = real(p)
            if p_real > 1e-12
                S -= p_real * log(p_real)
            end
        end
        S_avg += S 
    end
    return S_avg / L 
end
#=
    L - chain size 
    N - Dou
=#
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
#=
    Returns the density of up and down 
    electrons.
=#
function density_operators(N, psi)
    upd = fill(0.0, N)
    dnd = fill(0.0, N)
    updn = fill(0.0, N) 
    for j in 1:N
    orthogonalize!(psi, j)
    psidag_j = dag(prime(psi[j], "Site"))
    upd[j] = scalar(psidag_j * op(sites, "Nup", j) * psi[j])
    dnd[j] = scalar(psidag_j * op(sites, "Ndn", j) * psi[j])
    updn[j] = scalar(psidag_j * op(sites, "Nupdn", j) * psi[j])
    end
    return upd, dnd, updn
end
