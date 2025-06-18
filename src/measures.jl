using LinearAlgebra
# Von Neumann entropy in bits = configuration space.
function von_neumann_entropy(rho; atol=1e-12)
    vals = eigen(Hermitian(rho)).values
    vals = vals[vals .> atol]
    return -sum(vals .* log.(vals)), -sum(vals .* log2.(vals))
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
#=
    Returns the density of up and down 
    electrons.
=#
function density_operators(N, psi, sites)
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
