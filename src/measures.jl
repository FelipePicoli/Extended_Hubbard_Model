using LinearAlgebra
# Von Neumann entropy in bits - configuration space.#
function S(rdm, L)

    println("trace (rho) = ", tr(rdm))
    println("ishermitian(rho) = ", ishermitian(rdm))

    lambdas = eigvals(rdm)
    S_bit = 0
    S = 0

    #=
    println("eigenvals = ", lambdas)
    println("eigenvals / L= ", lambdas / L)
    println("sum(eigvals) = ", sum(real(lambdas)))
    =#

    for lambda in lambdas
        lambda = real(lambda) 
        # S -= lambda * log2(lambda)
        # CRITICAL: Handle lambda <= 0 for log2
        if lambda > 1e-15 # A small threshold to avoid errors with log2
            S_bit -= lambda * log2(lambda)
            S -= lambda * log(lambda)
        end
    end 
    return S, S_bit
end

function quantum_coherence(rho, L)
    C = 0.0
    for i in 1:(2*L)
        for j in (i+1):(2*L)
            C += abs(rho[i, j])
        end
    end
    return C 
end
#=
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


function m_sdw(L, Sj)
    m_sdw = 0.0 
    for j in 1:L 
        m_sdw += (-1)^(j) * Sj[j]
    end
    return m_sdw / L
end
function m_cdw(L, nj)
    m_cdw = 0 
    for j in 1:L 
        m_cdw += (-1)^(j) * (nj[j] - 1)
    end
    return m_cdw / L 
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
