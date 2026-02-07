#=
    Various functions for creation of operators in ItensorMPS.jl.
    - Extended Hubbard model Hamiltonian.
    - Density operators
    - Random ansatz for initial states in phases such as
        - Metallic
        - Charge density wave
        - Spin density wave
=#
using ITensors
using ITensorMPS
using Random
using BlockDiagonals

#=
    Generates the MPO for the EHM Hamiltonian
    with strengths J, U and V.
    Requires a SiteType sites.
=#
function H_EHM(N, J, U, V, sites)
    os = OpSum()
    for i in 1:(N - 1)
      # Hopping
      os -= J, "Cdagup", i, "Cup", i + 1
      os -= J, "Cdagup", i + 1, "Cup", i
      os -= J, "Cdagdn", i, "Cdn", i + 1
      os -= J, "Cdagdn", i + 1, "Cdn", i

      # Nearest-neighbours intearction
      os += V, "Ntot", i, "Ntot", i + 1
    end
    # On-site Coulomb interaction
    for i in 1:N
      os += U, "Nupdn", i
    end
    # Trying to use the cutoff to optimize things.
    return MPO(os, sites; cutoff=1e-12)
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
#=
    Builds the 1-particle reduced density matrix.

    This matrix depends on just two correlators so it is
    easily implemented by the correlation_matrix from ITensorMPS.
=#
function get_1_particle_rdm(state, sites)
    L = length(sites)
    #=
        N and not L since any site can have spin up or down
    =#
    rho_1 = zeros(ComplexF64, 2*L, 2*L)

    Cupup = correlation_matrix(state, "Cdagup", "Cup")
    Cdndn = correlation_matrix(state, "Cdagdn", "Cdn")
    Cupdn = correlation_matrix(state, "Cdagup", "Cdn")

    for i in 1:L
        for j in 1:L
            rho_1[i, j] = Cupup[i,j]
            rho_1[i, j + L] = Cupdn[i,j]
            rho_1[i + L, j] = Cupdn[i, j]'
            rho_1[i + L, j + L] = Cdndn[i,j]
        end
    end
    return rho_1 = rho_1 / L
end

## Methods for building the two-particle reduced density matrix.
function get_all_modes_pairs(modes)
    return [(m1, m2) for m2 in modes for m1 in modes]
end

function select_unique_pairs(sites_pairs; ordered = false)
    if ordered
        # return [(i, j) for (j, i) in sites_pairs if id(i) < id(j)]
        return [(i, j) for (j, i) in sites_pairs if i < j]
    else
        return unique(sites_pairs)
    end
end

function handler_disjoint_spins(all_mode_pairs)

    pairs = []
    for (mode_m, mode_n) in all_mode_pairs

        m_site, m_spin = mode_m
        n_site, n_spin = mode_n

        if m_spin != n_spin
            push!(pairs, (m_site, n_site))
        end
    end
    return select_unique_pairs(pairs, ordered = false)
end
function handler_joint_spins(all_mode_pairs)

    pairs_spin_up = []
    pairs_spin_dn = []

    for (mode_m, mode_n) in all_mode_pairs

        m_site, m_spin = mode_m
        n_site, n_spin = mode_n

        if m_site != n_site && m_spin == n_spin
            target = (m_spin == "up") ? pairs_spin_up : pairs_spin_dn
            push!(target, (m_site, n_site))
        end
    end
    return [select_unique_pairs(pairs_spin_up, ordered = true),
            select_unique_pairs(pairs_spin_dn, ordered = true)]
end
function two_fermions_basis_sites(modes; disjoint_spins = true)
    all_mode_pairs = get_all_modes_pairs(modes)
    if disjoint_spins
        return handler_disjoint_spins(all_mode_pairs)
    end
    return handler_joint_spins(all_mode_pairs)
end
#=
    parameters:
    phi - state
    sites - selected sites sorted

    Returns a two-particle reduced density matrix.

    This function computes only the part of the matrix in which
    all of the correlators are different.

    Tolerance of the julia language packages are very low, so in general
    this computation gives various non-hermitian matrices.
=#
function compute_2rdm_block(phi, sites, pairs, spins)

    dim = length(pairs)
    rho_rdm = zeros(ComplexF64, dim, dim)

    s1, s2, s3, s4 = spins

    # site_number(site_index) = parse(Int, match(r"n=(\d+)", string(tags(site_index))).captures[1])

    for p in 1:dim
        i, j = pairs[p]
        for q in p:dim
            k, l = pairs[q]

            # @show i, j, k, l
            # @show s1, s2, s3, s4

            os = OpSum()
            os -= "Cdag$s1", i,
                  "Cdag$s2", j,
                  "C$s3", k,
                  "C$s4", l

            O_mpo = MPO(os, sites; cutoff=1e-15)

            val = inner(phi', O_mpo, phi)
            rho_rdm[p, q] = val
            rho_rdm[q, p] = conj(val)
        end
    end
    return rho_rdm
end
# Create list of all modes (site, spin) for selected sites
function site_spin_modes(sites)
    site_number(site_index) = parse(Int, match(r"n=(\d+)", string(tags(site_index))).captures[1])
    return [(site_number(site), spin) for site in sites for spin in ("up", "dn")]
end
function get_2_particle_RDM(phi, sites)
    L = length(sites)

    modes = site_spin_modes(sites)

    sites_pairs_equal_spins_sector = two_fermions_basis_sites(modes, disjoint_spins = false)

    pairs_spin_up = sites_pairs_equal_spins_sector[1]
    pairs_spin_dn = sites_pairs_equal_spins_sector[2]

    pairs_mixed_spins = two_fermions_basis_sites(modes, disjoint_spins = true)

    rho_upup = compute_2rdm_block(phi, sites, pairs_spin_up, ["up","up", "up", "up"])
    rho_dndn = compute_2rdm_block(phi, sites, pairs_spin_dn, ["dn","dn", "dn", "dn"])
    rho_updn = compute_2rdm_block(phi, sites, pairs_mixed_spins, ["up","dn", "up", "dn"])

    rho_2 = BlockDiagonal([rho_upup, rho_dndn, rho_updn])
    rho_2 *= (2.0 / (L*(L-1)))
    return rho_2
end
