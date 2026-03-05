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

using PrettyTables

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
    Generates the MPOs for the EHM Hamiltonian
        H_J, H_U, U_V
    One can build the full Hamiltonian with
    potentials as
    H = J * H_J + U * H_U + V * H_V
    and adding
    truncate!(H; cutoff= eps)
    increases performance.
=#
function get_EHM_Hamiltonian_components(N, sites)
    # Hopping term
    os_J = OpSum()
    for i in 1:(N - 1)
        os_J -= 1.0, "Cdagup", i, "Cup", i + 1
        os_J -= 1.0, "Cdagup", i + 1, "Cup", i
        os_J -= 1.0, "Cdagdn", i, "Cdn", i + 1
        os_J -= 1.0, "Cdagdn", i + 1, "Cdn", i
    end
    H_J = MPO(os_J, sites)

    # On-site Coulomb interaction
    os_U = OpSum()
    for i in 1:N
        os_U += 1.0, "Nupdn", i
    end
    H_U = MPO(os_U, sites)

    # Nearest-neighbor Interaction
    os_V = OpSum()
    for i in 1:(N - 1)
        os_V += 1.0, "Ntot", i, "Ntot", i + 1
    end
    H_V = MPO(os_V, sites)

    return H_J, H_U, H_V
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
function get_1_particle_rdm(psi, sites)
    L = length(sites)
    #=
        N and not L since any site can have spin up or down
    =#
    rho_1 = zeros(ComplexF64, 2*L, 2*L)

    rho_upup = correlation_matrix(psi, "Cdagup", "Cup")
    rho_dndn = correlation_matrix(psi, "Cdagdn", "Cdn")
    rho_updn = correlation_matrix(psi, "Cdagup", "Cdn")

    for i in 1:L
        for j in 1:L
            rho_1[i, j] = rho_upup[i,j]
            rho_1[i, j + L] = rho_updn[i,j]
            rho_1[i + L, j] = rho_updn[i, j]'
            rho_1[i + L, j + L] = rho_dndn[i,j]
        end
    end
    return rho_1 = rho_1 / L
end

#=
    Given an initial MPS phi and list of sites, returns
    the two-particle reduced density matrix in block diagonal form.
=#
function get_2_particle_rdm(phi, sites)

    L = length(sites)
    pairs_equal_spins, pairs_mixed_spins = get_pairs_spins_sites(sites)

    @show pairs_equal_spins
    @show pairs_mixed_spins

    rho_upup = compute_block(phi, sites, pairs_equal_spins, ["up", "up"])
    # If L is even, then both blocks are equal. Still a must check.
    # rho_dndn = (L%2 == 0) ? rho_upup : compute_block(phi, sites, pairs_equal_spins, ["dn", "dn"])
    rho_dndn = compute_block(phi, sites, pairs_equal_spins, ["dn", "dn"])
    rho_updn = compute_block(phi, sites, pairs_mixed_spins, ["up", "dn"])

    @show size(rho_upup)
    pretty_table(rho_upup)
    @show size(rho_dndn)
    pretty_table(rho_dndn)
    @show size(rho_updn)
    pretty_table(rho_updn)
    # rho_nup_nup = correlation_matrix(phi, "Nup", "Nup")
    # pretty_table(rho_nup_nup)
    rho_2 = BlockDiagonal([rho_upup, rho_dndn, rho_updn])
    @show size(rho_2)
    rho_2 *= (2.0 / (L*(L-1)))
    return rho_2
end
#=
    Builds a block of the two particle reduced density matrix given an initial
    MPS phi, a list of sites, pairs - can be mixed or same spins, and which spins.
    The spins variables can be either
    ["up", "up"], ["dn", "dn"],["up", "dn"]
=#
function compute_block(phi, sites, pairs, spins)
    dim = length(pairs)
    rho = zeros(ComplexF64, dim, dim)

    s1, s2 = spins

    for p in 1:dim
        i, j = pairs[p]
        for q in p:dim
            k, l = pairs[q]

            os = OpSum()
            os -= "Cdag$s1", i, "Cdag$s2", j, "C$s1", k, "C$s2", l

            O_mpo = MPO(os, sites) # ; cutoff=1e-15)
            val = inner(phi', O_mpo, phi)

            rho[p, q] = val
            rho[q, p] = conj(val)
        end
    end
    return rho
end
#=
    Given a list of sites returns the unique list of pairs of sites for same or mixed spins blocks
    of the rho_2rdm.
=#
function get_pairs_spins_sites(sites)

    # Create list of all modes (site, spin) for selected sites
    site_number(site_index) = parse(Int, match(r"n=(\d+)", string(tags(site_index))).captures[1])
    modes = [(site_number(site), spin) for site in sites for spin in ("up", "dn")]

    all_mode_pairs = [(m1, m2) for m2 in modes for m1 in modes]

    # Combinations of pairs of sites for equal spins and mixed spins.
    # Pairs of same spins corresponds to blocks of the matrix where i < j, k < l
    # the constraint does not apply to the mixed spins case.
    pairs_same_spins = []
    pairs_mixed_spins = []

    for (mode_m, mode_n) in all_mode_pairs
        m_site, m_spin = mode_m
        n_site, n_spin = mode_n
        # Case of same spins, different sites
        if m_site != n_site && m_spin == n_spin
            push!(pairs_same_spins, (n_site, m_site))
        end
        # Case mixed spins, can be same sites
        if m_spin != n_spin
            push!(pairs_mixed_spins, (n_site, m_site))
        end
    end
    # For the equal spins case the sites must be ordered, and we only need the unique cases, no repeated combinatiosn:
    pairs_same_spins = unique([(i, j) for (j, i) in pairs_same_spins if i < j])
    # For the mixed spins case we only select the unique combinations, with no restrictions on site ordering.#
    pairs_mixed_spins = unique(pairs_mixed_spins)
    return pairs_same_spins, pairs_mixed_spins
end
