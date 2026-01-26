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
#=
    Generates the MPO for the EHM Hamiltonian
    with strengths J, U and V.
    Requires a SiteType sites.
=#
function H_EHM(N, J, U, V, sites)
    os = OpSum()
    for i in 1:(N - 1)
      # Knetic
      os -= J, "Cdagup", i, "Cup", i + 1
      os -= J, "Cdagup", i + 1, "Cup", i
      os -= J, "Cdagdn", i, "Cdn", i + 1
      os -= J, "Cdagdn", i + 1, "Cdn", i
      # Nearest-neighbours
      os += V, "Ntot", i, "Ntot", i + 1
    end
    # on-site
    for i in 1:N
      os += U, "Nupdn", i
    end
    return MPO(os, sites)
end
#=
    Returns the density of up and down
    electrons.
=#
function fermionic_density_operators(N, psi, sites)
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
function build_1_particle_rdm(state, sites)
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



#=
    Create the basis for pairs of spins
    1 = (1, up) , 2 = (1, dn), 3 = (2, up), ...
=#
function two_fermion_basis_pairs(modes)
    pairs = []
    for m1 in 1:length(modes)
        for m2 in (m1 + 1):length(modes)
            push!(pairs, (m1, m2))
        end
    end
    return pairs
end

#=
    Gets the indexed site and spin.
=#
function mode_to_site_spin(m::Int)
    site = (m + 1) ÷ 2
    spin = isodd(m) ? "up" : "dn"
    return site, spin
end
# Create list of all modes (site, spin) for selected sites
function site_spin_modes(sites)
    modes = []
    for s in sites
        push!(modes, (s, "up"))
        push!(modes, (s, "dn"))
    end
    return modes
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
function build_2_particle_rdm(phi, sites)

    L = length(sites) # length(phi)

    @show typeof(sites)

    modes = site_spin_modes(sites)
    pairs = two_fermion_basis_pairs(modes)

    # Basis for pairs of fermions.
    # pairs = two_fermion_basis_pairs(L)

    dim = length(pairs)

    # @show dim

    rho_2 = zeros(ComplexF64, dim, dim)

    # Ugly gambiarra necessary to obtain the site number.
    site_number(site_index) = parse(Int, match(r"n=(\d+)", string(tags(site_index))).captures[1])
    #=
        Loops through all the configurations with no repeated indices and spins,
        stores the elements rho_2[p, q] = <psi' | O | psi>.
    =#
    for (p, (i, j)) in enumerate(pairs)
        # @show (p, (i, j))

        i_site, i_spin = modes[i]
        j_site, j_spin = modes[j]

        i_site = site_number(i_site)
        j_site = site_number(j_site)
        # @show i_site, j_site

        for (q, (k, l)) in enumerate(pairs)
            # @show (q, (k, l))

            os = OpSum()

            k_site, k_spin = modes[k]
            l_site, l_spin = modes[l]

            k_site = site_number(k_site)
            l_site = site_number(l_site)
            # @show k_site l_site

            os += "Cdag$i_spin", i_site, "Cdag$j_spin", j_site, "C$l_spin", l_site, "C$k_spin", k_site

            O_mpo = MPO(os, sites)

            rho_2[p, q] = inner(phi', O_mpo, phi)
        end
    end
    # Force hermiticity
    rho_2 = (rho_2 + rho_2') / 2.0
    # Normalize by remaining number of particles
    rho_2 = (2.0 / (L*(L-1))) * rho_2
    @show rho_2
    return rho_2
end


# Still need fixing.
function sites_to_modes(sites)
    return [(s, spin), for s in sites for spin in ('up', 'dn')]
end

function select_unique_pairs()
end
function get_all_modes_pairs()
end
function handler_disjoint_spins()
end
function handler_joint_spins()
end
function two_fermions_basis_sites_pairs_spin_sectors
end
function compute_block_2rdm()
end
function gbuild_2_particle_rdm()
end
