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
    Builds the 1-particle reduced density matrix. 
    
    This matrix depends on just two correlators so it is 
    easily implemented by the correlation_matrix from ITensorMPS. 

=#
function build_1_particle_rdm(state, sites) 
    chain_size = length(sites)
    L = chain_size # length(state)
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
#=
    Create the basis for pairs of spins 
    1 = (1, up) , 2 = (1, dn), 3 = (2, up), ...
=#
function two_fermion_basis_pairs(L)
    num_modes = 2 * L
    pairs = []
    for m1 in 1:num_modes
        for m2 in (m1 + 1):num_modes
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
    
    # Basis for pairs of fermions.
    pairs = two_fermion_basis_pairs(L)

    dim = length(pairs)

    # @show dim

    rho_2 = zeros(ComplexF64, dim, dim)
    #= 
        Loops through all the configurations with no repeated indices and spins, 
        stores the elements rho_2[p, q] = <psi' | O | psi>. 
    =#
    for (p, (i, j)) in enumerate(pairs)
        # @show (p, (i, j))

        i_site, i_spin = mode_to_site_spin(i)
        j_site, j_spin = mode_to_site_spin(j)

        # @show (mode_to_site_spin(i), mode_to_site_spin(j))

        for (q, (k, l)) in enumerate(pairs)
            # @show (q, (k, l))
           
            os = OpSum() 

            k_site, k_spin = mode_to_site_spin(k)
            l_site, l_spin = mode_to_site_spin(l)

            # @show (mode_to_site_spin(k), mode_to_site_spin(l))

            os += "Cdag$i_spin", sites[i_site], "Cdag$j_spin", sites[j_sit], "C$l_spin", sites[l_sit], "C$k_spin", sites[k_site]

            O_mpo = MPO(os, sites)

            rho_2[p, q] = inner(phi', O_mpo, phi)
        end
    end
    # Force hermiticity
    rho_2 = (rho_2 + rho_2') / 2.0

    # Normalize by remaining number of particles 
    rho_2 = (2.0 / (L*(L-1))) * rho_2
    return rho_2
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

