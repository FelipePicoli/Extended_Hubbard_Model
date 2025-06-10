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
function build_1_particle_rdm(state) 
    L = length(state)
    #=  
        N and not L since any site can have spin up or down
    =#
    rho_1 = zeros(ComplexF64, 2*L, 2*L) 

    Cupup = correlation_matrix(state, "Cdagup", "Cup")
    Cdndn = correlation_matrix(state, "Cdagdn", "Cdn")
    Cupdn = correlation_matrix(state, "Cdagup", "Cdn")
    Cdnup = correlation_matrix(state, "Cdagdn", "Cup")   

    for i in 1:L
        for j in 1:L
            # Blocks of the correlation matrix.
            i_up = 2 * (i-1) + 1 
            i_dn = 2 * (i-1) + 2
            j_up = 2 * (j-1) + 1
            j_dn = 2 * (j-1) + 2

            rho_1[i_up, j_up] = Cupup[i,j]
            rho_1[i_up, j_dn] = Cupdn[i,j]
            rho_1[i_dn, j_up] = Cdnup[i,j]
            rho_1[i_dn, j_dn] = Cdndn[i,j]
        end 
    end
    return rho_1 = rho_1 / L
end
=#


function build_1_particle_rdm(state, up, dn, bulk_range=nothing)
    
    L = length(state) 

    if isnothing(bulk_range)
        bulk_range = (L ÷ 4 + 1):(3*L ÷ 4)
    end

    Cupup = correlation_matrix(state, "Cdagup", "Cup")
    Cdndn = correlation_matrix(state, "Cdagdn", "Cdn")
    Cupdn = correlation_matrix(state, "Cdagup", "Cdn")
    Cdnup = correlation_matrix(state, "Cdagdn", "Cup")   

    bulk_sites = collect(bulk_range)
    L_b = length(bulk_sites)

    rho_bulk = zeros(ComplexF64, 2*L_b, 2*L_b)

    for (bi, i) in enumerate(bulk_sites)
        for (bj, j) in enumerate(bulk_sites)
            i_up = 2 * (bi - 1) + 1
            i_dn = 2 * (bi - 1) + 2
            j_up = 2 * (bj - 1) + 1
            j_dn = 2 * (bj - 1) + 2

            rho_bulk[i_up, j_up] = Cupup[i, j]
            rho_bulk[i_dn, j_dn] = Cdndn[i, j]
            rho_bulk[i_up, j_dn] = Cupdn[i, j]
            rho_bulk[i_dn, j_up] = Cdnup[i, j]
        end
    end

    # Total number of particles in the bulk
    N_bulk = sum(up[bulk_sites]) + sum(dn[bulk_sites])

    # Normalize so Tr[rho_bulk] = 1
    rho_bulk ./= N_bulk

    return rho_bulk, N_bulk
end

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
function mode_to_site_spin(m::Int)
    site = (m + 1) ÷ 2
    spin = isodd(m) ? "up" : "dn"
    return site, spin
end

function build_2_particle_rdm(phi, sites)

    L = length(phi)

    pairs = two_fermion_basis_pairs(L)

    dim = length(pairs)

    @show dim

    rho_2 = zeros(ComplexF64, dim, dim)

    for (p, (i, j)) in enumerate(pairs)
        @show (p, (i, j))

        i_site, i_spin = mode_to_site_spin(i)
        j_site, j_spin = mode_to_site_spin(j)

        # @show (mode_to_site_spin(i), mode_to_site_spin(j))

        for (q, (k, l)) in enumerate(pairs)

            os = OpSum() 

            @show (q, (k, l))

            k_site, k_spin = mode_to_site_spin(k)
            l_site, l_spin = mode_to_site_spin(l)

            # @show (mode_to_site_spin(k), mode_to_site_spin(l))

            os -= "Cdag"*i_spin, i_site, "Cdag"*j_spin, j_site, "C"*k_spin, k_site, "C"*l_spin, l_site

            O_mpo = MPO(os, sites)

            rho_2[p, q] = inner(phi', O_mpo, phi)
        end
    end

    is_hermitian = true
    tolerances = [1e-12, 1e-13, 1e-14, 1e-16, 1e-18, 1e-22]
    for tolerance in tolerances
        for k in 1:dim
            for l in k:dim
                if abs(rho_2[k,l] - conj(rho_2[l,k])) > tolerance
                    println("Non-Hermitian element found at ($k, $l)")
                    is_hermitian = false
                end
            end
        end
        if is_hermitian
            println("rho_2 Hermitian for tolerance $tolerance.")
        else
            println("rho_2 NOT Hermitian for tolerance $tolerance.")
        end
    end
    # @pt rho_2
    @show tr(rho_2)
    rho_2 = (2 / (L*(L-1))) * rho_2
    @show tr(rho_2)
    return rho_2
end

#= 
    A collection of functions to generate 
    to generate random initial product states  
    using the random_mps method. 
    
    These states corresponds to diferent 
    phases of matter and are used according 
    to the phase diagram of the EHM.
=#

function random_metallic_state(L, Nup, Ndn)
    state = fill("Emp", L)
    p = Nup + Ndn
    for i in 1:L
        j = L - i
        if(p > j)
            state[j] = "UpDn"
            p -= 2
        elseif (p > 0) 
            state[j] = j % 2 == 1 ? "Up" : "Dn"
            p -= 1
        end
    end
    return state
end

function random_cdw_state(L, Nup, Ndn)
    state = fill("Emp", L)
    Nup_extra = Int.(Nup % Ndn)
    for i in 1:2:(L - Nup_extra)
        state[i] = "UpDn"
        Nup-=1
    end
    if Nup != 0
        state[L] = "Up"
    end
    return state
end 
function random_sdw_state(L, Nup, Ndn)
    state = fill("Emp", L)
    Nup_extra = Int.(Nup % Ndn)
    for i=1:2:(L-Nup_extra)
        state[i] = "Up"
        state[i+1] = "Dn"
    end
    if Nup != 0 
        state[L] = "Up"
    end
    return state 
end


function random_ps_state(L, Nup, Ndn)
    state = fill("Emp", L)
    
    Ndbl = min(Nup, Ndn)
    Nup -= Ndbl
    Ndn -= Ndbl
    
    for i in 1:Ndbl
        state[i] = "UpDn"
    end
    
    next_site = Ndbl + 1
    
    if Nup > 0
        state[next_site] = "Up"
        next_site += 1
    elseif Ndn > 0
        state[next_site] = "Dn"
        next_site += 1
    end
    return state
end
