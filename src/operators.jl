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
=#
function build_1_particle_rdm(psi) 
    L = length(psi)
    #=  
        N and not L since any site can have spin up or down
    =#
    rho_1 = zeros(ComplexF64, 2*L, 2*L) 

    Cupup = correlation_matrix(psi, "Cdagup", "Cup")
    Cdndn = correlation_matrix(psi, "Cdagdn", "Cdn")
    Cupdn = correlation_matrix(psi, "Cdagup", "Cdn")
    Cdnup = correlation_matrix(psi, "Cdagdn", "Cup")   

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
    rho_1 = rho_1 / L
end

#=
      os -= J, "Cdagup", i, "Cup", i + 1
      os -= J, "Cdagup", i + 1, "Cup", i
      os -= J, "Cdagdn", i, "Cdn", i + 1
=#
function build_2_particle_rdm(psi)
    
    os = OpSum()

    spins = ["up", "dn"]
    op = ["Cdag", "C"]

    for i in 1:L 
        for j in 1:L 
            for k in 1:L 
                for l in 1:L 
                end
            end
        end
    end

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
