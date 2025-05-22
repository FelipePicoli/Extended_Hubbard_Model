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
    A collection of functions to generate 
    to generate random initial product states  
    using the random_mps method. 
    
    These states corresponds to diferent 
    phases of matter and are used according 
    to the phase diagram of the EHM.  
=#
function random_metallic_state(L, Nup, Ndn)
    state = fill("Emp", L)
    # Random order of site indices
    up_sites = shuffle(1:L)[1:Nup]
    for i in up_sites
        state[i] = "Up"
    end
    dn_sites = shuffle(1:L)
    added_dn = 0
    for i in dn_sites
        if added_dn == Ndn
            break
        end
        if state[i] == "Emp"
            state[i] = "Dn"
        elseif state[i] == "Up"
            state[i] = "UpDn"
        end
        added_dn += 1
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

