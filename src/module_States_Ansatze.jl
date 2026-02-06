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


#=
    Returns a state for the point (U, V) of the Extended Hubbard Model.
=#
function state_ehm_diagram(L, Nup, Ndn, U, V)

    state = fill("Emp", L)

    region = get_region(U, V)

    # Weak coupling = metallic
    if region == "METALLIC"
        state = random_metallic_state(L, Nup, Ndn)
    # CDW
    elseif region == "CDW"
        state = random_cdw_state(L, Nup, Ndn)
    elseif region == "SDW"
        state = random_sdw_state(L, Nup, Ndn)
    else
        state = random_ps_state(L, Nup, Ndn)
    end
    return state
end
#=
    Trying to separate the phase-diagram states not only
    in big squared blocks.
    See notebook ploting_ehm_diagram_selection_of_states.ipynb for an example of plot
    using this function. I try to made it look like the plot schematic
    diagram for the EHM.
=#
function get_region(u, v)
    alpha = 0.5

    # Boundary functions
    get_metallic_lower_boundary(u) = u <= 0 ? -exp(alpha * u) : -alpha * u - 1.0
    get_metallic_upper_boundary_h(u) = -0.1 * u
    get_metallic_upper_boundary_negative(u) = -2.5 * alpha * u

    lower = get_metallic_lower_boundary(u)
    upper_h = get_metallic_upper_boundary_h(u)
    upper_neg = get_metallic_upper_boundary_negative(u)

    if v <= lower
        return "PS"
    elseif (u < 0 && v < 0 && v > lower) ||
           (u > 0 && v > 0 && v < upper_h) ||
           (u > 0 && v < 0 && v > lower && v < upper_neg)
        return "METALLIC"
    elseif (u < 0 && v > 0) ||
           (u > 0 && v > 0 && v >= 2 * u)
        return "CDW"
    elseif (u > 0 && v < 0 && v > lower) ||
           (u > 0 && v > 0 && v >= upper_h && v < 2 * u)
        return "SDW"
    end
    return "PS"
end
