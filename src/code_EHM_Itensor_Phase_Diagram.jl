#= 
    DMRG for various U and V for the exteded Hubbard model.
=#
#=
    Implement argparser inputs 
=#
using ITensors
using ITensorMPS

using DelimitedFiles
using Random

using ArgParse

include("operators.jl")
include("measures.jl")
include("parser.jl")

parser = parse_commandline() 

# model parameters
L = parser["L"]
model = parser["model"]
results = parser["results"]
J = parser["J"]
U0 = parser["U0"]
Uf = parser["Uf"]
V0 = parser["V0"]
Vf = parser["Vf"]

Npoints = parser["Npoints"]

U_values = range(U0, stop=Uf, length=Npoints)
V_values = range(V0, stop=Vf, length=Npoints)

# dmrg parameters 
nsweeps = parser["nsweeps"]
sites = siteinds("Electron", L; conserve_qns=true)

maxdim = [50, 100, 200, 400, 800, 800] # , 1000]
cutoff = [1E-12]

Npart = floor(Int, L/2) 
Nup = Npart + L % 2 
Ndn = L - Nup 

# Correction to BOW phase which i'm not considering.
epsilon = 10e-2 

time_initial = time() 

for (i, U) in enumerate(U_values)
    for (j, V) in enumerate(V_values)
        #=
            Add check if files with Uf and Vf already exists. 
        =#
        println("U = $U, V = $V")
        # for running cases still not computed.
        H = H_EHM(L, J, U, V, sites)
        state = fill("Emp", L)

        if (V == 0 && U == 0)
          state = random_metallic_state(L, Nup, Ndn)
        elseif (V >= 2 * U + epsilon) 
          state = random_cdw_state(L, Nup, Ndn)
        elseif (V < 2 * U + epsilon) 
          state = random_sdw_state(L, Nup, Ndn)
        end

        psi0 = random_mps(sites, state; linkdims=10)

        # Start DMRG calculation:
        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        upd, dnd, updn = density_operators(L, psi)
        
        println("upd = $upd")
        println("dnd = $dnd")
        println("updn = $updn")

        #=
            Can reconstruct the densities from these three 
            results.
        =#
        charge_density = upd .+ dnd
        # removed the 1/2 factor
        magnetization = upd .- dnd
        doublon = updn

        #= 
            Implementation of the single-site entanglement 
            to be compared with the form of S = 1 - (1/L) \sum_{i} Tr (rho_i^2)
        =#
        S = average_single_site_entanglement(L, upd, dnd, doublon)
        E_p = average_entanglement_entropy(psi, L)

        filename = joinpath(results, "charge_density_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
        writedlm(filename, charge_density)

        filename = joinpath(results, "magnetization_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
        writedlm(filename, magnetization)
        
        filename = joinpath(results, "doublons_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
        writedlm(filename, doublon)

        filename = joinpath(results, "E_GS_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
        writedlm(filename, energy)

        filename = joinpath(results, "E_p_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
        writedlm(filename, E_p)

        filename = joinpath(results, "S_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
        writedlm(filename, S)
    end
end

filename = joinpath(results, "U_vals_NPoints=$(Npoints).txt")
writedlm(filename, U_values)
filename = joinpath(results, "V_vals_NPoints=$(Npoints).txt")
writedlm(filename, V_values)
