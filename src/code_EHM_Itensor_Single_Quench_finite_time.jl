#= 
   Quench in U or V  
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
# model = parser["model"]
# results = parser["results"]
J = 1.0  # parser["J"]
U = -3.06122  # parser["U0"]
Uf = 8.0 # parser["Uf"]
V = -0.08163 # parser["V0"]
Vf = 8.0 # parser["Vf"]

# dmrg parameters 
nsweeps = parser["nsweeps"]
sites = siteinds("Electron", L; conserve_qns=true)

maxdim = [150, 100, 200, 400, 800, 800 , 1000, 1200]
cutoff = [1E-12]

Npart = floor(Int, L/2) 
Nup = Npart + L % 2 
Ndn = L - Nup 

# Correction to BOW phase which i'm not considering.
epsilon = 10e-2 

time_initial = time() 
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

psi0 = random_mps(sites, state; linkdims=4)

# Start DMRG calculation:
energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

upd, dnd, updn = density_operators(L, psi)

# println("upd = $upd")
# println("dnd = $dnd")
# println("updn = $updn")

#=
    Can reconstruct the densities from these three 
    results.
=#

upd, dnd, updn = density_operators(L, psi)

rho_1 = build_1_particle_rdm(psi)
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

# S = average_single_site_entanglement(L, upd, dnd, doublon)
S, S_bit = Ep(rho_1, L)
corr = quantum_coherence(rho_1, L)

println("<ψ | ψ > = ", inner(psi,psi)) # scalar(dag(psi)*psi))

# println("S = ", S)
println("E_p_bit = ", S_bit - log2(L))
println("E_p = ", S - log(L))

#= 
    Implementation of the single-site entanglement 
    to be compared with the form of S = 1 - (1/L) \sum_{i} Tr (rho_i^2)
=#
#=
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

filename = joinpath(results, "U_vals_NPoints=$(Npoints).txt")
writedlm(filename, U_values)
filename = joinpath(results, "V_vals_NPoints=$(Npoints).txt")
writedlm(filename, V_values)
=#
