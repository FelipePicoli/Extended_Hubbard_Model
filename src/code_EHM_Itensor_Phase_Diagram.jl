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

maxdim = [50, 100, 200, 400, 800, 800]
cutoff = [1E-12]

Npart = floor(Int, L/2) 
Nup = Npart + L % 2 
Ndn = L - Nup 

# Correction to BOW phase which i'm not considering.
epsilon = 10e-2 

for (i, Uf) in enumerate(U_values)
    for (j, Vf) in enumerate(V_values)
    
        #=
            Add check if files with Uf and Vf already exists. 
        =#
        println("Uf = $Uf, Vf = $Vf")
#        if((Uf >= 6.0 && Uf <= Umax) && (Vf <= 2 * Uf))
        H = H_EHM(L, J, Uf, Vf, sites)
        state = fill("Emp", L)
        if (Vf == 0 && Uf == 0)
          state = random_metallic_state(L, Nup, Ndn)
        elseif (Vf >= 2 * Uf + epsilon) 
          state = random_cdw_state(L, Nup, Ndn)
        elseif (Vf < 2 * Uf + epsilon) 
          state = random_sdw_state(L, Nup, Ndn)
        end

        psi0 = random_mps(sites, state; linkdims=10)

        # Start DMRG calculation:
        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

        upd, dnd = density_operators(L, psi)

        charge_density = upd .+ dnd
        magnetization = (upd .- dnd) / 2

        filename = joinpath(results, "charge_density_$(model)_L=$(L)_U=$(Uf)_V=$(Vf)_NPoints=$(Npoints).txt")
        writedlm(filename, charge_density)

        filename = joinpath(results, "magnetization_$(model)_L=$(L)_U=$(Uf)_V=$(Vf)_NPoints=$(Npoints).txt")
        writedlm(filename, magnetization)

        filename = joinpath(results, "E_GS_$(model)_L=$(L)_U=$(Uf)_V=$(Vf)_NPoints=$(Npoints).txt")
        writedlm(filename, energy)
#        end
    end
end
