#=
    Reproduces results of Fig. 5 of PhysRevB.92.075423.
    The list of E_p values are stored in .csv files in the
    folder where the code are run named
    "result_U=$(U_Value)...csv".
=#


using ITensors, ITensorMPS
using CSV
using DataFrames
using LinearAlgebra

#=
    Includes functions that generate initial
    states for each reagion of the phase diagram  used
    as initial mps ansatz.
    The states follows the ones shown in (3), (5) and (6) from
    the aforementioned paper.
=#
include("module_States_Ansatze.jl")
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

let

    N = 20

    Npart = floor(Int, N/2)
    Nup = Npart + N % 2
    Ndn = N - Nup

    sites = siteinds("Electron", N; conserve_qns = true)

    J = 1.0

    # Cuts along the extended Hubbard model phase diagram.
    U_vals = [-1.0]
    V_vals =[1.0] #range(-3, 3, length=50)


    #=
        Generates OpSum objects for each block of
        the EHM Hamiltonian.
    =#
    os_J = OpSum()
    for i in 1:(N-1)
        os_J -= 1.0, "Cdagup", i, "Cup", i + 1
        os_J -= 1.0, "Cdagup", i + 1, "Cup", i
        os_J -= 1.0, "Cdagdn", i, "Cdn", i + 1
        os_J -= 1.0, "Cdagdn", i + 1, "Cdn", i
    end
    os_U = OpSum()
    for i in 1:N
        os_U += 1.0, "Nupdn", i
    end
    os_V = OpSum()
    for i in 1:(N - 1)
        os_V += 1.0, "Ntot", i, "Ntot", i + 1
    end

    H_J, H_U, H_V = MPO(os_J ,sites), MPO(os_U ,sites), MPO(os_V ,sites)

    nsweeps = 25
    cutoff = 1E-8
    maxdim = [2,2,2,2,2,8,8,8,8,8,20,20,20,20,20,50,50,50,50,50,100,100,200,200,400]
    # bound dimension
    m = 10

    for U in U_vals
        @show U
        E_p = zeros(length(V_vals))
        for (i, V) in enumerate(V_vals)

            state = state_ehm_diagram(N, Nup, Ndn, U, V)
            psi0 = random_mps(sites, state; linkdims=m)

            H = J*H_J + U*H_U + V*H_V
            energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

            rho_1 = get_1_particle_rdm(psi, sites)
            # diagonalize rho
            eigenvalues = eigen(Hermitian(rho_1)).values
            eigenvalues = max.(0.0, eigenvalues)
            shifted_von_neumann = -sum(eigenvalues.*log.(eigenvalues)) - log(N)
            E_p[i] = shifted_von_neumann
        end
        df = DataFrame(
                       V = [V for V in V_vals],
                       E_p = [E_p[i] for i=1:length(V_vals)])
        CSV.write("result_U=$(round(U, digits=2)).csv", df)
    end
end
