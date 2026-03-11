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
include("module_Fermionic_Operators.jl")
include("module_Quantum_Info_Tools.jl")

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
            rho_1[i, j + L] = rho_updn[i, j]
            rho_1[i + L, j] = rho_updn[i, j]
            rho_1[i + L, j + L] = rho_dndn[i,j]
        end
    end
    return rho_1 = rho_1 / L
end

let
    N = 4

    Npart = floor(Int, N/2)
    Nup = Npart + N % 2
    Ndn = N - Nup

    sites = siteinds("Electron", N; conserve_qns = true)
    J = 1.0
    # Cuts along the extended Hubbard model phase diagram.
    U_vals = [0.0] # 1.6] # , -0.95, -1.76, -2.85, -7.73]
    V_vals = [0.0] # range(0, -1.7, length=2)
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

    nsweeps = 12
    cutoff = 1E-8
    maxdim = [2,2,2,2,2,8,8,8,8,8,20,20,20,20,20,50,50,50,50,50,100,100,200,200,400]
    state = [isodd(n) ? "Up" : "Dn" for n in 1:N]
    psi0 = random_mps(sites, state; linkdims = 10)

    for U in U_vals
        E_p1 = zeros(length(V_vals))
        E_p2 = zeros(length(V_vals))
        central_site_entanglement = zeros(length(V_vals))

        for (i, V) in enumerate(V_vals)
            @show i, U,  V

            H = J*H_J + U*H_U + V*H_V
            energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

            up, dn, updn = density_operators(N, psi, sites)
            central_site_entanglement[i] = single_site_entanglement(N, floor(Int, N/2), up, dn, updn)

            @show central_site_entanglement[i]

            cutoff_eigenvals = 1e-16

            rho_1 = get_1_particle_rdm(psi, sites)
            # diagonalize rho
            eigenvalues = eigen(Hermitian(rho_1)).values
            eigenvalues = max.(0.0, eigenvalues)
            eigenvalues = eigenvalues[eigenvalues .> cutoff_eigenvals]
            shifted_von_neumann_1 = -sum(eigenvalues.*log2.(eigenvalues)) - log2(N)
            E_p1[i] = shifted_von_neumann_1

            rho_2 = get_2_particle_rdm(psi, sites)
            # diagonalize rho
            eigenvalues = eigen(Hermitian(rho_2)).values
            eigenvalues = max.(0.0, eigenvalues)
            eigenvalues = eigenvalues[eigenvalues .> cutoff_eigenvals]
            shifted_von_neumann_2 = -sum(eigenvalues.*log2.(eigenvalues)) - log2(N * (N-1.0) / 2.0)
            E_p2[i] = shifted_von_neumann_2

            @show E_p1[i]
            @show E_p2[i]
        end
        df = DataFrame(
                       V = [V for V in V_vals],
                       E_p1 = [E_p1[i] for i=1:length(V_vals)],
                       E_p2 = [E_p2[i] for i=1:length(V_vals)],
                       single_site_entanglemenet = [central_site_entanglement[i] for i=1:length(V_vals)])
        CSV.write("result_U=$(round(U, digits=2))_N=$(N).csv", df)
    end
end
