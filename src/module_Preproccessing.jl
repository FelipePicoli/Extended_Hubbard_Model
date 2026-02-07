include("module_States_Ansatze.jl")

function preprocessing_simulation(L, J, U, V, Nup, Ndn, m, paths)
    if isfile(paths.sites)
        loaded_sites = JLD2.load(paths.sites)["sites"]
        println(length(loaded_sites))
        if length(loaded_sites) == L
            println("-- Loading stored siteinds for L = $(L)")
            sites = loaded_sites
        else
            # Case it exists but for different L value.
            println("-- Building siteinds")
            sites = siteinds("Electron", L; conserve_qns=true)
            jldsave(paths.sites; sites)
        end
    else
        println("-- Building siteinds.")
        sites = siteinds("Electron", L; conserve_qns=true)
        println("-- Storing siteinds.")
        jldsave(paths.site; sites)
    end

    println("Obtaining random MPS")
    if isfile(paths.random_mps)
        println("-- Using previous point random MPS.")
        psi0 = JLD2.load(paths.random_mps)["psi"]
        @show psi0
    else
        println("Building random MPS.")
        state = state_ehm_diagram(L, Nup, Ndn, U, V)
        psi0 = random_mps(sites, state; linkdims=m)
        @show psi0
    end
    println("Building H")
    if isfile(paths.hamiltonian_mpos)
        println("-- Loading Hamiltonian MPOS.")
        hamiltonian_mpos = JLD2.load(paths.hamiltonian_mpos)["dict"]
        H_J, H_U, H_V = hamiltonian_mpos["H_J"], hamiltonian_mpos["H_U"], hamiltonian_mpos["H_V"]
        @show hamiltonian_mpos
        H_hamiltonian = J * H_J + U*H_U + V*H_V
    else
        println("-- Building Hamiltonian MPOS.")
        H_J, H_U, H_V = get_EHM_Hamiltonian_components(L, sites)
        H_hamiltonian = J * H_J + U*H_U + V*H_V
        dict = Dict("H_J" => H_J, "H_U" =>  H_U, "H_V" => H_V)
        jldsave(stored_hamiltonian_mpos; dict)
    end
    truncate!(H_hamiltonian; cutoff=1e-15)

    return H_hamiltonian, psi0, sites
end
