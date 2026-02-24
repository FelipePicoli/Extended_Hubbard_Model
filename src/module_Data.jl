using JLD2
using DataFrames
using CSV


function data_structure(L)
   return Dict(
        "scalar" =>         Dict(
                                "energy_GS" => 0.0,
                                "von_neumann_1rdm_GS" => 0.0,
                                "quantum_coherence_1rdm_GS" => 0.0,
                                "quantum_coherence_2rdm_GS" => 0.0,
                            ),
        "vector" =>         Dict(
                                "single_site_entanglement_GS" => zeros(L),
                                "charge_density_GS" => zeros(L),
                                "magnetization_GS" => zeros(L),
                                "doublons_GS" => zeros(L),

                                "entanglement_spectrum_1rdm_GS" => zeros(2*L),
                                "entanglement_spectrum_2rdm_GS" => zeros(L*(2*L-1)),
                            )
        )
end


function store_results(results, info_output, path_results)
    if !isdir(path_results)
        mkpath(path_results)
    end
    full_path = joinpath(path_results, info_output)

    println("full_path results = ", full_path)

    jldsave(full_path; results)

    println("Results saved to: $full_path")
end

fmt(x) = @sprintf("%.2f", x)

"""
    Transforms a collection of .jld2 files generated point-by-point on a phase-diagram mesh
    into consolidated CSV files.

    Expected JLD2 structure:
        Dict(
            "scalar" => Dict(String => Number),
            "vector" => Dict(String => AbstractArray)
        )

    Creates:
        *_scalar.csv
        *_vectors.csv
"""
function organize_data_phase_diagram(L, U0, Uf, V0, Vf, N_points, path; model="EHM_ITensor", output_name="results")

    # Parameter grids
    U_values = range(U0, Uf; length=N_points)
    V_values = range(V0, Vf; length=N_points)

    scalar_rows = Dict{String,Any}[]

    vector_rows = DataFrame(
        U=Float64[],
        V=Float64[],
        observable=String[],
        site=Int[],
        value=Float64[])

    for U in U_values, V in V_values

        path_jld2 = joinpath(path, "jld2_files")
        fname = joinpath(path_jld2,"$(output_name)_$(model)_L=$(L)_U=$(fmt(U))_V=$(fmt(V))_Npoints=$(N_points).jld2")

        if !isfile(fname)
            println("File not found: $fname. Exiting.")
            return
        end

        # Load data
        data = load(fname, "results")

        # ---------- scalar ----------
        scalar_dict = data["scalar"]
        scalar_row = Dict{String,Any}(
            "U" => U,
            "V" => V
        )

        for (k, v) in scalar_dict
            scalar_row[k] = v
        end

        push!(scalar_rows, scalar_row)

        # ---------- vector ----------
        for (obs, vec) in data["vector"]
            vec = vec[:]   # flatten
            for (site, value) in enumerate(vec)
                push!(vector_rows, (
                    U=U,
                    V=V,
                    observable=String(obs),
                    site=site,
                    value=value
                ))
            end
        end
    end

    # Output filenames
    output_file_name =
        "$(path)/$(model)_L=$(L)" *
        "_U0=$(fmt(U0))_Uf=$(fmt(Uf))" *
        "_V0=$(fmt(V0))_Vf=$(fmt(Vf))" *
        "_NPoints=$(N_points)"

    CSV.write("$(output_file_name)_scalar.csv", DataFrame(scalar_rows))
    CSV.write("$(output_file_name)_scalar.csv", scalar_rows)
    CSV.write("$(output_file_name)_vectors.csv", vector_rows)

    return nothing
end
# organize_data_phase_diagram(5, -6.0, 6.0, -6.0, 6.0, 100, "../results/")

