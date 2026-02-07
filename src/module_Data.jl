using JLD2


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






















