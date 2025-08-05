#= 
    All my simulations temporarily store results per simulation in .txt 
    files. When all iterations of a given simulation is finalized, auxiliary 
    scripts are run to store all these results from .txt files into .csv files. 
    
    This form of organization works for me so that i don't have to worry too much in 
    stoping simulations without having finished them. 

    These functions are used to store results from simulations and to organize them into 
    .csv files.
=#


#= 
    Used to store measures for the GS of the EHM at U, V points of the 
    phase diagram.
    This stores results for 
    - Energy 
    - Magnetization 
    - Charge density 
    - Doublons 

    - Single site entanglement
    
    - Order parameters 
        - M_SDW 
        - MCDW

    - 1RDM measures
        - Entanglment of particles
        - Coherence 1RDM 
        - Entanglement GAP
    - 2RDM measures 
        - Q_2 : quantum correlations
        - Coherence 2RDM 
        - Entanglement GAP
        - Cumulant matrix
=#
using Printf 
#=
    Stores measures for U, V point in .txt files.
=#
function store_EHM_GS_measures_results(results, model, L, U, V, charge_density, 
                                                magnetization,
                                                doublons,
                                                energy,
                                                E_p, 
                                                E_p_bits,
                                                single_site_entanglement,
                                                coherence_1rdm,
                                                gaps_1rdm, 
                                                xis_1rdm,
                                                coherence_2rdm,
                                                Q_2, 
                                                Q_2_bits,
                                                gaps_2rdm, 
                                                xis_2rdm,
                                                m_sdw,
                                                m_cdw, 
                                                Npoints)
    
    info_input = @sprintf("XXXXX_%s_L=%d_U=%.2f_V=%.2f_NPoints=%d.txt", model, L, U, V, Npoints)

    filename = joinpath(results, replace(info_input, "XXXXX" => "charge_density"))
    writedlm(filename, charge_density)

    filename = joinpath(results, replace(info_input, "XXXXX" => "magnetization"))
    writedlm(filename, magnetization)

    filename = joinpath(results, replace(info_input, "XXXXX" => "doublons"))
    writedlm(filename, doublons)

    filename = joinpath(results, replace(info_input, "XXXXX" => "E_GS"))
    writedlm(filename, energy)

    filename = joinpath(results, replace(info_input, "XXXXX" => "E_p"))
    writedlm(filename, E_p)

    filename = joinpath(results, replace(info_input, "XXXXX" => "E_p_bits"))
    writedlm(filename, E_p_bits)

    filename = joinpath(results, replace(info_input, "XXXXX" => "S"))
    writedlm(filename, single_site_entanglement)

    filename = joinpath(results, replace(info_input, "XXXXX" => "coh_1rdm"))
    writedlm(filename, coherence_1rdm)

    filename = joinpath(results, replace(info_input, "XXXXX" => "gaps_1rdm"))
    writedlm(filename, coherence_1rdm)

    filename = joinpath(results, replace(info_input, "XXXXX" => "xis_1rdm"))
    writedlm(filename, coherence_1rdm)

    filename = joinpath(results, replace(info_input, "XXXXX" => "coh_2rdm"))
    writedlm(filename, coherence_2rdm)

    filename = joinpath(results, replace(info_input, "XXXXX" => "Q_2"))
    writedlm(filename, Q_2)

    filename = joinpath(results, replace(info_input, "XXXXX" => "Q_2_bits"))
    writedlm(filename, Q_2_bits)

    filename = joinpath(results, replace(info_input, "XXXXX" => "gaps_2rdm"))
    writedlm(filename, coherence_1rdm)

    filename = joinpath(results, replace(info_input, "XXXXX" => "xis_2rdm"))
    writedlm(filename, coherence_1rdm)

    filename = joinpath(results, replace(info_input, "XXXXX" => "m_sdw"))
    writedlm(filename, m_sdw)

    filename = joinpath(results, replace(info_input, "XXXXX" => "m_cdw"))
    writedlm(filename, m_cdw)
end

#= 
    Organize a list of vectorial results a.k.a magnetization, occupation, ..., into 
    DataFrames to create .csv files
=#
function add_vectorial_measures_to_data_frames(vectors, vector_measures, info_input, results, L, U, V)

    row_data = Dict{Symbol, Any}()
    row_data[:U] = U
    row_data[:V] = V

    for measure in vector_measures

        filename = joinpath(results, replace(info_input, "XXXXX" => measure))               

                println("Searching for file: ", filename)
        # println(filename)

        if isfile(filename)
            vec = try
                readdlm(filename)
            catch e
                @warn "Failed to read $filename: $e"
                fill(NaN, L)
            end
        else
            @warn "File not found: $filename"
            vec = fill(NaN, L)
        end
        vec = string(join(vec, ", "))
        row_data[Symbol(measure)] = vec
    end
    push!(vectors, NamedTuple(row_data))
end


#= 
    Organize a list of scalar results to DataFrames for .csv files
=#
function add_scalar_measures_to_data_frames(scalars, measures, info_input, results, U, V)
    scalar_data = Dict{Symbol, Float64}()
    scalar_data[:U] = U
    scalar_data[:V] = V

    for measure in measures

        filename = joinpath(results, replace(info_input, "XXXXX" => measure))               

        if isfile(filename)
            val = try
                parse(Float64, strip(read(filename, String)))
            catch e
                @warn "Failed to parse $filename: $e"
                NaN
            end
        else
            @warn "File not found: $filename"
            val = NaN
        end
        scalar_data[Symbol(measure)] = val
    end
    push!(scalars, NamedTuple(scalar_data))
end

#= 
    Stores all results for range of U_values and V_values into .csv files that 
    can be read into DataFrames.
=#
function store_to_CSV_files(model, results, code, L, Npoints, U_values, V_values; is_only_pairs = false)

    scalar_measures = ["E_p", "E_p_bits", "S", "Q_2", "Q_2_bits", "gaps_2rdm", 
                        "xis_2rdm", "gaps_2rdm", "xis_2rdm",
                        "E_GS", "m_cdw", "m_sdw", "coh_1rdm", "coh_2rdm"]
    scalars = []

    vector_measures = ["magnetization", "charge_density", "doublons"]
    vectors = []
    
    for (i, U) in enumerate(U_values)
        for (j, V) in enumerate(V_values)
            info_input = @sprintf("XXXXX_%s_L=%d_U=%.2f_V=%.2f_NPoints=%d.txt", model, L, U, V, Npoints)

            add_scalar_measures_to_data_frames(scalars, scalar_measures, info_input, results, U, V)

            add_vectorial_measures_to_data_frames(vectors, vector_measures, info_input, results, L, U, V)
        end 
    end

    df = DataFrame(scalars)
    println(first(df, 30))
    info_output = joinpath(results, "scalars_$(model)_$(code)_L=$(L)_NPoints=$(Npoints).csv")
    CSV.write(info_output, df)

    df_vec = DataFrame(vectors)

    println(first(df_vec, 10))
    info_output_vec = joinpath(results, "vector_measures_$(model)_$(code)_L=$(L)_NPoints=$(Npoints).csv")
    CSV.write(info_output_vec, df_vec)
end
