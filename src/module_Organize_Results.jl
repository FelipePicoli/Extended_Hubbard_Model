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
function store_EHM_GS_measures_results(results, model, L, U, V, Npoints, dict_results)
    info_template = @sprintf("XXXXX_%s_L=%d_U=%.2f_V=%.2f_NPoints=%d.txt", model, L, U, V, Npoints)
    for (key, value) in dict_results
        filename = joinpath(results, replace(info_template, "XXXXX" => key))
        writedlm(filename, value)
    end
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

    scalar_measures = [
                       "E_p", "E_p_bits","coh_1rdm",
                       "Q_2", "Q_2_bits", "coh_2rdm",
                        "average_single_site_entanglement", "energy", "m_cdw", "m_sdw", 
                        ]
    scalars = []

    vector_measures = ["magnetization", "charge_density", "doublons", "Omega_1rdm", "Omega_2rdm"]
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
