#=
    Creates a DataFrame object and stores it in a .csv file. 
=#
using DelimitedFiles
using DataFrames
using CSV

using Printf

include("parser.jl")

let

    parser = parse_commandline()

    # Model parameters
    L = parser["L"]
    model = parser["model"]
    code = parser["code"]

    results = parser["results"]
    J = parser["J"]

    U0 = parser["U0"]
    Uf = parser["Uf"]
    V0 = parser["V0"]
    Vf = parser["Vf"]

    Npoints = parser["Npoints"]

    U_values = range(U0, stop=Uf, length=Npoints)
    V_values = range(V0, stop=Vf, length=Npoints)


    measures = ["E_p", "E_p_bits", "S", "coherence", "E_GS", "m_cdw", "m_sdw"]
    scalars = []

    for U in U_values
        for V in V_values
            
            scalar_data = Dict{Symbol, Float64}()
            scalar_data[:U] = U
            scalar_data[:V] = V

            for measure in measures

                info_input = @sprintf("XXXXX_%s_L=%d_U=%.2f_V=%.2f_NPoints=%d.txt", model, L, U, V, Npoints)
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
    end

    df = DataFrame(scalars)
    println(first(df, 30))
    info_output = joinpath(results, "scalars_$(model)_$(code)_L=$(L)_NPoints=$(Npoints).csv")
    CSV.write(info_output, df)

    vec_measures = ["magnetization", "charge_density"]
    vector_data = []

    for U in U_values
        for V in V_values

            row_data = Dict{Symbol, Any}()
            row_data[:U] = U
            row_data[:V] = V

            for measure in vec_measures
                info_input = @sprintf("XXXXX_%s_L=%d_U=%.2f_V=%.2f_NPoints=%d.txt", model, L, U, V, Npoints)
                filename = joinpath(results, replace(info_input, "XXXXX" => measure))               

                println("Searching for file: ", filename) # ADD THIS LINE FOR DEBUGGING
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
            push!(vector_data, NamedTuple(row_data))
        end
    end

    df_vec = DataFrame(vector_data)

    println(first(df_vec, 10))
    # Optionally: save as CSV (note that storing arrays directly in CSV may not be very usable)
    info_output_vec = joinpath(results, "vector_measures_$(model)_$(code)_L=$(L)_NPoints=$(Npoints).csv")
    CSV.write(info_output_vec, df_vec)
end
