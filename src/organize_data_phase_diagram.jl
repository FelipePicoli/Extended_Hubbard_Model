#=
    Creates a DataFrame object and stores it in a .csv file. 
=#
using DelimitedFiles
using DataFrames
using CSV

include("parser.jl")

parser = parse_commandline()

# Model parameters
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

measures = ["E_p", "E_p_bits", "S", "Coh", "E_GS", "m_cdw", "m_sdw"]
scalars = []

for U in U_values
    for V in V_values

        scalar_data = Dict{Symbol, Float64}()
        scalar_data[:U] = U
        scalar_data[:V] = V

        for measure in measures
            file_template = joinpath(results, "XXXXX_$(model)_L=$(L)_U=$(U)_V=$(V)_NPoints=$(Npoints).txt")
            filename = replace(file_template, "XXXXX" => measure)

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
info_output = joinpath(results, "scalars_$(model)_L=$(L)_NPoints=$(Npoints).csv")
CSV.write(info_output, df)

