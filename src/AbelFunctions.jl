module AbelFunctions
using Compat, ValidatedNumerics, IntervalRootFinding

include("general.jl")
include("neutral.jl")
include("newton.jl")
include("pomeaumanneville.jl")

end # module
