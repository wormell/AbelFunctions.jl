module AbelFunctions
using Compat, ValidatedNumerics, IntervalRootFinding, FastTransforms

include("general.jl")
include("taylortransform.jl")
include("neutral.jl")
include("neutral_newton.jl")
include("pomeaumanneville.jl")
include("abelfunction.jl")
# include("abelfunction_newton.jl")

end # module
