module AbelFunctions
using Compat, IntervalArithmetic, IntervalRootFinding, FastTransforms
#using IntervalRootFinding

include("general.jl")
include("taylortransform.jl")
include("neutral.jl")
include("neutral_newton.jl")
include("pomeaumanneville.jl")
include("abelfunction.jl")
include("abelfunction_newton.jl")

export AbelFunction, mapinv, mapinv_trans, map_trans, mapD, NeutralRecurrence

end # module
