module AbelFunctions
using Compat, IntervalArithmetic, FastTransforms
#using IntervalRootFinding

import Base.^

^(a::Complex{Interval{T}}, b::Interval{T}) where T = exp(b*log(a))

include("general.jl")
include("taylortransform.jl")
include("neutral.jl")
include("neutral_newton.jl")
include("pomeaumanneville.jl")
include("abelfunction.jl")
include("abelfunction_newton.jl")

export AbelFunction, mapinv, mapinv_trans, map_trans, mapD, NeutralRecurrence

end # module
