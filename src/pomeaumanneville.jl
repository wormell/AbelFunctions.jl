@compat struct PomFA{T<:Real}
  α::T
  pow2α::T
end
PomFA(α) = PomFA(α,exp2(α))
@compat (p::PomFA)(x) = p.pow2α*x

@compat struct PomDFA{T<:Real}
  α::T
  pow2α::T
end
PomDFA(α) = PomDFA(α,exp2(α))
@compat (p::PomDFA)(x) = p.pow2α

function Pom(α,T=typeof(float(α)))
  αT = convert(T,α)
  rad = one(T)/3 # may need to be smaller for Newton convergence
  NeutralRecurrence(PomFA(αT),PomDFA(αT),αT,rad,T(π)/2,exp2(αT+1)*αT;p=zero(T),sgn=1)
end
