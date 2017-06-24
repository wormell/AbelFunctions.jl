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
  pow2αT = exp2(αT)
  Ψ = T(π)/4max(1,α)
  rad = exp2(-α-1)# for convergence of derivatives of w, log(f̂)-1 
  NeutralRecurrence(PomFA(αT),PomDFA(αT),αT,rad,T(π)/2,2pow2αT*αT,pow2αT,(1+2rad*pow2αT)^(-αT)/sqrt(T(2));p=zero(T),sgn=1)
end
