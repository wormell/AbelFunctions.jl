type NeutralRecurrence{T<:Real,ffa,dffa}
  fa::ffa
  dfa::dffa
  α::T
  rad::T # radius of sector of analyticity
  Ψ::T # angle of sector of analyticity
  dlogdfmax::T # max value of log-derivative of dfh/dx = d/dx (x*(1+fa(x))^α)) for x in sector of analyticity
  dfa0::T
  p::T
  sgn::Int
end
function NeutralRecurrence(fa,dfa,α,rad,Ψ,dlogdfmax;p=zero(typeof(α)),sgn=1,
    dfa0=dfa(zero(promote_type(p,α,rad,Ψ,dlogdfmax))))
    @assert sgn^2 == 1
  (Tα,Trad,TΨ,Tdlogdfmax,Tdfa0,Tp) = promote(α,rad,Ψ,dlogdfmax,dfa0,p)
  NeutralRecurrence{T,typeof(fa),typeof(dfa)}(fa,dfa,Tα,Trad,TΨ,Tdlogdfmax,Tdfa0,Tp,Int(sgn))
end

# type AbelFunction{T<:Real,ffa,dffa}
#   r::NeutralRecurrence{T,ffa,dffa}
#   cfs::Vector{T}
#   r0::T # radius of
