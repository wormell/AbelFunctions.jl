type NeutralRecurrence{T<:Real,ffa,dffa}
  fa::ffa
  dfa::dffa
  α::T
  rad::T # radius of sector of analyticity
  Ψ::T # angle of sector of analyticity from real line in each direction
  dlogdfmax::T # max value of log-derivative of dfh/dx = d/dx (x*(1+fa(x))^α)) for x in sector of analyticity
  dfa0::T
  p::T
  sgn::Int
end
function NeutralRecurrence(fa,dfa,α,rad,Ψ,dlogdfmax;p=zero(typeof(α)),sgn=1)
  @assert sgn^2 == 1
  (Tα,Trad,TΨ,Tdlogdfmax,Tp) = promote(α,rad,Ψ,dlogdfmax,p)
  T = typeof(Tα)
  Tdfa0=dfa(zero(T))
  NeutralRecurrence{T,typeof(fa),typeof(dfa)}(fa,dfa,Tα,Trad,TΨ,Tdlogdfmax,Tdfa0,Tp,Int(sgn))
end

@compat (r::NeutralRecurrence)(x) = x*(1+r.fa(x^r.α))
map_trans(r::NeutralRecurrence,x) = x*(1+r.fa(x))^r.α


# type AbelFunction{T<:Real,ffa,dffa}
#   r::NeutralRecurrence{T,ffa,dffa}
#   cfs::Vector{T}
#   r0::T # radius of
