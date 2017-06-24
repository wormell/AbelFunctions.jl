@compat struct NeutralRecurrence{T<:Real,ffa,dffa}
  fa::ffa
  dfa::dffa
  α::T
  rad::T # radius of sector of analyticity
  Ψ::T # angle of sector of analyticity from real line in each direction
  dlogdfmax::T # max value of log-derivative of df̂/dẑ = d/dẑ (ẑ*(1+fa(ẑ))^α)) for ẑ in sector of analyticity
  F::T # maxabs of derivative of fa
  M::T # inf Re(d/dẑ (1+fa(ẑ))^-α)/dfa0 which should be between 0 and 1
  dfa0::T
  p::T
  sgn::Int
end


function NeutralRecurrence(fa,dfa,α,rad,Ψ,dlogdfmax,F,M;p=zero(typeof(α)),sgn=1)
  @assert sgn^2 == 1
  @assert α>0
  (Tα,Trad,TΨ,Tdlogdfmax,TF,TM,Tp) = promote(α,rad,Ψ,dlogdfmax,F,M,p)
  T = typeof(Tα)
  Tdfa0=dfa(zero(T))
  NeutralRecurrence{T,typeof(fa),typeof(dfa)}(fa,dfa,Tα,Trad,TΨ,Tdlogdfmax,TF,TM,Tdfa0,Tp,Int(sgn))
end

hat(z,r) = ((z-r.p)*r.sgn)^r.α
unhat(ẑ,r) = r.sgn*ẑ^(1/r.α)+r.p

@compat (r::NeutralRecurrence)(x) = (x-r.p)*(1+r.fa(hat(x,r)))
map_trans(r::NeutralRecurrence,x) = x*(1+r.fa(x))^r.α



# type AbelFunction{T<:Real,ffa,dffa}
#   r::NeutralRecurrence{T,ffa,dffa}
#   cfs::Vector{T}
#   r0::T # radius of
