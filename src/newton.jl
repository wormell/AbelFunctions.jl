# REAL INTERVAL NEWTON

@compat struct FHat{T}
  r::NeutralRecurrence{T}
end
@compat (d::FHat)(x) = x*(1+d.r.fa(x))^d.r.α

@compat struct DFHat{T}
  r::NeutralRecurrence{T}
end
@compat function (d::DFHat)(x)
  fa = 1+d.r.fa(x)
  faα = fa^d.r.α
  (1+d.r.α*x/fa)*faα
end

mapinv_trans{T}(r::NeutralRecurrence{Interval{T}},y::Interval{T}) = newton(FHat(r),DFHat(r),y)
mapinv{T}(r::NeutralRecurrence{Interval{T}},y::Interval{T}) = mapinv_trans(r,(y-r.p)^r.α*r.sgn)^(1/r.α)*r.sgn + r.p


# COMPLEX/GENERAL NEWTON

function neutral_newtown_step{T}(r::NeutralRecurrence{T},y::T,x::T)
  fa = 1+r.fa(x)
  faα = fa^r.α
  x = x-(x-y/faα)/(1+r.α*x/fa)
  x,x*faα-y
end


function neutral_newton_trans{T}(r::NeutralRecurrence{T},y::T,tol=20eps(abs(y)))
  x = zero(T)
  K = log(r.dlogdfmax/2)
  for n = 1:2ceil(Int,log((log(r.rad)+K)/(log(tol)+K)))
    x,rm = neutral_newtown_step(r,y,x)
    abs(rm) < tol && break
  end
  # check nothing's wrong
  x
end

neutral_newton{T<:Real}(r::NeutralRecurrence{T},y::T,tol=20eps(abs(y-r.p)^α)) =
  neutral_newton_trans(r,(y-r.p)^r.α *r.sgn,tol)^(1/r.α) * r.sgn+r.p

for op in (:∩,:∪)
  @eval @compat ($op){T<:Interval,S<:Interval}(x::Complex{T},y::Complex{S}) =
    complex(($op)(real(x),real(y)),($op)(imag(x),imag(y)))
end

function neutral_newton_trans{T}(r::NeutralRecurrence{Interval{T}},y,tol=20eps(abs(y).hi))  @assert
  @assert abs(y) < r.rad
  @assert abs(arg(y)) < r.Ψ
  x = zero(Interval(T))
  K = log(r.dlogdfmax/2) # make exact
  for n = 1:ceil(Int,log((log(r.rad)+K)/(log(tol)+K)))
    x = Interval(mid(x))
    x,rm = neutral_newton_step(r,y,x)
    abs(rm).hi < tol && break
  end
  neutral_newton_step(r,y,x) ∩ x
  end
end

neutral_newton{T}(r::NeutralRecurrence{T},y::T,tol=20eps(abs(y-r.p))) =
  neutral_newton_trans(r,(y-r.p)^r.α *r.sgn,tol^r.α)^(1/r.α) * r.sgn+r.p

mapinv_trans{T}(r::NeutralRecurrence{T},y) = neutral_newton_trans(r,T(y))
mapinv{T}(r::NeutralRecurrence{T},y) = neutral_newton(r,T(y))
