# REAL INTERVAL NEWTON

@compat struct FHat{T}<:Function
  r::NeutralRecurrence{T}
end
@compat (d::FHat)(x) = x*(1+d.r.fa(x))^d.r.α

@compat struct DFHat{T}<:Function
  r::NeutralRecurrence{T}
end
@compat function (d::DFHat)(x)
  fa = 1+d.r.fa(x)
  faα = fa^d.r.α
  (1+d.r.α*x/fa)*faα
end

function mapinv_trans{T}(r::NeutralRecurrence{Interval{T}},y::Real)
  nr = newton(FHat(r),DFHat(r),Interval{T}(y))
  @assert length(nr) == 1
  nr[1].interval
end
mapinv{T}(r::NeutralRecurrence{Interval{T}},y::Real) = mapinv_trans(r,(Interval{T}(y)-r.p)^r.α*r.sgn)^(1/r.α)*r.sgn + r.p


# COMPLEX/GENERAL NEWTON

function neutral_newton_step{T}(r::NeutralRecurrence{T},y,x)
  fa = 1+r.fa(x)
  faα = fa^r.α
  xx = x-(x-y/faα)/(1+r.α*r.dfa(x)*x/fa)
  xx,abs(x*faα-y)
end


function neutral_newton_trans{T}(r::NeutralRecurrence{T},y,tol=100eps(abs(y)))
  x = T(r.rad)
  K = log(r.dlogdfmax/2)
  rm = zero(T)
  for n = 1:3ceil(Int,log2((log(tol)+K)/(log(r.rad)+K)))
    x,rm = neutral_newton_step(r,y,x)
    rm < tol && break
  end
  rm > tol && error("Newton iteration: failure to converge: x=$x, y=$y, rm=$rm")
  # check nothing's wrong
  x
end

neutral_newton{T}(r::NeutralRecurrence{T},y,tol=20eps(abs(y-r.p))) =
  neutral_newton_trans(r,(y-r.p)^r.α *r.sgn,tol^r.α)^(1/r.α) * r.sgn+r.p


function neutral_newton_trans{T}(r::NeutralRecurrence{Interval{T}},y,tol=20eps(abs(y).hi))
  @assert abs(y) < r.rad
  @assert abs(angle(y)) < r.Ψ
  x = copy(y)
  K = log(r.dlogdfmax/2) # make exact
  for n = 1:roundupbound(log2((log(tol)+K)/(log(r.rad)+K)))
    x,rm = neutral_newton_step(r,y,midinterval(x))
    rm.hi < tol && break
  end
  neutral_newton_step(r,y,x)[1] ∩ x
end

neutral_newton{T}(r::NeutralRecurrence{Interval{T}},y,tol=20eps((abs(y-r.p).hi)^r.α)) =
  neutral_newton_trans(r,(y-r.p)^r.α *r.sgn,tol)^(1/r.α) * r.sgn+r.p

mapinv_trans{T}(r::NeutralRecurrence{T},y::Real) = neutral_newton_trans(r,T(y))
mapinv{T}(r::NeutralRecurrence{T},y::Real) = neutral_newton(r,T(y))
mapinv_trans{T}(r::NeutralRecurrence{T},y::Complex) = neutral_newton_trans(r,Complex{T}(y))
mapinv{T}(r::NeutralRecurrence{T},y::Complex) = neutral_newton(r,Complex{T}(y))
