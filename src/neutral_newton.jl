# REAL INTERVAL NEWTON

@compat struct FHat{T}
  r::NeutralRecurrence{T}
  val::T
end
@compat (d::FHat)(ẑ) = ẑ*(1+d.r.fa(ẑ))^d.r.α - d.val

@compat struct DFHat{T}
  r::NeutralRecurrence{T}
end
@compat (d::DFHat)(ẑ) = mapD_trans(d.r,ẑ)

function mapinv_trans(r::NeutralRecurrence{Interval{T}},y::Interval) where T
  nr = IntervalRootFinding_newton(FHat(r,y),DFHat(r),Interval{T}(0,y.hi))
  # @assert length(nr) == 1
  # nr[1].interval
  nr
end
mapinv(r::NeutralRecurrence{Interval{T}},y::Real) where T = unhat(mapinv_trans(r,hat(Interval{T}(y),r)),r)

function IntervalRootFinding_newton(f,df,x::Interval{T}) where T#,tol=20eps(T))
  ctr = 0
  while true #diam(x) > tol
    xm = mid(x)
    xn = xm - f(xm)./df(x)
    x,xo = x ∩ xn,x
    xo == x && break
    ctr += 1
    ctr > 100 && error("Failure to converge")
  end
  x
end




# COMPLEX/GENERAL NEWTON

function neutral_newton_step(r::NeutralRecurrence{T},y,x) where T
  fa = 1+r.fa(x)
  faα = fa^r.α
  xx = x-(x-y/faα)/(1+r.α*r.dfa(x)*x/fa)
  xx,abs(x*faα-y)
end


function neutral_newton_trans(r::NeutralRecurrence{T},y,tol=100eps(abs(y))) where T
  x = T(r.rad)
  K = log(r.dlogdfmax/2)
  rm = zero(T)
  for n = 1:3ceil(Int,log2((log(tol)+K)/(log(r.rad)+K)))
    x,rm = neutral_newton_step(r,y,x)
    rm < tol && break
  end
  rm > tol && error("Newton iteration: failure to converge: x=$x, y=$y, rm=$rm")
  # check nothing's wrong
  x += erroradjustment(typeof(x),tol)
  x
end

neutral_newton(r::NeutralRecurrence{T},y,tol=20eps(abs(y-r.p))) where T =
  unhat(neutral_newton_trans(r,hat(y,r),tol^r.α),r)


function neutral_newton_trans(r::NeutralRecurrence{Interval{T}},y,tol=20eps(abs(y).hi)) where T
  @assert abs(y) < r.rad
  @assert abs(angle(y)) < r.Ψ
  x = copy(y)
  K = log(r.dlogdfmax/2) # make exact
  for n = 1:roundupbound(log2((log(tol)+K)/(log(r.rad)+K)))
    x,rm = neutral_newton_step(r,y,midinterval(x))
    rm.hi < tol && break
  end
  x += erroradjustment(typeof(x),tol)
  neutral_newton_step(r,y,x)[1] ∩ x
end

neutral_newton(r::NeutralRecurrence{Interval{T}},y,tol=20eps((abs(y-r.p).hi)^r.α)) where T =
  unhat(neutral_newton_trans(r,hat(y,r),tol),r)

mapinv_trans(r::NeutralRecurrence{T},y::Real) where T = neutral_newton_trans(r,convert(T,y))
mapinv(r::NeutralRecurrence{T},y::Real) where T = neutral_newton(r,convert(T,y))
mapinv_trans(r::NeutralRecurrence{T},y::Complex) where T = neutral_newton_trans(r,convert(Complex{T},y))
mapinv(r::NeutralRecurrence{T},y::Complex) where T = neutral_newton(r,convert(Complex{T},y))
