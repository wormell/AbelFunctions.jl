function abelfunction_newton_step(a::AbelFunction{T},y,x) where T
  R,dR = mapP_trans(a,x)
  rem = R-y
  xx = x-rem/dR
  xx,abs(rem)
end

abelfunction_newton_initialguess(a::AbelFunction{T},y) where T = a.coeffs[1]/y
function abelfunction_newton_trans(a::AbelFunction{T},y,tol=100eps(abs(y))) where T
  # @assert abs(y) < r.rad
  # @assert abs(angle(y)) < r.Ψ
  x = abelfunction_newton_initialguess(a,y)
  K = log(T(3)) # TODO: better
  for n = 1:3ceil(Int,log2((log(tol)+K)/(log(a.noptrad)+K)))
    x,rm = abelfunction_newton_step(a,y,x)
    # rm < tol && break
  end
  # rm > tol && error("Newton iteration: failure to converge: x=$x, y=$y, rm=$rm")
  # check nothing's wrong
  x += erroradjustment(typeof(x),tol)
  x
end

function abelfunction_newton_trans(a::AbelFunction{Interval{T}},y,tol=20eps((1/abs(y)).hi)) where T
  # @assert abs(y) < r.rad
  # @assert abs(angle(y)) < r.Ψ
  x = abelfunction_newton_initialguess(a,y)
  K = log(T(20)) # TODO: better
  for n = 1:2roundupbound(log2((log(tol)+K)/(log(a.callrad)+K)))
    x,rm = abelfunction_newton_step(a,y,x)
    # rm.hi < tol && break
  end
  x += erroradjustment(typeof(x),tol)
  abelfunction_newton_step(a,y,midinterval(x))[1] ∩ x
end

abelfunction_newton(a::AbelFunction{Interval{T}},y,tol=20eps((abs(y-a.r.p).hi)^(a.r.α))) where T =
  unhat(abelfunction_newton_trans(a,y,tol),a.r)

mapinv_trans(a::AbelFunction{T},y::Real) where T = abelfunction_newton_trans(a,T(y))
mapinv(a::AbelFunction{T},y::Real) where T = abelfunction_newton(a,T(y))
mapinv_trans(a::AbelFunction{T},y::Complex) where T = abelfunction_newton_trans(a,Complex{T}(y))
mapinv(a::AbelFunction{T},y::Complex) where T = abelfunction_newton(a,Complex{T}(y))
