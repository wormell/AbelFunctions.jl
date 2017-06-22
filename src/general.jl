for op in (:intersect,:union)
  @eval @compat (Base.$op){T<:Interval,S<:Interval}(x::Complex{T},y::Complex{S}) =
    complex((Base.$op)(real(x),real(y)),(Base.$op)(imag(x),imag(y)))
end
# mid{T}(x::Complex{Interval{T}}) = complex(mid(real(x)),mid(imag(x)))

function midinterval(x::Interval)
  md = mid(x);
  Interval(md,md)
end
midinterval{T<:Interval}(z::Complex{T}) = complex(midinterval(real(z)),midinterval(imag(z)))


#@compat typealias RealComplex{T} Union{T,Complex{T}}

# from ValidatedTransfer
function upbound(x)
  xu = _upbound(x)
  assert(xu≥0)
  xu
end
_upbound(x::Real) = x
_upbound(x::Interval) = x.hi

roundupbound(x) = ceil(Int,upbound(x))
