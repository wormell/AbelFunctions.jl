# import Base.intersect, Base.union, Base.log, Base.hypot

# for op in (:intersect,:union)
#   @eval @compat (Base.$op){T<:Interval,S<:Interval}(x::Complex{T},y::Complex{S}) =
#     complex((Base.$op)(real(x),real(y)),(Base.$op)(imag(x),imag(y)))
# end
# mid{T}(x::Complex{Interval{T}}) = complex(mid(real(x)),mid(imag(x)))
# Base.log(z::Complex{T}) where {T<:Interval} = complex(log(abs(z)),atan2(imag(z),real(z)))

# VERSION < v"0.5" && (Base.hypot(x::Interval,y::Interval) = sqrt(x*x+y*y))


function midinterval(x::Interval)
  md = mid(x);
  Interval(md,md)
end
midinterval(z::Complex{T}) where {T<:Interval} = complex(midinterval(real(z)),midinterval(imag(z)))

#@compat typealias RealComplex{T} Union{T,Complex{T}}

# from ValidatedTransfer
function upbound(x)
  xu = _upbound(x)
  @assert xu≥0
  xu
end
_upbound(x::Real) = x
_upbound(x::Interval) = x.hi

roundupbound(x) = ceil(Int,upbound(x))

prec(::Type{Interval{T}}) where T = eps(T)
prec(T) = eps(T)

erroradjustment(T,err) = zero(T)
erroradjustment(::Type{<:Interval},err) = Interval(-upbound(err),upbound(err))
erroradjustment(::Type{Complex{IT}},err) where {IT<:Interval} =
  complex(erroradjustment(IT,err),erroradjustment(IT,err))
erroradjustment(err) = erroradjustment(typeof(err),err)
