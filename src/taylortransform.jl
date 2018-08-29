# from FastTransforms
function backend_fft_pow2!(x::Vector{T}) where T
    n,big2,bigpi=length(x),2one(T),convert(T,pi)
    nn,j=n÷2,1
    for i=1:2:n-1
        if j>i
            x[j], x[i] = x[i], x[j]
            x[j+1], x[i+1] = x[i+1], x[j+1]
        end
        m = nn
        while m ≥ 2 && j > m
            j -= m
            m = m÷2
        end
        j += m
    end
    logn = 2
    while logn < n
        piθ=-bigpi*big2/logn
        wtemp = sin(piθ/2)
        wpr, wpi = -2wtemp^2, sin(piθ)
        wr, wi = one(T), zero(T)
        for m=1:2:logn-1
            for i=m:2logn:n
                j=i+logn
                mixr, mixi = wr*x[j]-wi*x[j+1], wr*x[j+1]+wi*x[j]
                x[j], x[j+1] = x[i]-mixr, x[i+1]-mixi
                x[i], x[i+1] = x[i]+mixr, x[i+1]+mixi
            end
            wr = (wtemp=wr)*wpr-wi*wpi+wr
            wi = wi*wpr+wtemp*wpi+wi
        end
        logn = logn << 1
    end
    return x
end

function fft_pow2(x::Vector)
    y = FastTransforms.interlace(real(x),imag(x))
    backend_fft_pow2!(y)
@compat    return complex.(y[1:2:end],y[2:2:end])
end

taylorinterpolationerror(r,Ninterp::Int,M,R) = M*(r/R)^Ninterp/(1-r/R)

function taylortransform(vals::Vector{Complex{T}},r::T,M::T,R::T) where T
  Ninterp = length(vals)
  cfs = fft_pow2(vals)/Ninterp
  err= erroradjustment(Complex{T},taylorinterpolationerror(r,Ninterp,M,R))
  powr = one(T)
  for i = 1:Ninterp
    cfs[i] *= powr
    cfs[i] += err
    err /= R
    powr /= r
  end
  cfs
end
