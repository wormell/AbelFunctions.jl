@compat immutable AbelFunction{T<:Real,ffa,dffa}
  r::NeutralRecurrence{T,ffa,dffa}
  # rad::T # radius of accuracy
  coeffs::Vector{T}
  offset::T
  n::Int # order of approximation (≥ncall, TODO: assert this)
  nrad::T # radius of there being bounds given n
  noptrad::T # radius of exponential accuracy given n
  ncall::Int # (=length(coeffs) -2)
  callrad::T # radius of callability (≥noptrad, TODO: assert this)
  compresserr::T # compression error on radius of callability
  cr::T
  E1::T
  E2::T
end
# useful constants
cr(r::NeutralRecurrence) = min(1,r.dfa0/r.F)/4
function abelerrorconstant1(r::NeutralRecurrence)
  @assert 1-r.F*r.rad > 0
  r.α*r.F/r.M*(1-r.F*r.rad)^(-r.α-1)
end
abelerrorconstant2(r::NeutralRecurrence) = (2r.F*cr(r)/3r.dfa0)

## CONSTRUCTORS

# constructing given coefficients and offset
function AbelFunction{T,ffa,dffa}(r::NeutralRecurrence{T,ffa,dffa},coeffs::Vector{T},offset::T;n=length(coeffs)-2,compress=true)
  # n = length(coeffs)-2
  nrad = min(log(T(2))/(r.F*(r.α*n-1)),r.rad)

  noptrad = nrad*cr(r)/Base.e
  if compress
    callrad=4noptrad/3 # default??
    callrad_xtra=3callrad/2 # for abelerrorD
    @compat cs = reverse(cumsum(reverse((1:n).*abs.(coeffs[3:n+2]).*callrad_xtra.^(1:n))))
    ncall = findlast(cs.>prec(T))
    compresserr = cs[ncall+1]/callrad_xtra^(ncall+1)/(ncall+1)
    # coeffss = view(coeffs,1:ncall+2)
  else
    callrad=nrad
    ncall = n
    compresserr = zero(T)
    # coeffss = coeffs
  end

  AbelFunction{T,ffa,dffa}(r,coeffs[1:ncall+2],offset,n,nrad,noptrad,ncall,callrad,compresserr,
  cr(r),abelerrorconstant1(r),abelerrorconstant2(r))
end

# constructing given coefficients and basepoint
function AbelFunction{T,ffa,dffa}(r::NeutralRecurrence{T,ffa,dffa},coeffs::Vector{T};basepoint::Real=unhat(r.rad,r),compress=true)
  @assert basepoint > 0
  bp = hat(basepoint,r)
  a = AbelFunction(r,coeffs,zero(T),compress=false)
  ctr = 0
  while ~(bp < a.noptrad)
    bp = mapinv_trans(r,bp)
    ctr += 1
  end
  offset = -(map_trans(a,bp)-ctr)
  AbelFunction(r,coeffs,offset,compress=compress)
end

# constructing given order and basepoint
function AbelFunction{T,ffa,dffa}(r::NeutralRecurrence{T,ffa,dffa},n::Int;basepoint::Real=unhat(r.rad,r),compress=true)
  s = min(log(T(2))/(r.F*(r.α*n-1)),r.rad/2)
  Rs = max(r.rad,2s)

  xn = n*Rs*r.F*r.α*(1+r.F*Rs)^(r.α*n-1)
  Ninterp = nextpow2(max(n+2,roundupbound(-log(prec(T)/xn)/log(Rs/s))))
  wvals = Array{Complex{T}}(Ninterp)
  Wvals = Array{Complex{T}}(Ninterp,n+1)
  ẑ = complex(s);
  twid = exp(im*2T(pi)/Ninterp)
  for i = 1:Ninterp
    faẑ = (1+r.fa(ẑ))^r.α
    wvals[i] = ((1-1/faẑ)/(r.dfa0*r.α*ẑ) - 1)/ẑ
    # wvals[Ninterp-i+1] = conj(wvals[i])
    Wvals[i,1] = log(faẑ)/ẑ
    # Wvals[Ninterp-i+1,1] = conj(Wvals[i,1])

    faẑpow = one(Complex{T})
    for j = 1:n
      faẑpow *= faẑ
      Wvals[i,j+1] = (faẑpow-1)/ẑ
      # Wvals[Ninterp-i+1,j+1] = conj(faẑpow-1)
    end
    ẑ *= twid
  end
  MRsw = Rs*r.F/2r.dfa0*(1-r.F*Rs)^(-r.α-1)
  w = real(taylortransform(wvals,s,MRsw,Rs)[1:n+1])

  W = zeros(T,n+1,n+1)

  MRsW0 = Rs*r.F/(1-r.F*Rs)
  W[:,1] = real(taylortransform(Wvals[:,1],s,MRsW0,Rs)[1:n+1])

  MRsWtwid = (1+r.F*Rs)^r.α
  MRsW = Rs*r.F*r.α/(1+r.F*Rs)
  for i = 1:n
    MRsW *= MRsWtwid
    W[i+1:end,i+1] = real(taylortransform(Wvals[:,i+1],s,i*MRsW,Rs)[1:n+1-i])
  end
  AbelFunction(r,[1/(r.dfa0*r.α);W\w],basepoint=basepoint,compress=compress)#,W,w
end

# constructing given only basepoint
function AbelFunction{T,ffa,dffa}(r::NeutralRecurrence{T,ffa,dffa};basepoint::Real=unhat(r.rad,r),tol=prec(T),compress=true)
  ec = abelerrorconstant1(r)*(1+abelerrorconstant2(r))
  n = roundupbound(-log(tol/ec)-1)
  AbelFunction(r,n;basepoint=basepoint,compress=compress)
end

# TODO: #PRUNING COEFFICIENTS
# function prunecoeffs!(a::AbelFunction{T,ffa,dffa},tol=prec(T))
#   cv = @. a.coeffs[3:end]*a.callablerad^(1:a.n)

## CALLING THE FUNCTION

# # TODO: ascertain 1. if this is right 2. decide if need rad to be the same as well
# (==)(a1::AbelFunction,a2::AbelFunction) = (a1.offset==a2.offset)&&(a1.r == a2.r)

function abelerror{T}(a::AbelFunction{T},ẑ,n=a.n)
  Reẑii = real(ẑ)*(1+(imag(ẑ)/real(ẑ))^2)
  a.E1*(1+a.E2)*(Reẑii/(a.cr*a.nrad))^(n+1) + abs(ẑ)^(a.ncall+1)*a.compresserr
end
function abelerrorD{T}(a::AbelFunction{T},ẑ,n=a.n)
  absẑ = abs(ẑ)
  ẑanalytic_edge = 3/(4/ẑ-2/absẑ)
  2/absẑ*erroradjustment(T,abelerror(a,ẑanalytic_edge,n))
  #includes: + (a.ncall+1)*abs(ẑ)^(a.ncall+1)*a.compresserr
end

@compat (a::AbelFunction)(z) = map_trans(a,hat(z,a.r))
function map_trans(a::AbelFunction,ẑ)
  @assert abs(ẑ) ≤ a.callrad
  R = a.coeffs[1]/ẑ + a.coeffs[2]*log(ẑ)
  ẑpow = copy(ẑ)
  for i = 1:a.ncall
    R += a.coeffs[i+2]*ẑpow
    # println(a.coeffs[i+2]*ẑpow)
    ẑpow *= ẑ
  end
  R += erroradjustment(abelerror(a,ẑ))
  R += a.offset
  R
end

function mapD{T}(a::AbelFunction{T},z)
  ẑ = hat(z,a.r)
  mapD_trans(a,ẑ)*a.r.α*ẑ/(z-a.r.p) # sgn?
end

function mapD_trans{T}(a::AbelFunction{T},ẑ)
  @assert abs(ẑ) ≤ a.callrad
  dR = (-a.coeffs[1]/ẑ + a.coeffs[2])/ẑ
  ẑpow = one(typeof(ẑ))
  for i = 1:a.ncall
    dR += i*a.coeffs[i+2]*ẑpow
    ẑpow *= ẑ
  end
  dR += erroradjustment(T,abelerrorD(a,ẑ))
  dR
end

function mapP{T}(a::AbelFunction{T},z)
  zt = a.r.sgn*(z-a.r.p)
  ẑ = zt^a.r.α
  R,dR = map_transP(a,ẑ)
  R,dR*a.r.α*ẑ/zt
end
function mapP_trans{T}(a::AbelFunction{T},ẑ)
  @assert abs(ẑ) ≤ a.callrad
  R = a.coeffs[1]/ẑ + a.coeffs[2]*log(ẑ)
  dR = -a.coeffs[1]/ẑ
  dR = (dR + a.coeffs[2])/ẑ
  ẑpow = one(typeof(ẑ))
  for i = 1:a.ncall
    pe = a.coeffs[i+2]*ẑpow
    dR += i*pe
    R += pe*ẑ
    ẑpow *= ẑ
  end
  dR += erroradjustment(T,abelerrorD(a,ẑ))
  R += erroradjustment(T,abelerror(a,ẑ))
  R += a.offset
  R,dR
end
