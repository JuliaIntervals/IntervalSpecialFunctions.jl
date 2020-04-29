# This code was due to Jeffrey Sarnoff. See
# https://github.com/JuliaIntervals/IntervalArithmetic.jl/issues/372

using IntervalArithmetic
import SpecialFunctions
import ArbNumerics
using ArbNumerics: ArbFloat, ArbReal, setinterval, lowerbound, upperbound

# number of bits that suffices to alter the least significant digit
#   use: prevfloat(x, DigitBits), nextfloat(x, DigitBits)
const DigitBits = ceil(Int, log2(10))

# Interval{T} --> ArbReal, thin enclosure
function arbreal(x::Interval{T}) where T
     setprecision(ArbFloat, T)
     lo = prevfloat(ArbFloat(x.lo), DigitBits) # assure least significant    
     hi = nextfloat(ArbFloat(x.hi), DigitBits) #   decimal digit changes
     return setinterval(lo, hi)
end

# T, ArbReal --> Interval{T}
#    P is the working precision from ArbNumerics
function tinterval(::Type{T}, x::ArbReal{P}) where {T, P}
     lo, hi = lowerbound(x), upperbound(x)
     return Interval{T}(T(lo), T(hi))
end

function specialfunction(fn, a::Interval{T}) where T
     x = fn(arbreal(a))
     return tinterval(T, x)
end

function specialfunction(fn, a::Interval{T}, b::Interval{T}) where T
     x = fn(arbreal(a), arbreal(b))
     return tinterval(T, x)
end

function specialfunction(fn, a::Interval{T}, b::Interval{T},
                             c::Interval{T}) where T
     x = fn(arbreal(a), arbreal(b), arbreal(c))
     return tinterval(T, x)
end

function specialfunction(fn, a::Interval{T}, b::Interval{T},
                             c::Interval{T}, d::Interval{T}) where T
     x = fn(arbreal(a), arbreal(b), arbreal(c), arbreal(d))
     return tinterval(T, x)
end
     
Base.setprecision(::Type{ArbFloat}, ::Type{Float64}) =
    setprecision(ArbFloat, 53)
Base.setprecision(::Type{ArbFloat}, ::Type{Float32}) =
    setprecision(ArbFloat, 24)
Base.setprecision(::Type{ArbFloat}, ::Type{BigFloat}) =
    setprecision(ArbFloat, precision(BigFloat))

# Don't import these functions
# functions of 1 arg from Base.Math (for Interval{BigFloat})
#for F in (:acos, :acosh, :acot, :acoth, :acsc, :acsch, :asec,
#          :asech, :asin, :asinh, :atan, :atanh, :cbrt, :cos,
#          :cosc, :cosh, :cospi, :cot, :coth, :csc, :csch,
#          :exp, :expm1, :log, :log10, :log1p, :log2, :sec,
#          :sech, :sin, :sinc, :sinh, :sinpi, :sqrt, :tan, :tanh)
#    @eval Base.$F(x::A) where {T, A<:Interval{BigFloat}} =
#          specialfunction($F,x)
#end

# functions of 2 arg from Base.Math (for Interval{BigFloat})
for F in (:^, :atan, :hypot)
  @eval begin
    Base.$F(x::A, y::A) where {T, A<:Interval{BigFloat}} =
        specialfunction($F, x, y)
    Base.$F(x::A, y::Real) where {T, A<:Interval{BigFloat}} =
        specialfunction($F, x, A(y))
    Base.$F(x::Real, y::A) where {T, A<:Interval{BigFloat}} =
        specialfunction($F, A(x), y)
  end
end
  
# functions of 1 arg exported by SpecialFunctions and ArbNumerics
for F in (:airyai, :airyaiprime, :airybi, :airybiprime,
          :besselj0, :besselj1, :bessely0, :bessely1,
          #:erf, :erfc, :erfi,
          :gamma, :lgamma, :digamma, :trigamma,
          :eta, :zeta)
    @eval SpecialFunctions.$F(x::A) where {T, A<:Interval{T}} =
          specialfunction($F,x)
end

# functions of 1 arg exported by both with distinct names
for (SF, AF) in ((:cosint, :ci), (:sinint, :si),
                 (:ellipk, :elliptic_k), (:ellipe, :elliptic_e))
    @eval SpecialFunctions.$SF(x::A) where {T, A<:Interval{T}} =
          specialfunction($AF,x)
end

# functions of 1 arg exported by ArbNumerics only
for F in (:chi, :shi, :ei, :xi,  :lambertw)
    @eval ArbNumerics.$F(x::A) where {T, A<:Interval{T}} =
          specialfunction($F,x)
end

# functions of a 2 args exported by ArbNumerics only
for F in (:elliptic_pi, :elliptic_pi2, # the 2 arg versions 
          :polylog,
          :hypergeometric_0F1, :regular_hypergeometric_0F1)
    @eval ArbNumerics.$F(x::A, y::A) where {T, A<:Interval{T}} =
          specialfunction($F, x, y)
end

# functions of a 3 args exported by ArbNumerics only
for F in (:elliptic_pi, :elliptic_pi2, # the 3 arg versions
          :hypergeometric_1F1, :regular_hypergeometric_1F1)
    @eval ArbNumerics.$F(a::A, b::A, z::A) where {T, A<:Interval{T}} =
          specialfunction($F, a, b, z)
end

# functions of a 4 args exported by ArbNumerics only
for F in (:hypergeometric_2F1, :regular_hypergeometric_2F1)
    @eval ArbNumerics.$F(a::A, b::A, c::A, z::A) where
                                   {T, A<:Interval{T}} =
          specialfunction($F, a, b, c, z)
end

function __init__()
    # provide padding for prevfloat, nextfloat in arbreal()
    extrabits = ArbNumerics.ExtraBits[] + DigitBits
    ArbNumerics.setextrabits(extrabits)
    return nothing
end
