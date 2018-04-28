for f in (:erf, :erfc)

    @eval function($f)(x::BigFloat, r::RoundingMode)
        setrounding(BigFloat, r) do
            ($f)(x)
        end
    end

    @eval ($f)(a::Interval{Float64}) = convert(Interval{Float64}, ($f)(big53(a)))

end

function erf(a::Interval{T}) where T
    isempty(a) && return a
    @round( erf(a.lo), erf(a.hi) )
end

function erfc(a::Interval{T}) where T
    isempty(a) && return a
    @round( erfc(a.hi), erfc(a.lo) )
end

for f in (:erfinv, :erfcinv)

    @eval function($f)(x::BigFloat, r::RoundingMode)
        setrounding(BigFloat, r) do
            ($f)(x)
        end
    end
end

erfinv(a::BigFloat) = mid(_erfinv(a))
erfcinv(a::BigFloat) = mid(_erfcinv(a))

function _erfinv(a::T) where T
    domain = Interval{T}(-1, 1)
    a ∉ domain && return DomainError("$a is not in [-1, 1]")
    f = x -> erf(x) - a
    fp = x->2/sqrt(pi_interval(T)) * exp(-x^2)
    rts = roots(f, Interval{T}(-Inf, Inf), deriv=fp)
    @assert length(rts) == 1 # && rts[1].status == :unique

    rts[1].interval
end

function erfinv(a::Interval{T}) where T
    domain = Interval{T}(-1, 1)
    a = a ∩ domain

    isempty(a) && return a
    hull(_erfinv(a.lo), _erfinv(a.hi))
end

function _erfcinv(a::T) where T
    domain = Interval{T}(0, 2)
    a ∉ domain && return DomainError("$a is not in [0, 2]")
    f = x -> erfc(x) - a
    fp = x -> -2/sqrt(pi_interval(T)) * exp(-x^2)
    rts = roots(f, Interval{T}(-Inf, Inf), deriv=fp)
    @assert length(rts) == 1 # && rts[1].status == :unique

    rts[1].interval
end

function erfcinv(a::Interval{T}) where T
    domain = Interval{T}(0, 2)
    a = a ∩ domain

    isempty(a) && return a
    hull(_erfcinv(a.hi), _erfcinv(a.lo))
end
