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

function erfinv(a::BigFloat)
    domain = Interval{BigFloat}(-1, 1)
    a ∉ domain && return DomainError("$a is not in [-1, 1]")
    f = x -> erf(x) - a
    fp = x->2/sqrt(pi_interval(BigFloat)) * exp(-x^2)
    rts = roots(f, Interval{BigFloat}(-Inf, Inf), Newton(f, fp))
    @assert length(rts) == 1 && rts[1].status == :unique

    mid(rts[1].interval)
end

function erfinv(a::Interval{T}) where T
    domain = Interval{T}(-1, 1)
    a = a ∩ domain

    isempty(a) && return a
    Interval{T}(erfinv(a.lo), erfinv(a.hi))
end

function erfcinv(a::BigFloat)
    domain = Interval{BigFloat}(0, 2)
    a ∉ domain && return DomainError("$a is not in [0, 2]")
    f = x -> erfc(x) - a
    fp = x -> -2/sqrt(pi_interval(BigFloat)) * exp(-x^2)
    rts = roots(f, Interval{BigFloat}(-Inf, Inf), Newton(f, fp))
    @assert length(rts) == 1 && rts[1].status == :unique

    rts[1].interval
end

function ercfinv(a::Interval{T}) where T
    domain = Interval{T}(0, 2)
    a = a ∩ domain

    isempty(a) && return a
    Interval{T}(erfinv(a.hi), erfinv(a.lo))
end
