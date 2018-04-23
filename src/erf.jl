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

function _erfinv(a::Float64)
    domain = -1..1
    a ∉ domain && return DomainError("$a is not in [-1, 1]")
    rts = roots(x -> erf(x) - a, -Inf..Inf,
        Newton(x->2/sqrt(pi_interval(Float64)) * exp(-x^2)))
    @assert length(rts) == 1 && rts[1].status == :unique

    rts[1].interval
end

function erfinv(a::BigFloat)
    domain = Interval{BigFloat}(-1, 1)
    a ∉ domain && return DomainError("$a is not in [-1, 1]")
    rts = roots(x -> erf(x) - a, Interval{BigFloat}(-Inf, Inf),
        Newton(x->2/sqrt(pi_interval(BigFloat)) * exp(-x^2)))
    @assert length(rts) == 1 && rts[1].status == :unique

    rts[1].interval
end

function erfinv(a::Interval{BigFloat})
    domain = Interval{BigFloat}(-1, 1)
    a = a ∩ domain

    isempty(a) && return a
    hull(erfinv(a.lo), erfinv(a.hi))
end

function erfinv(a::Interval{Float64})
    domain = Interval{Float64}(-1, 1)
    a = a ∩ domain

    isempty(a) && return a
    hull(_erfinv(a.lo), _erfinv(a.hi))
end

function _erfcinv(a::Float64)
    domain = 0..2
    a ∉ domain && return DomainError("$a is not in [-1, 1]")
    rts = roots(x -> erfc(x) - a, -Inf..Inf,
        Newton(x->-2/sqrt(pi_interval(Float64)) * exp(-x^2)))
    @assert length(rts) == 1 && rts[1].status == :unique

    rts[1].interval
end

function erfcinv(a::BigFloat)
    domain = Interval{BigFloat}(0, 2)
    a ∉ domain && return DomainError("$a is not in [-1, 1]")
    rts = roots(x -> erfc(x) - a, Interval{BigFloat}(-Inf, Inf),
        Newton(x->-2/sqrt(pi_interval(BigFloat)) * exp(-x^2)))
    @assert length(rts) == 1 && rts[1].status == :unique

    rts[1].interval
end

function erfcinv(a::Interval{BigFloat})
    domain = Interval{BigFloat}(0, 2)
    a = a ∩ domain

    isempty(a) && return a
    hull(erfcinv(a.hi), erfcinv(a.lo))
end

function erfcinv(a::Interval{Float64})
    domain = Interval{Float64}(0, 2)
    a = a ∩ domain

    isempty(a) && return a
    hull(_erfcinv(a.hi), _erfcinv(a.lo))
end
