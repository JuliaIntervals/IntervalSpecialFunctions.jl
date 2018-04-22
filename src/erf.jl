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

function erfinv(a::Interval{T}) where T
    domain = Interval{T}(-1, 1)
    a = a ∩ domain

    isempty(a) && return a
    roots(x -> erf(x) - a, domain, Newton(x->2/√(pi_interval(T)) * exp(-x^2)))
end

function erfcinv(a::Interval{T}) where T
    domain = Interval{T}(0, 2)
    a = a ∩ domain

    isempty(a) && return a
    roots(x -> erfc(x) - a, domain, Newton(x->-2/√(pi_interval(T)) * exp(-x^2)))
end
