module IntervalSpecialFunctions

using IntervalArithmetic
using SpecialFunctions
using IntervalRootFinding

import IntervalArithmetic: @round, big53

import SpecialFunctions:
    erf, erfc,
    erfinv, erfcinv


include("erf.jl")



end # module
