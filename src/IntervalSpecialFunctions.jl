module IntervalSpecialFunctions

using IntervalArithmetic
using SpecialFunctions

import IntervalArithmetic: @round, big53

import SpecialFunctions:
    erf, erfc


include("erf.jl")



end # module
