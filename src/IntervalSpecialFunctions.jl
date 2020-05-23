module IntervalSpecialFunctions

using IntervalArithmetic
using SpecialFunctions
using IntervalRootFinding

import IntervalArithmetic: @round, big53

import SpecialFunctions:
    erf, erfc, erfinv, erfcinv

export 
    # functions of 1 arg exported by SpecialFunctions and ArbNumerics
    airyai, airyaiprime, airybi, airybiprime,
    besselj0, besselj1, bessely0, bessely1,
    erf, erfc, erfinv, erfcinv,
    gamma, lgamma, digamma, trigamma,
    eta, zeta,

    # functions of 1 arg exported by both with distinct names
    cosint, ci, sinint, si, ellipk, elliptic_k, ellipe, elliptic_e,

    # functions of 1 arg exported by ArbNumerics only
    chi, shi, ei, xi,  lambertw,

    # functions of a 2 args exported by ArbNumerics only
    elliptic_pi, elliptic_pi2, # the 2 arg versions 
    polylog,
    hypergeometric_0F1, regular_hypergeometric_0F1,

    # functions of a 3 args exported by ArbNumerics only
    elliptic_pi, elliptic_pi2, # the 3 arg versions
    hypergeometric_1F1, regular_hypergeometric_1F1,

    # functions of a 4 args exported by ArbNumerics only
    hypergeometric_2F1, regular_hypergeometric_2F1

# Cannot make erfinv and erfcinv work
include("erf.jl")

include("arb.jl")


end # module
