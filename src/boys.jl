# TODO: this probably goes to a separate package.
raw"""
    boys(n::Int, x)

Compute Boys function $F_n (x)$.
"""
function boys(n::Int, x)
    if x < 3.7
        return _boys_small(n, x)
    elseif x > 14.9
        return _boys_large(n, x)
    end

    return _boys_medium(n, x)
end

# Implementation of the Boys function suited for x < 3.7.
# This uses three different Taylor expansions up to order 18 at x = 0.
# This accounts for ~42.9% of the cases.
function _boys_small(n::Int, x)
    if x < 6e-4
        return _boys_taylor0(n, x)
    elseif x < 7e-1
        return _boys_taylor8(n, x)
    end

    return _boys_taylor18(n, x)
end

# Implementation of the Boys function suited for x > 14.9.
# This uses an asymptotic implementation.
# This accounts for ~36.4% of the cases.
function _boys_large(n::Int, x)
    if n > 3 && x > 1e2
        return zero(float(x))
    end

    return _boys_asymptotic(n, x)
end

# Implementation of the Boys function suited for all x.
# This uses implementations depending on n.
# This accounts for ~20.7% of the cases.
function _boys_medium(n::Int, x)
    if n == 0
        return _boys0_erf(x)
    end

    _boys_gamma(n, x)
end

# Implementation of the Boys function using the incomplete gamma function.
function _boys_gamma(n::Int, x::Float64)
    # gamma_inc is much slower for Float16 or Float32
    w = n + 0.5
    return gamma_inc(w, x, 0)[1] * gamma(w) / (2 * x^w)
end
_boys_gamma(n::Int, x::AbstractFloat) = oftype(x, _boys_gamma(n, convert(Float64, x)))
_boys_gamma(n::Int, x) = _boys_gamma(n, float(x))

# Implementation of the Boys function for n = 0 using the error function.
@fastmath _boys0_erf(x) = erf(sqrt(x)) * sqrt(π / x) / 2

# Some Taylor expansions of the Boys function at x = 0
@inline @fastmath _boys_taylor0(n::Int, x) = one(x) / (2n + 1)
@inline @fastmath _boys_taylor1(n::Int, x) = _boys_taylor0(n, x) - x / (2n + 3)
@inline @fastmath _boys_taylor2(n::Int, x) = _boys_taylor1(n, x) + x^2 / (2 * (2n + 5))
@inline @fastmath _boys_taylor3(n::Int, x) = _boys_taylor2(n, x) - x^3 / (6 * (2n + 7))
@inline @fastmath _boys_taylor4(n::Int, x) = _boys_taylor3(n, x) + x^4 / (24 * (2n + 9))
@inline @fastmath _boys_taylor5(n::Int, x) = _boys_taylor4(n, x) - x^5 / (120 * (2n + 11))
@inline @fastmath _boys_taylor6(n::Int, x) = _boys_taylor5(n, x) + x^6 / (720 * (2n + 13))
@inline @fastmath _boys_taylor7(n::Int, x) = _boys_taylor6(n, x) - x^7 / (5040 * (2n + 15))
@inline @fastmath _boys_taylor8(n::Int, x) = _boys_taylor7(n, x) + x^8 / (40320 * (2n + 17))
@inline @fastmath _boys_taylor9(n::Int, x) = _boys_taylor8(n, x) - x^9 / (362880 * (2n + 19))
@inline @fastmath _boys_taylor10(n::Int, x) = _boys_taylor9(n, x) + x^10 / (3628800 * (2n + 21))
@inline @fastmath _boys_taylor11(n::Int, x) = _boys_taylor10(n, x) - x^11 / (39916800 * (2n + 23))
@inline @fastmath _boys_taylor12(n::Int, x) = _boys_taylor11(n, x) + x^12 / (479001600 * (2n + 25))
@inline @fastmath _boys_taylor13(n::Int, x) = _boys_taylor12(n, x) - x^13 / (6227020800 * (2n + 27))
@inline @fastmath _boys_taylor14(n::Int, x) = _boys_taylor13(n, x) + x^14 / (87178291200 * (2n + 29))
@inline @fastmath _boys_taylor15(n::Int, x) = _boys_taylor14(n, x) - x^15 / (1307674368000 * (2n + 31))
@inline @fastmath _boys_taylor16(n::Int, x) = _boys_taylor15(n, x) + x^16 / (20922789888000 * (2n + 33))
@inline @fastmath _boys_taylor17(n::Int, x) = _boys_taylor16(n, x) - x^17 / (355687428096000 * (2n + 35))
@inline @fastmath _boys_taylor18(n::Int, x) = _boys_taylor17(n, x) + x^18 / (6402373705728000 * (2n + 37))
@inline @fastmath _boys_taylor19(n::Int, x) = _boys_taylor18(n, x) - x^19 / (121645100408832000 * (2n + 39))
@inline @fastmath _boys_taylor20(n::Int, x) = _boys_taylor19(n, x) + x^20 / (2432902008176640000 * (2n + 41))

# An asymptotic approximation of the Boys function at infinity
# http://shivupa.github.io/blog/efficient-evaluation-of-the-boys-function
# TODO: this could be rewritten in terms of the gamma function, but is it worth it?
@fastmath _boys_asymptotic(n::Int, x) = doublefactorial(2n - 1) * sqrt(π / (x^(2n + 1))) / (2^(n + 1))
