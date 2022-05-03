"""
    doublefactorial(n)

Compute the exact
[double factorial n!!](https://en.wikipedia.org/wiki/Double_factorial).

This implementation is not suited for large n (n <= 33) as it might overflow.
This implementation might be memoized in the future.

See [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) for
another implementation.
"""
# NOTE: consider increasing the precomputations to 20
@inline function doublefactorial(n::Int)
    if n < -1
        return zero(n)
    elseif n < 1
        return 1
    elseif n < 4
        return n
    elseif n == 4
        return 8
    elseif n == 5
        return 15
    elseif n == 6
        return 48
    elseif n == 7
        return 105
    elseif n == 8
        return 384
    elseif n == 9
        return 945
    elseif n == 10
        return 3840
    end

    return n * doublefactorial(n - 2)
end

"""
    dist(x, y)

Compute the distance between vectors or Gaussian basis functions.

The base implementation might change in the future, or use anoter package for
it.
"""
@inline function dist(x, y)
    s = zero(x[1])
    for i in 1:length(x)
        s += (x[i] - y[i])^2
    end
    return sqrt(s)
end
@inline dist(a::GaussianBasisFunction, b::GaussianBasisFunction) = dist(a.coord, b.coord)

# Find the triangular root of x
# See also this trick: https://math.stackexchange.com/a/699050/408160
@inline _trrt(x) = (sqrt(8x + one(x)) - one(x)) / 2

# Create a Gaussian basis function with modified momenta
@inline adjust_momenta(g, Δl, Δm, Δn) = GaussianBasisFunction(
    g.coord, g.alphas, g.coeffs, g.l + Δl, g.m + Δm, g.n + Δn
)
