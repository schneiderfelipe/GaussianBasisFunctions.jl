"""
	doublefactorial(n)

Computes the exact
[double factorial n!!](https://en.wikipedia.org/wiki/Double_factorial).

This implementation is not suited for large n (n <= 33) as it might overflow.
This implementation might be memoized in the future.

See [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) for
another implementation.
"""
function doublefactorial(n)
	if n < -1
		return zero(n)
	elseif n < 1
		return one(n)
	elseif n < 4
		return n
	end

	return n * doublefactorial(n - one(n) - one(n))
end

"""
	dist(x, y)

Compute the distance between vectors or Gaussian basis functions.

The base implementation might change in the future, or use anoter package for
it.
"""
dist(x, y) = sqrt(sum((xi - yi)^2 for (xi, yi) in zip(x, y)))
dist(a::GaussianBasisFunction, b::GaussianBasisFunction) = dist(a.coord, b.coord)
