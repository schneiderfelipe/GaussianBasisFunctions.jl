"""
    gaussian_norm(α, l, m, n)

Compute the normalization constant for a Gaussian primitive basis function.
"""
@inline function gaussian_norm(α, l::Int, m::Int, n::Int)
    N = sqrt((2α / π)^3) * (4α)^(l + m + n)
    κ = doublefactorial(2l - 1) * doublefactorial(2m - 1) * doublefactorial(2n - 1)

    return sqrt(N / κ)
end

"""
    gpt!(p_coord, α, α_coord, β, β_coord)

Compute the relevant quantities from the Gaussian product theorem for two Gaussian primitives.

`p_coord` receives the coordinate of the new Gaussian center. `γ = α + β` and
`Kₚ`, the center from the Gaussian product theorem, are returned as values.
"""
@inline function gpt!(p_coord::AbstractVector, α::Real, α_coord::AbstractVector, β::Real, β_coord::AbstractVector)
    γ = α + β

    # TODO: the following has an extra unneeded ^2
    Kₚ = exp(-(α * β / γ) * dist(α_coord, β_coord)^2)

    @. p_coord = (α * α_coord + β * β_coord) / γ
    return Kₚ, γ
end

raw"""
    c(k, l, m, a, b)

Compute the $c_k$ auxiliary coefficient.
"""
@inline @fastmath function c(k, l, m, a, b)
    I = zero(a)

    # TODO: this can be implemented in a single loop, instead of
    # a double loop plus if, which could be faster.
    for i in 0:l
        l_choose_i = binomial(l, i)
        a_at_l_minus_i = a^(l - i)

        for j in 0:m
            if i + j == k
                m_choose_j = binomial(m, j)
                b_at_m_minus_j = b^(m - j)

                I += l_choose_i * a_at_l_minus_i * m_choose_j * b_at_m_minus_j
            end
        end
    end

    return I
end

"""
    Sᵢ(la, lb, pax, pbx, γ)

Compute the Sᵢ auxiliary term.
"""
@inline @fastmath function Sᵢ(la, lb, pax, pbx, γ)
    I = zero(pax)

    for k in 0:fld(la + lb, 2)
        double_k = 2k
        κ = sqrt(π / γ) * doublefactorial(double_k - 1) / (2γ)^k
        I += κ * c(double_k, la, lb, pax, pbx)
    end

    return I
end

raw"""
    θ(l, la, lb, a, b, γ)

Compute the $\theta$ auxiliary coefficient.
"""
@inline @fastmath function θ(l, la, lb, a, b, r, γ)
    κ = factorial(l) / (factorial(r) * factorial(l - 2r))
    return κ * c(l, la, lb, a, b) * γ^(r - l)
end

raw"""
    v(indices::Tuple{Int, Int, Int}, la, lb, pax, pbx, pcx, ϵ)

Compute the $v_{lri}$ auxiliary coefficient.
"""
@inline @fastmath function v(indices::Tuple{Int, Int, Int}, la, lb, pax, pbx, pcx, ϵ)
    l, r, i = indices

    m = (-1)^(l + i)

    t = r + i
    k = l - 2t
    p = pcx^k * ϵ^t

    κ = factorial(l) / (factorial(r) * factorial(i) * factorial(k))
    f = oftype(ϵ, κ)

    return m * p * f * c(l, la, lb, pax, pbx)
end

raw"""
    g(indices::Tuple{Int, Int, Int, Int, Int}, la, lb, lc, ld, pax, pbx, qcx, qdx, pqx, γp, γq)

Compute the $g_{l_P l_Q r_P r_Q i}$ auxiliary coefficient.
"""
@inline @fastmath function g(indices::Tuple{Int, Int, Int, Int, Int}, la, lb, lc, ld, pax, pbx, qcx, qdx, pqx, γp, γq)
    lp, lq, rp, rq, i = indices

    δ = 1 / 4γp + 1 / 4γq

    ℓ = lp + lq
    m₁ = 2 * (rp + rq)
    m₂ = ℓ - m₁
    m₃ = m₂ - 2i

    κ = factorial(m₂) / (factorial(i) * factorial(m₃))
    I = (-1)^(lp + i) * θ(lp, la, lb, pax, pbx, rp, γp) * θ(lq, lc, ld, qcx, qdx, rq, γq)
    return (κ * δ^(m₁ + i - ℓ) * pqx^m₃ / 2^(2ℓ - m₁)) * I
end
