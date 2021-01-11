function _compute_primitive_norm_constant(alpha, l, m, n)
    N = sqrt((2alpha / π)^3)
    N *= (4alpha)^(l + m + n)
    N /= doublefactorial(2l - 1) * doublefactorial(2m - 1) * doublefactorial(2n - 1)

    return sqrt(N)
end

function _compute_primitive_third_center(alpha, alpha_coord, beta, beta_coord)
    gamma = alpha + beta

    return (alpha * alpha_coord + beta * beta_coord) / gamma
end

# TODO: I think this function can be implemented in a single loop, instead of
# a double loop plus if
function _compute_ck(k, l, m, a, b)
    res = zero(a)

    for i in 0:l
        l_choose_i = binomial(l, i)
        a_at_l_minus_i = a^(l - i)

        for j in 0:m
            if i + j == k
                m_choose_j = binomial(m, j)
                b_at_m_minus_j = b^(m - j)

                res += l_choose_i * a_at_l_minus_i * m_choose_j * b_at_m_minus_j
            end
        end
    end

    return res
end

function _compute_Si(i_a, i_b, pa_i, pb_i, gamma)
    res = zero(pa_i)

    for k in 0:floor(Int, (i_a + i_b) / 2)
        factor = sqrt(π / gamma) * doublefactorial(2k - 1) / (2gamma)^k

        res += factor * _compute_ck(2k, i_a, i_b, pa_i, pb_i)
    end

    return res
end

# This assumes real basis functions, should we have
# `RealBasisFunction <: AbstractBasisFunction`?
function overlap(basis::Vector{<:AbstractBasisFunction})
    n = length(basis)
    # TODO: generalize types
    S = Matrix{Float64}(undef, n, n)

    # TODO: avoid calculating the diagonal
    for j in 1:n, i in 1:j
        # We assume overlap(a, a) returns 1
        S[i, j] = overlap(basis[i], basis[j])
    end

    # We defined the lower triangular part, so we return a symmetric matrix
    return Symmetric(S)
end

# TODO: create a generic, numerically integrated, fallback for any callable
# objects
function overlap(a::GaussianBasisFunction, b::GaussianBasisFunction)
    res = zero(a.coeffs[1])

    for i in 1:length(a)
        alpha = a.alphas[i]
        N_a = _compute_primitive_norm_constant(alpha, a.l, a.m, a.n)

        term = zero(res)
        for j in 1:length(b)
            beta = b.alphas[j]
            N_b = _compute_primitive_norm_constant(beta, b.l, b.m, b.n)

            gamma = alpha + beta

            # TODO: the following has an extra unneeded ^2
            K_p = exp(-(alpha * beta / gamma) * dist(a, b)^2)

            p_coord = _compute_primitive_third_center(alpha, a.coord, beta, b.coord)
            pa = p_coord - a.coord
            pb = p_coord - b.coord
            S_x = _compute_Si(a.l, b.l, pa[1], pb[1], gamma)
            S_y = _compute_Si(a.m, b.m, pa[2], pb[2], gamma)
            S_z = _compute_Si(a.n, b.n, pa[3], pb[3], gamma)

            term += b.coeffs[j] * N_b * K_p * S_x * S_y * S_z
        end

        res += a.coeffs[i] * N_a * term
    end

    return res
end
