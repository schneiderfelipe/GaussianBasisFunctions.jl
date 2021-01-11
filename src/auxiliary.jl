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