"""
    _normalization_constant(alpha, l, m, n)

Compute a Gaussian primitive normalization constant.
"""
function _normalization_constant(alpha, l, m, n)
    N = sqrt((2alpha / π)^3)
    N *= (4alpha)^(l + m + n)
    N /= doublefactorial(2l - one(l)) * doublefactorial(2m - one(m)) * doublefactorial(2n - one(n))

    return sqrt(N)
end

"""
    _third_center(alpha, alpha_coord, beta, beta_coord, gamma)

Compute a weighted center for two Gaussian primitives.

This is the center from the Gaussian product theorem. `gamma` must be equal to
`alpha + beta`.
"""
@inline function _third_center(alpha, alpha_coord, beta, beta_coord, gamma)
    return @. (alpha * alpha_coord + beta * beta_coord) / gamma
end

# TODO: this could be merged into the above function and work inplace?
"""
    _gaussian_prod_factor(alpha, alpha_coord, beta, beta_coord, gamma)

Compute the factor from the Gaussian product theorem.
"""
function _gaussian_prod_factor(alpha, alpha_coord, beta, beta_coord, gamma)
    # TODO: the following has an extra unneeded ^2
    return exp(-(alpha * beta / gamma) * dist(alpha_coord, beta_coord)^2)
end

"""
    _ck(k, l, m, a, b)

Compute the c_k auxiliary coefficient.
"""
@fastmath function _ck(k, l, m, a, b)
    res = zero(a)

    # TODO: this can be implemented in a single loop, instead of
    # a double loop plus if
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

"""
    _Si(la, lb, pax, pbx, gamma)

Compute the _Si auxiliary term.
"""
function _Si(la, lb, pax, pbx, gamma)
    res = zero(pax)

    for k in 0:fld(la + lb, 2)
        double_k = 2k
        factor = sqrt(π / gamma) * doublefactorial(double_k - one(k)) / (2gamma)^k
        res += factor * _ck(double_k, la, lb, pax, pbx)
    end

    return res
end

raw"""
    _vlri(l, r, i, la, lb, pax, pbx, pcx, epsilon)

Compute the $v_{lri}$ auxiliary coefficient.
"""
function _vlri(l, r, i, la, lb, pax, pbx, pcx, epsilon)
    m = (-one(l))^(l + i)

    t = r + i
    k = l - 2t
    p = pcx^k * epsilon^t
    f =  factorial(l) / (factorial(r) * factorial(i) * factorial(k))

    return m * p * f * _ck(l, la, lb, pax, pbx)
end
