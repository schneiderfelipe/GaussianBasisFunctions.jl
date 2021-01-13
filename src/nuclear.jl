# TODO: this might benefic from merging with oei(..., ..., ::KineticOperator) and
# calling a kernel.
"""
    oei(a::GaussianBasisFunction, b::GaussianBasisFunction, operator::NuclearOperator)

Compute a nuclear matrix element for a pair of basis functions.
"""
function oei(a::GaussianBasisFunction, b::GaussianBasisFunction, operator::NuclearOperator)
    res = zero(a.coeffs[1])
    for i in 1:length(a)
        alpha = a.alphas[i]
        N_a = _normalization_constant(alpha, a.l, a.m, a.n)

        term = zero(res)
        for j in 1:length(b)
            beta = b.alphas[j]
            N_b = _normalization_constant(beta, b.l, b.m, b.n)

            term += b.coeffs[j] * N_b * _nuclear_kernel(
                operator, alpha, a.coord, a.l, a.m, a.n, beta, b.coord, b.l, b.m, b.n
            )
        end

        res += a.coeffs[i] * N_a * term
    end

    return res
end

"""
    _nuclear_kernel(operator::NuclearOperator)

Compute a nuclear matrix element between two primitive Gaussian basis functions.
"""
function _nuclear_kernel(operator::NuclearOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)
    gamma = alpha + beta
    epsilon = 1 / (4gamma)

    # TODO: the three lines below are bottlenecks and might benefit from
    # inplace operations (in order to avoid allocations). This could be
    # accomplished by reusing intermediate vectors by the caller.
    p_coord = _third_center(alpha, a_coord, beta, b_coord, gamma)
    pa = @. p_coord - a_coord
    pb = @. p_coord - b_coord

    K_p = _gaussian_prod_factor(alpha, a_coord, beta, b_coord, gamma)

    B = zero(gamma)
    for z in 1:length(operator)
        c_charge = operator.charges[z]
        c_coord = @view operator.coords[:, z]

        # TODO: the following has an extra unneeded ^2
        x = gamma * dist(p_coord, c_coord)^2

        # TODO: I think this causes an allocation.
        pc = @. p_coord - c_coord

        c_term = zero(B)
        for l in 0:(la + lb), r in 0:fld(l, 2), i in 0:fld(l - 2r, 2)
            vlri = _vlri(l, r, i, la, lb, pa[1], pb[1], pc[1], epsilon)

            for m in 0:(ma + mb), s in 0:fld(m, 2), j in 0:fld(m - 2s, 2)
                vmsj = _vlri(m, s, j, ma, mb, pa[2], pb[2], pc[2], epsilon)

                for n in 0:(na + nb), t in 0:fld(n, 2), k in 0:fld(n - 2t, 2)
                    vntk = _vlri(n, t, k, na, nb, pa[3], pb[3], pc[3], epsilon)

                    f = boys(l + m + n - 2 * (r + s + t) - (i + j + k), x)
                    c_term += vlri * vmsj * vntk * f
                end
            end
        end
        B += c_charge * c_term
    end

    return -2π * K_p * B / gamma
end
