"""
    integral_kernel(::OverlapOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)

Compute the overlap between two primitive Gaussian basis functions.
"""
# TODO: should we compare inputs and return one if they are the same?
@inline function integral_kernel(::OverlapOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)
    gamma = alpha + beta

    # It is possible to do a screening based on K_p
    K_p = _gaussian_prod_factor(alpha, a_coord, beta, b_coord, gamma)
    # if K_p < 1e-6
    #     println("Screened $alpha, $beta with $K_p")
    #     return zero(K_p)
    # end

    # TODO: the three lines below are bottlenecks and might benefit from
    # inplace operations (in order to avoid allocations). This could be
    # accomplished by reusing intermediate vectors by the caller.
    p_coord = _third_center(alpha, a_coord, beta, b_coord, gamma)
    pa = @. p_coord - a_coord
    pb = @. p_coord - b_coord

    S_x = _Si(la, lb, pa[1], pb[1], gamma)
    S_y = _Si(ma, mb, pa[2], pb[2], gamma)
    S_z = _Si(na, nb, pa[3], pb[3], gamma)

    return K_p * S_x * S_y * S_z
end

"""
    integral_kernel(::KineticOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)

Compute the kinetic energy matrix element between two primitive Gaussian basis functions.
"""
function integral_kernel(::KineticOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)
    overlap_op = OverlapOperator()

    l0 = (2 * (lb + mb + nb) + 3) * integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb
    )
    lplus2 = integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb + 2, mb, nb
    ) + integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb + 2, nb
    ) + integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb + 2
    )
    lminus2 = lb * (lb - one(lb)) * integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb - 2, mb, nb
    ) + mb * (mb - one(mb)) * integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb - 2, nb
    ) + nb * (nb - one(nb)) * integral_kernel(
        overlap_op, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb - 2
    )

    return beta * (l0 - 2 * beta * lplus2) - lminus2 / 2
end

"""
    integral_kernel(operator::NuclearOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)

Compute a nuclear matrix element between two primitive Gaussian basis functions.
"""
function integral_kernel(operator::NuclearOperator, alpha, a_coord, la, ma, na, beta, b_coord, lb, mb, nb)
    gamma = alpha + beta

    # It is possible to do a screening based on K_p
    K_p = _gaussian_prod_factor(alpha, a_coord, beta, b_coord, gamma)
    # if K_p < 1e-5
    #     println("Screened $alpha, $beta with $K_p")
    #     return zero(K_p)
    # end

    # TODO: the three lines below are bottlenecks and might benefit from
    # inplace operations (in order to avoid allocations). This could be
    # accomplished by reusing intermediate vectors by the caller.
    p_coord = _third_center(alpha, a_coord, beta, b_coord, gamma)
    pa = @. p_coord - a_coord
    pb = @. p_coord - b_coord  # TODO: can we use @view on this?

    epsilon = 1 / (4gamma)
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

    return -K_p * B * 2π / gamma
end
