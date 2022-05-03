"""
    integral_kernel!((p_coord, pa, pb), ::OverlapOperator, a, b, i, j)

Compute the overlap between two primitive Gaussian basis functions.

This overrides `p_coord`, `pa` and `pb`.
"""
# TODO: should we compare inputs and return one if they are the same?
@inline function integral_kernel!((p_coord, pa, pb), ::OverlapOperator, a, b, i, j)
    α = a.alphas[i]
    β = b.alphas[j]

    # It is possible to do a screening based on Kₚ
    Kₚ, γ = gpt!(p_coord, α, a.coord, β, b.coord)
    # if Kₚ < 1e-6
    #     println("Screened $α, $β with $Kₚ")
    #     return zero(Kₚ)
    # end

    @. pa = p_coord - a.coord
    @. pb = p_coord - b.coord

    S₁ = Sᵢ(a.l, b.l, pa[1], pb[1], γ)
    S₂ = Sᵢ(a.m, b.m, pa[2], pb[2], γ)
    S₃ = Sᵢ(a.n, b.n, pa[3], pb[3], γ)

    return Kₚ * S₁ * S₂ * S₃
end

"""
    integral_kernel!(scratch, ::KineticOperator, a, b, i, j)

Compute the kinetic energy matrix element between two primitive Gaussian basis functions.

This overrides `scratch`.
"""
function integral_kernel!(scratch, ::KineticOperator, a, b, i, j)
    Ŝ = OverlapOperator()

    I₁ = (2 * (b.l + b.m + b.n) + 3) * integral_kernel!(scratch, Ŝ, a, b, i, j)
    I₂ = integral_kernel!(scratch, Ŝ, a, adjust_momenta(b, 2, 0, 0), i, j) +
        integral_kernel!(scratch, Ŝ, a, adjust_momenta(b, 0, 2, 0), i, j) +
        integral_kernel!(scratch, Ŝ, a, adjust_momenta(b, 0, 0, 2), i, j)
    I₃ = b.l * (b.l - 1) * integral_kernel!(scratch, Ŝ, a, adjust_momenta(b, -2, 0, 0), i, j) +
        b.m * (b.m - 1) * integral_kernel!(scratch, Ŝ, a, adjust_momenta(b, 0, -2, 0), i, j) +
        b.n * (b.n - 1) * integral_kernel!(scratch, Ŝ, a, adjust_momenta(b, 0, 0, -2), i, j)

    β = b.alphas[j]
    return β * (I₁ - 2β * I₂) - I₃ / 2
end

"""
    integral_kernel!((p_coord, pa, pb, pc), operator::NuclearOperator, a, b, i, j)

Compute a nuclear matrix element between two primitive Gaussian basis functions.

This overrides `p_coord`, `pa`, `pb` and `pc`.
"""
@views function integral_kernel!((p_coord, pa, pb, pc), operator::NuclearOperator, a, b, i, j)
    α = a.alphas[i]
    β = b.alphas[j]

    # It is possible to do a screening based on Kₚ
    Kₚ, γ = gpt!(p_coord, α, a.coord, β, b.coord)
    # if Kₚ < 1e-6
    #     println("Screened $α, $β with $Kₚ")
    #     return zero(Kₚ)
    # end

    @. pa = p_coord - a.coord
    @. pb = p_coord - b.coord

    ϵ = 1 / 4γ
    I = zero(γ)
    for u in 1:length(operator)
        Zᵤ = operator.charges[u]

        u_coord = @view operator.coords[:, u]
        @. pc = p_coord - u_coord

        # TODO: the following has an extra unneeded ^2
        x = γ * dist(p_coord, u_coord)^2
        Iᵤ = zero(I)
        for l in 0:(a.l + b.l), r in 0:fld(l, 2)
            double_max_i = l - 2r
            for i in 0:fld(double_max_i, 2)
                vₗᵣᵢ = v((l, r, i), a.l, b.l, pa[1], pb[1], pc[1], ϵ)

                for m in 0:(a.m + b.m), s in 0:fld(m, 2)
                    double_max_j = m - 2s
                    for j in 0:fld(double_max_j, 2)
                        vₘₛⱼ = v((m, s, j), a.m, b.m, pa[2], pb[2], pc[2], ϵ)

                        for n in 0:(a.n + b.n), t in 0:fld(n, 2)
                            double_max_k = n - 2t
                            for k in 0:fld(double_max_k, 2)
                                vₙₜₖ = v((n, t, k), a.n, b.n, pa[3], pb[3], pc[3], ϵ)

                                n = double_max_i + double_max_j + double_max_k - (i + j + k)
                                Fₙ = boys(n, x)
                                Iᵤ += vₗᵣᵢ * vₘₛⱼ * vₙₜₖ * Fₙ
                            end
                        end
                    end
                end
            end
        end
        I += Zᵤ * Iᵤ
    end

    return -(Kₚ / γ) * I
end

"""
    integral_kernel!((p_coord, q_coord, pq, pa, pb, qc, qd), operator::CoulombOperator, a, b, c, d, i, j, k, l)

Compute a electron repulsion matrix element between two primitive Gaussian basis functions.

This overrides `p_coord`, `q_coord`, `pq`, `pa`, `pb`, `qc` and `qd`.
"""
@inline function integral_kernel!((p_coord, q_coord, pq, pa, pb, qc, qd), operator::CoulombOperator, a, b, c, d, i, j, k, l)
    αa = a.alphas[i]
    αb = b.alphas[j]

    Kp, γp = gpt!(p_coord, αa, a.coord, αb, b.coord)
    @. pa = p_coord - a.coord
    @. pb = p_coord - b.coord

    αc = c.alphas[k]
    αd = d.alphas[l]

    Kq, γq = gpt!(q_coord, αc, c.coord, αd, d.coord)
    @. qc = q_coord - c.coord
    @. qd = q_coord - d.coord

    @. pq = p_coord - q_coord

    δ = 1 / 4γp + 1 / 4γq

    # TODO: precalculate a.l + b.l, etc.

    # TODO: the following has an extra unneeded ^2
    x = dist(p_coord, q_coord)^2 / 4δ
    I = zero(δ)
    for lp in 0:(a.l + b.l), rp in 0:fld(lp, 2), lq in 0:(c.l + d.l), rq in 0:fld(lq, 2)
        double_max_i = lp + lq - 2 * (rp + rq)
        for i in 0:fld(double_max_i, 2)
            gᵢ = g((lp, lq, rp, rq, i), a.l, b.l, c.l, d.l, pa[1], pb[1], qc[1], qd[1], pq[1], γp, γq)

            for mp in 0:(a.m + b.m), sp in 0:fld(mp, 2), mq in 0:(c.m + d.m), sq in 0:fld(mq, 2)
                double_max_j = mp + mq - 2 * (sp + sq)
                for j in 0:fld(double_max_j, 2)
                    gⱼ = g((mp, mq, sp, sq, j), a.l, b.l, c.l, d.l, pa[2], pb[2], qc[2], qd[2], pq[2], γp, γq)

                    for np in 0:(a.n + b.n), tp in 0:fld(np, 2), nq in 0:(c.n + d.n), tq in 0:fld(nq, 2)
                        double_max_k = np + nq - 2 * (tp + tq)
                        for k in 0:fld(double_max_k, 2)
                            gₖ = g((np, nq, tp, tq, k), a.l, b.l, c.l, d.l, pa[3], pb[3], qc[3], qd[3], pq[3], γp, γq)

                            n = double_max_i + double_max_j + double_max_k - (i + j + k)
                            Fₙ = boys(n, x)
                            I += gᵢ * gⱼ * gₖ * Fₙ
                        end
                    end
                end
            end
        end
    end

    return ((Kp / γp) * (Kq / γq) / sqrt(γp + γq)) * I
end
