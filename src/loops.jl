raw"""
    combine_gaussians(
        f::Function,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction
    )

Linearly combine the primitives of two Gaussian basis functions.

This function returns

    ``\Sum_i \Sum_j d_i d_j N_i N_j f(i, j)``

where ``N_i`` is a normalization constant, and `f` is a closure that receives
the indices of each primitive pair.
"""
@inline function combine_gaussians(
        f::Function,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction
    )
    I = zero(a.coeffs[1])
    for i in 1:length(a)
        Nᵢ = gaussian_norm(a.alphas[i], a.l, a.m, a.n)

        Iᵢ = zero(I)
        for j in 1:length(b)
            Nⱼ = gaussian_norm(b.alphas[j], b.l, b.m, b.n)

            Iᵢ += b.coeffs[j] * Nⱼ * f(i, j)
        end

        I += a.coeffs[i] * Nᵢ * Iᵢ
    end

    return I
end

raw"""
    combine_gaussians(
        f::Function,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction
    )

Linearly combine the primitives of four Gaussian basis functions.

This function returns

    ``\Sum_i \Sum_j \Sum_k \Sum_l d_i d_j d_k d_l N_i N_j N_k N_l f(i, j, k, l)``

where ``N_i`` is a normalization constant, and `f` is a closure that receives
the indices of each primitive pair.
"""
@inline function combine_gaussians(f::Function, a::GaussianBasisFunction, b::GaussianBasisFunction, c::GaussianBasisFunction, d::GaussianBasisFunction)
    return combine_gaussians(a, b) do i, j
        return combine_gaussians(c, d) do k, l
            return f(i, j, k, l)
        end
    end
end
