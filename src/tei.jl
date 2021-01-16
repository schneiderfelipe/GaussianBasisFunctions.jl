# TODO: create an inplace version that receives the reference of an allocated GTensor.
"""
    tei(
        basis::AbstractVector{<:AbstractBasisFunction},
        operator::AbstractTwoElectronOperator
    )

Compute a symmetric electron repulsion tensor of rank four for a basis set.
"""
function tei(
        basis::AbstractVector{<:AbstractBasisFunction},
        operator::AbstractTwoElectronOperator
    )
    n = length(basis)

    T = typeof(basis[1].coeffs[1])
    G = GTensor{T}(undef, n)

    # The following loop traverses the following order, because this is
    # column-wise in a upper triangular matrix (3 × 3 case means a
    # 2 × 2 × 2 × 2 tensor):
    #
    # (1, 1) => ((1, 1), (1, 1))     => (1, 1, 1, 1)
    # (1, 2) => ((1, 1), (1, 2))     => (1, 1, 1, 2)
    # (2, 2) => ((1, 2), (1, 2))     => (1, 2, 1, 2)
    # (1, 3) => ((1, 1), (2, 2))     => (1, 1, 2, 2)
    # (2, 3) => ((1, 2), (2, 2))     => (1, 2, 2, 2)
    # (3, 3) => ((2, 2), (2, 2))     => (2, 2, 2, 2)
    #
    # This is the memory order, and arguably the most efficient.

    # This loops over an upper triangular part of a matrix of size n(n + 1)/2 × n(n + 1)/2
    scratch = create_scratch(operator)
    constant = create_constant(operator)
    for kl in 1:fld(n * (n + 1), 2)
        for ij in 1:kl
            # TODO: make a function for this triangular transformation.
            #
            # Calculate the side of the triangle we are now filling.
            l = ceil(Int, _trrt(kl))
            # Our current number plus the height of the current column minus
            # the area of triangle above is equal to our row.
            k = kl + l - l * (l + 1) ÷ 2

            j = ceil(Int, _trrt(ij))
            i = ij + j - j * (j + 1) ÷ 2

            # TODO: In this case we know ij and kl, so we could also index
            # G[ij, kl] and avoid inner calculations!
            G[i, j, k, l] = constant * tei!(scratch, basis[i], basis[j], basis[k], basis[l], operator)
        end
    end

    return G
end

"""
    tei(
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        c::GaussianBasisFunction,
        d::GaussianBasisFunction,
        operator::AbstractTwoElectronOperator
    )

Compute an electron repulsion matrix element for a group of basis functions.

Please prefer the potentially more efficient `tei!` over this function.
"""
function tei(
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        c::GaussianBasisFunction,
        d::GaussianBasisFunction,
        operator::AbstractTwoElectronOperator
    )
    return create_constant(operator) * tei!(
        create_scratch(operator), a, b, c, d, operator
    )
end

"""
    tei!(
        scratch,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        c::GaussianBasisFunction,
        d::GaussianBasisFunction,
        operator::AbstractTwoElectronOperator
    )

Compute an electron repulsion matrix element for a group of basis functions.

This function works in-place and overrides `scratch`.
"""
function tei!(
        scratch,
        a::GaussianBasisFunction,
        b::GaussianBasisFunction,
        c::GaussianBasisFunction,
        d::GaussianBasisFunction,
        operator::AbstractTwoElectronOperator
    )
    return combine_gaussians(a, b, c, d) do i, j, k, l
        return integral_kernel!(scratch, operator, a, b, c, d, i, j, k, l)
    end
end
