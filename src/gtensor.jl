struct GTensor{T<:Real, S<:AbstractMatrix{T}} <: AbstractArray{T, 4}
    data::S  # n(n + 1)/2 × n(n + 1)/2 matrix with upper triangular stored
    n::Int
    uplo::Char

    function GTensor{T, S}(data, n, uplo) where {T<:Real, S<:AbstractMatrix{<:T}}
        Base.require_one_based_indexing(data)
        new{T, S}(data, n, uplo)
    end
end
function GTensor(A::AbstractMatrix{T}, n::Int, uplo::Symbol=:U) where {T<:Real}
    checksquare(A)
    return GTensor{T, typeof(A)}(A, n, char_uplo(uplo))
end
# TODO: define GTensor(::AbstractMatrix{T}, uplo::Symbol=:U) => needs to infer n! (same as above with implicit)
# TODO: define GTensor(::AbstractArray{T, 4}, uplo::Symbol=:U) (compacts a larger structure)

function GTensor{T}(::UndefInitializer, n::Int, uplo::Symbol=:U) where {T<:Real}
	d = fld(n * (n + 1), 2)
	A = Matrix{T}(undef, d, d)
	return GTensor(A, n, uplo)
end

Base.size(G::GTensor) = (G.n, G.n, G.n, G.n)

@inline function Base.getindex(G::GTensor, I::Vararg{Int, 4})
    i = _compact_indices(G.uplo, G.n, I[1], I[2])
    j = _compact_indices(G.uplo, G.n, I[3], I[4])

    @boundscheck checkbounds(G.data, i, j)
    @inbounds if i == j
        return hermitian(G.data[i, j], sym_uplo(G.uplo))::hermitian_type(eltype(G.data))
    elseif (G.uplo == 'U') == (i < j)
    return G.data[i, j]
    else
        return adjoint(G.data[j, i])
    end
end

function Base.setindex!(G::GTensor, v, I::Vararg{Int, 4})
    i = _compact_indices(G.uplo, G.n, I[1], I[2])
    j = _compact_indices(G.uplo, G.n, I[3], I[4])

    @boundscheck checkbounds(G.data, i, j)
    @inbounds if (G.uplo == 'U') == (i < j)
        Base.setindex!(G.data, v, i, j)
    else
        Base.setindex!(G.data, v, j, i)
    end
end

# Enumerate a triangular part of a square matrix.
# TODO: should we differentiate between 'U' and 'L'?
# NOTE: other possibilities are: diagonal enumeration and row enumeration
# (which is essentially 'L' column enumeration for 'U' enumeration and vice-versa)
function _compact_indices(uplo::Char, n::Integer, i::Integer, j::Integer)
    if uplo == 'U'
        if i > j
            j, i = i, j
        end

        # Linearization of the upper part of a square matrix
        return _compact_indices(i, j)
    else
        if i < j
            j, i = i, j
        end

        # Linearization of the lower part of a square matrix
        return _compact_indices(i, j) + j * (n - (j - 1)) - n
    end
end
_compact_indices(i::Integer, j::Integer) = i + fld((j - 1) * j, 2)
