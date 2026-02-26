using AbstractAlgebra
using DataStructures

include("TypeAKostantPartition.jl")

# This is an implementation of the algorithms described in the article [1]: Canonical bases via pairing monomials.

"""
    AHeckGL

Stores information about the geometry and representation theory of the affine Hecke algebra of GL with a given weight.

Fields:
- `weight::Vector{Int}`: Stores the weight for which the remaining fields are computed
- `kostParts::Vector{TypeAKostantPartition}`: Stores all Kostant partitions of the given weight. This is an indexing set for the irreducibles.
- `aVectors:: Dict{TypeAKostantPartition,Vector{Int}}`: Stores for each Kostant partition the corresponding a-vector.
- `coordSubspaces:: Dict{TypeAKostantPartition,CoordSubspace}`: Stores for each Kostant partition the corresponding corrdinate subpace R_{i, a}.
- `poincarePolys:: Dict{TypeAKostantPartition, AbstractAlgebra.Generic.RationalFunctionFieldElem}`: Stores the quantum factorial of the a-vector.
- `psi::AbstractAlgebra.Generic.MatSpaceElem`: Stores the psi-matrix of the given weight. This encodes the bilinear form on the monomial basis.
- `q::AbstractAlgebra.Generic.MatSpaceElem`: Stores the q-matrix of the given weight. This encodes the change of basis between the canonical and the monomial basis.
- `p::AbstractAlgebra.Generic.MatSpaceElem`: Stores the p-matrix of the given weight. This encodes the change of basis between the standard and the canonical basis or equivalently the stalks of perverse sheaves.
- `h::AbstractAlgebra.Generic.MatSpaceElem`: Stores the h-vector of the given weight. This encodes the bilinear form between the (s,q)-Springer sheaf and the monomial basis.
- `f::AbstractAlgebra.Generic.MatSpaceElem`: Stores the f-vector of the given weight which is defined as f= q^T * psi^(-1) * h.
- `dimOfIrreducibles::Vector{Int}`: Stores the dimensions of the irreducible representations of the affine Hecke algebra with the given weight.
"""
struct AHeckGL
    weight::Vector{Int}
    kostParts::Vector{TypeAKostantPartition}
    aVectors:: Dict{TypeAKostantPartition,Vector{Int}}
    coordSubspaces:: Dict{TypeAKostantPartition,CoordSubspace}
    poincarePolys:: Dict{TypeAKostantPartition, AbstractAlgebra.Generic.RationalFunctionFieldElem}
    psi::AbstractAlgebra.Generic.MatSpaceElem
    q::AbstractAlgebra.Generic.MatSpaceElem
    p::AbstractAlgebra.Generic.MatSpaceElem
    h::AbstractAlgebra.Generic.MatSpaceElem
    f::AbstractAlgebra.Generic.MatSpaceElem
    dimOfIrreducibles::Vector{Int}

    function AHeckGL(w::Vector{Int})
        weight = w
        kostParts = kostantPartitions(w)
        aVectors = Dict(c => aVector(c) for c in kostParts)
        coordSubspaces = Dict(c => CoordSubspace(weight,aVectors[c]) for c in kostParts)
        poincarePolys = Dict(c => poincarePoly(aVectors[c]) for c in kostParts)
        psi = psiMatrix(weight, kostParts, coordSubspaces, poincarePolys)
        l = ldltDecomposition(psi)[1]
        q,p = qpDecomposition(l)
        h=hVector(weight, kostParts, coordSubspaces, poincarePolys)
        f = transpose(q) * psi^(-1) * h
        new(weight, kostParts, aVectors, coordSubspaces, poincarePolys, psi, q, p, h, f, evalAtOne(f))
    end
end

"""
    psiMatrix(weight::Vector{Int},
                    kostParts::Vector{TypeAKostantPartition},
                    coordSubspaces:: Dict{TypeAKostantPartition,CoordSubspace},
                    poincarePolys)

Computes the psi matrix assoicated to the given `weight`.
This differs from the definition of the psi-matrix in [1] by a factor of (1-x^2)^n.
"""
function psiMatrix(weight::Vector{Int},
                    kostParts::Vector{TypeAKostantPartition},
                    coordSubspaces:: Dict{TypeAKostantPartition,CoordSubspace},
                    poincarePolys)::AbstractAlgebra.Generic.MatSpaceElem
    R, x = rational_function_field(QQ, "x")
    coxlens = coxeterLengths(weight)
    allconjugates = Dict(c => allConjugates(coordSubspaces[c], coxlens) for c in kostParts)
    psi = identity_matrix(R,length(kostParts))
    for i in eachindex(kostParts)
        for j in i:length(kostParts) 
            c1 = kostParts[i]
            c2 = kostParts[j]
            entry = innerProduct(allconjugates[c1], poincarePolys[c1], coordSubspaces[c2], poincarePolys[c2])
            psi[i,j] = entry
            psi[j,i] = entry
        end
    end
    return psi
end

"""
    ldltDecomposition(mat::AbstractAlgebra.Generic.MatSpaceElem)


Computes the LDL^T-decomposition of the matrix `mat` where L is lower unitriangular und D is diagonal.
"""
function ldltDecomposition(mat::AbstractAlgebra.Generic.MatSpaceElem)::Tuple{AbstractAlgebra.Generic.MatSpaceElem,AbstractAlgebra.Generic.MatSpaceElem}
    R, x = rational_function_field(QQ, "x")
    n = nrows(mat)
    D = copy(mat)
    L = identity_matrix(R,n)
    for i in 1:n
        for j in i+1:n
            lambda = D[i,j]/D[i,i]
            add_column!(D, -lambda, i, j)
            add_row!(D, -lambda, i, j)
            add_column!(L, lambda, j, i)
        end
    end
    return(L,D)
end

"""
    nonNegativePart(f::AbstractAlgebra.Generic.RationalFunctionFieldElem)

Returns the non-negative part of rational function `f` that is a Laurent polynomial.
"""
function nonNegativePart(f::AbstractAlgebra.Generic.RationalFunctionFieldElem)::AbstractAlgebra.Generic.Poly
    if !is_monomial(denominator(f))
        error("$f is not a Laurent polynomial")
    end
    return divrem(numerator(f), denominator(f))[1]
end

"""
    symmetrize(f::AbstractAlgebra.Generic.Poly)

Symmetrizes the polynomial `f` to a Laurent polynomial that is
invariant under x -> x^(-1) and has the same non-negative part as `f`.
"""
function symmetrize(f::AbstractAlgebra.Generic.Poly)::AbstractAlgebra.Generic.RationalFunctionFieldElem
    R, x = rational_function_field(QQ, "x")
    if f == 0
        return f/x^0
    else
        return f/x^0 + reverse(divrem(f, polynomial_ring(QQ,"x")[2])[1])/x^degree(f)
    end
end

"""
    qpDecomposition(L::AbstractAlgebra.Generic.MatSpaceElem)

Computes the QP-decomposition of a lower triangular matrix `L` where P and Q are lower unitriangular,
Q has entries in Z[x+x^{-1}] and P has entries in x^{-1} Z[x^{-1}] below the diagonal.
"""
function qpDecomposition(L::AbstractAlgebra.Generic.MatSpaceElem)::Tuple{AbstractAlgebra.Generic.MatSpaceElem,AbstractAlgebra.Generic.MatSpaceElem}
    R, x = rational_function_field(QQ, "x")
    n = nrows(L)
    P = copy(L)
    Q = identity_matrix(R,n)
    for i in 1:n
        for j in i-1:-1:1
            lambda = symmetrize(nonNegativePart(P[i,j]))
            add_row!(P, -lambda, j, i)
            add_column!(Q, lambda, i, j)
        end
    end
    return(Q,P)
end

"""
innerProduct(V1conjugates:: Vector{Tuple{CoordSubspace, Int}}, a1::Vector{Int}, p1::AbstractAlgebra.Generic.RationalFunctionFieldElem, V2::coordSubspace, a2::Vector{Int}, p2::AbstractAlgebra.Generic.RationalFunctionFieldElem)

Computes the inner product between `V1` and `V2`.
The function requires all conjugates of `V1` for this as well as the Poincare polynomials associated to the a-vectors of `V_1` and `V_2`.
"""
function innerProduct(V1conjugates:: Vector{Tuple{CoordSubspace, Int}},
                    p1::AbstractAlgebra.Generic.RationalFunctionFieldElem,
                    V2::CoordSubspace,
                    p2::AbstractAlgebra.Generic.RationalFunctionFieldElem)::AbstractAlgebra.Generic.RationalFunctionFieldElem
    R, x = rational_function_field(QQ, "x")
    sigma = 0*x
    exponents = DefaultDict{Int, Int}(0)
    for con in V1conjugates
        exponents[2*(dimOfIntersection(con[1], V2) + con[2])] += 1
    end
    for (e, k) in exponents
        sigma += k*x^(e)
    end
    scalar = x^(-V1conjugates[1][1].dim - V2.dim)/(p1*p2)
    return scalar*sigma
end


"""
    hProduct(nspaces:: Vector{Tuple{CoordSubspace, Int}}, V::CoordSubspace, a::Vector{Int})

Computes the entry of the h-vector associated to `V`.
This differs from the definition of the h-vector in [1] by a factor of (1-x^2)^n.
"""
function hProduct(nspaces:: Vector{Tuple{CoordSubspace, Int}},
                V::CoordSubspace,
                p::AbstractAlgebra.Generic.RationalFunctionFieldElem)::AbstractAlgebra.Generic.RationalFunctionFieldElem
    R, x = rational_function_field(QQ, "x")
    sigma = 0*x
    exponents = DefaultDict{Int, Int}(0)
    for space in nspaces
        exponents[2*(dimOfIntersection(space[1], V) + space[2])] += 1
    end
    for (e, k) in exponents
        sigma += k*x^(e)
    end
    scalar = x^(- V.dim)/p
    return scalar*sigma
end

"""
    hVector(weight:: Vector{Int},
            kostParts::Vector{TypeAKostantPartition},
            aVectors:: Dict{TypeAKostantPartition,Vector{Int}},
            coordSubspaces:: Dict{TypeAKostantPartition,CoordSubspace},
            poincarePolys:: Dict{TypeAKostantPartition, AbstractAlgebra.Generic.RationalFunctionFieldElem})

Computes the h-vector of a given weight.
This differs from the definition of the h-vector in [1] by a factor of (1-x^2)^n.
"""
function hVector(weight:: Vector{Int},
                kostParts::Vector{TypeAKostantPartition},
                coordSubspaces:: Dict{TypeAKostantPartition,CoordSubspace},
                poincarePolys)::AbstractAlgebra.Generic.MatSpaceElem
    R, x = rational_function_field(QQ, "x")
    nsubspaces = allnSubspaces(weight)
    H = zero_matrix(R, length(kostParts), 1)
    for i in eachindex(kostParts)
        c = kostParts[i]
        entry = hProduct(nsubspaces, coordSubspaces[c], poincarePolys[c])
        H[i,1] = entry
    end
    return H
end

"""
    poincarePoly(n::Int)

Computes the normalized Poincare polynomial of `n`
"""
function poincarePoly(n::Int)::AbstractAlgebra.Generic.RationalFunctionFieldElem
    R, x = rational_function_field(QQ, "x")
    temp = 1*x^0
    res = 1*x^0
    for i in 1:n
        res = res*temp
        temp = temp + x^(2*i)
    end
    return x^(-n*(n-1)÷2)*res
end

"""
    poincarePoly(a::Vector{Int})

Computes the product of the Poincare polynomials of the entries of `a`.
"""
function poincarePoly(a::Vector{Int})::AbstractAlgebra.Generic.RationalFunctionFieldElem
    R, x = rational_function_field(QQ, "x")
    res = 1*x^0
    for i in 1:length(a)
        res = res*poincarePoly(a[i])
    end
    return res
end

"""
    evalAtOne(f::AbstractAlgebra.Generic.MatSpaceElem)

Computes for each entry in a matrix of rational functions the evaluation of the numerator at 1.
"""
function evalAtOne(f::AbstractAlgebra.Generic.MatSpaceElem)::Vector{Int}
    ans = []
    for p in f
        push!(ans, Int(evaluate(numerator(p), 1)))
    end
    return ans
end