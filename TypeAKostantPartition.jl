using Permutations

""""
    TriangularArray

Stores triangular arrays in the sense of Achar-Kulkarni-Matherne.
"""
struct TriangularArray
    size::Int
    weight::Vector{Int}
    entries::Dict{Tuple{Int,Int}, Int}
end

"""
    zeroTriangularArray(n::Int)

Creates the zero triangular array of size `n`.
"""
function zeroTriangularArray(n::Int)::TriangularArray
    entries = Dict((i,j) => 0 for i in 1:n, j in 1:n if j <= n-i+1)
    return TriangularArray(n, zeros(Int, n), entries)
end

"""
    addOneToEntry(i::Int, j::Int, T::TriangularArray)

Creates the triangular array with the same entries as `T`,
except that we add 1 to the (`i`, `j`)-th entry.
"""
function addOneToEntry(i::Int, j::Int, T::TriangularArray)::TriangularArray
    n = T.size
    if i < 1 || i  > n || j < 1 || j > n-i+1 
        error("Index ($i,$j) not in range for size $n triangular array") 
    end
    entries = copy(T.entries)
    entries[(i,j)] = entries[(i,j)] +1
    w = copy(T.weight)
    w[i] += 1
    return TriangularArray(n, w, entries)
end

"""
    addZeroChute(T::TriangularArray)

Creates a triangular array which is `T`
with a zero chute added on top.
"""
function addZeroChute(T::TriangularArray)::TriangularArray
    n = T.size
    new_entries = Dict((i+1,j) => T.entries[(i,j)] for i in 1:n, j in 1:n if j <= n-i+1)
    new_entries = merge(new_entries, Dict((1,j) => 0 for j in 1:n+1))
    return TriangularArray(n+1, pushfirst!(T.weight,0), new_entries)
end

"""
    notContained(list::Vector{TriangularArray}, T:TriangularArray)

Checks if `T` is not contained in `list`.
"""
function notContained(list::Vector{TriangularArray}, T::TriangularArray)::Bool
    for T2 in list
        if T.entries == T2.entries
            return false
        end
    end
    return true
end

"""
    triangularArrays(w::Vector{Int})

Computes all triangular arrays of weight `w`.
The output is a list that is ordered increasingly
with respect to the chute-wise dominance order.
"""
function triangularArrays(w::Vector{Int})::Vector{TriangularArray}
    n = length(w)
    if n == 0
        # If n == 0 there is only the zero triangular array.
        return [zeroTriangularArray(0)]
    elseif w[1] == 0
        # If w[1] == 0 triangular arrays of weight w are obtained from those of weight w[2:n] by adding a zero chute on top.
        return [addZeroChute(T) for T in triangularArrays(w[2:n])]
    else
        # Computes all triangular with weight w from those with weight w2 = w - (1,0,...,0).
        w2 = copy(w)
        w2[1] -= 1
        result = Vector{TriangularArray}()
        for T in triangularArrays(w2)
            # If T has weight w2, adding 1 to the (1,1)-entry is always a valid triangular array of weight w.
            T2 = addOneToEntry(1,1,T)
            if notContained(result, T2)
                push!(result, T2)
            end
            for i in 2:n
                if T.entries[(1,i)] < T.entries[(2,i-1)]
                    # If T has weight w2 and i > 1, adding 1 to the (1,i)-th entry of T
                    # is a valid triangular array of weight w iff T.entries[(1,i)] < T.entries[(2,i-1)].
                    T2 = addOneToEntry(1,i,T)
                    if notContained(result, T2)
                        push!(result,T2)
                    end
                end
            end
        end
        return result
    end
end

"""
    TypeAKostantPartition

Stores Kostant partitions for the type A root systems of a given `rank`
"""
struct TypeAKostantPartition
    rank::Int
    weight::Vector{Int}
    values::Dict{Tuple{Int,Int}, Int}
end

function Base.show(io::IO, ::MIME"text/plain", kp::TypeAKostantPartition)
    println(io, "Type A Kostant Partition")
    println(io, "  Weight : ", kp.weight)
    println(io, "  Values :")
    
    for ((i, j), v) in sort(collect(kp.values))
        println(io, "    alpha(", i, ",", j, ") ↦ ", v)
    end
end

"""
    KostantPartition(T::TriangularArray)

Computes the Kostant partition associated to the triangular array `T`.
"""
function kostantPartition(T::TriangularArray)::TypeAKostantPartition
    rank = T.size
    values = Dict{Tuple{Int,Int}, Int}()
    for j in 1:rank
        values[(1,j)] = T.entries[(1,j)]
    end
    for i in 2:rank
        for j in i:rank
            values[(i,j)] = T.entries[(i,j-i+1)] - T.entries[(i-1,j-i+2)]
        end
    end
    return TypeAKostantPartition(rank, T.weight, values)
end

"""
    kostantPartitions(w::Vector{Int})

Computes all type A Kostant partitions of weight `w`.
The output is a list that is ordered increasingly
with respect to the orbit closure order.
"""
function kostantPartitions(w::Vector{Int})::Vector{TypeAKostantPartition}
    return [kostantPartition(T) for T in triangularArrays(w)]
end

"""
    aVector(c::TypeAKostantPartition)

Computes the a-vector associated to the Kostant partition `c`.
"""
function aVector(c::TypeAKostantPartition)::Vector{Int}
    rank = c.rank
    a = Vector{Int}()
    for i in rank:-1:1
        for j in i:rank
            push!(a, sum([c.values[(i,k)] for k in j:rank]))
        end
    end
    return a
end



"""
    CoordSubspace

Stores coordinate subspaces inside a space of equioriented type A quiver representations,
i.e. coordinate subspaces inside the space of matrices M_{w2,w1} x M_{w3,w2} x ... x M_{wn, wn-1}.
"""
struct CoordSubspace
    n::Int
    weight::Vector{Int}
    dim::Int
    coordinates::Vector{Set{Tuple{Int,Int}}}
end



"""
    dimOfIntersection(V1::CoordSubspace, V2::CoordSubspace)

Computes the dimension of the intersection of `V1` and `V2`.
"""
function dimOfIntersection(V1::CoordSubspace, V2::CoordSubspace)::Int
    return sum([length(intersect(V1.coordinates[i], V2.coordinates[i])) for i in 1:V1.n-1])
end

"""
    permuteCoordinates(p1::Permutation, p2::Permutation, mat::Set{Tuple{Int,Int}})

Computes the product `p1` * `mat` * (`p2`)^-1.
"""
function permuteCoordinates(p1::Permutation, p2::Permutation, mat::Set{Tuple{Int,Int}})::Set{Tuple{Int,Int}}
    newmat = Set{Tuple{Int,Int}}()
    for t in mat
        push!(newmat, (p1(t[1]), p2(t[2])) )
    end
    return newmat
end

"""
    allPermutationsOfCoordinates(n::Int, m::Int, mat::Set{Tuple{Int,Int}})

Computes the products `p1` * `mat` * (`p2`)^-1 for all possible permutations `p1` and `p2`.
"""
function allPermutationsOfCoordinates(n::Int, m::Int, mat::Set{Tuple{Int,Int}})::Dict{Tuple{Permutation, Permutation}, Set{Tuple{Int,Int}}}
    return Dict((p1,p2) => permuteCoordinates(p1,p2,mat) for p1 in PermGen(n), p2 in PermGen(m))
end

"""
    coxeterLength(p::Permutation)

Computes the Coxeter length of `p`.
"""
function coxeterLength(p::Permutation)::Int
    return length(CoxeterDecomposition(p).terms)
end

"""
    coxeterLengths(weight::Vector{Int})

Computes the Coxeter lengths of all permutations on all sets with `weight`[i] elements for any i.
"""
function coxeterLengths(weight::Vector{Int})::Dict{Permutation, Int}
    coxlens = Dict{Permutation, Int}()
    for i in 1:length(weight)
        for p in PermGen(weight[i])
            coxlens[p] = coxeterLength(p)
        end
    end
    return coxlens
end

"""
    multiPermGen(weight::Vector{Int})

Creates an iterator over all permutations in the Weyl group of a given `weight`,
i.e. the permutations in the product of symmetric groups S_w1 x ... x S_wn.
"""
function multiPermGen(weight::Vector{Int})
    return  Iterators.product((PermGen(weight[i]) for i in 1:length(weight))...)
end

"""
    allConjugates(V::CoordSubspace, coxlens::Dict{Permutation, Int})

Computes all possible conjugates of the given coordinate subspace `V`
together with the Coxeter length of the element that was conjugated by.
"""
function allConjugates(V::CoordSubspace, coxlens::Dict{Permutation, Int})::Vector{Tuple{CoordSubspace, Int}}
    n = V.n
    conjlist = Vector{Tuple{CoordSubspace, Int}}()
    permsvec = [allPermutationsOfCoordinates(V.weight[i+1],V.weight[i], V.coordinates[i]) for i in 1:n-1]
    for perms in multiPermGen(V.weight)
        coordinates = Vector{Set{Tuple{Int,Int}}}()
        for i in 1:n-1
            push!(coordinates, permsvec[i][(perms[i+1], perms[i])])
        end
        newspace = CoordSubspace(V.n, V.weight, V.dim, coordinates)
        coxlen = sum([coxlens[perms[i]] for i in 1:n])
        push!(conjlist, (newspace, coxlen))
    end
    return conjlist
end

"""
    coordSubspace(rank::Int, a::Vector{Int})

Computes the coordinate subspace associated to `a`.
"""
function coordSubspace(weight::Vector{Int}, a::Vector{Int})::CoordSubspace
    n = length(weight)
    coordinates = [Set{Tuple{Int,Int}}() for i in 1:n-1]
    flagdim = zeros(Int, n)
    cycle = 1
    index = n
    dim = 0
    for i in length(a):-1:1
        if index < n
            for rowind in 0:flagdim[index+1]-1, colind in 0:a[i]-1
                push!(coordinates[index], (weight[index+1]-rowind, weight[index]-flagdim[index]-colind) )
            end
            dim += flagdim[index+1]*a[i]
        end
        flagdim[index] += a[i]
        if index > cycle
            index -= 1
        else
            index = n
            cycle +=1
        end
    end
    return CoordSubspace(n, weight, dim, coordinates)
end

"""
    nSubspace(weight::Vector{Int}, p::Permutation)

Computes the subspsace p*n*p^{-1} intersected with the space of type A quiver reps
where n is the positive part of the Lie algebra of GL.
"""
function nSubspace(weight::Vector{Int}, p::Permutation)::CoordSubspace
    n = length(weight)
    coordinates = [Set{Tuple{Int,Int}}() for i in 1:n-1]
    rowind=0
    colind=0
    dim = 0
    for i in 1:n-1
        rowind += weight[i]
        for a in 1:weight[i+1], b in 1:weight[i]
            if p(a+rowind) > p(b+colind)
                dim += 1
                push!(coordinates[i], (a,b))
            end
        end
        colind += weight[i]
    end
    return CoordSubspace(n, weight, dim, coordinates)
end

"""
    allnSubspaces(weight::Vector{Int})(V::CoordSubspace, coxlens::Dict{Permutation, Int})

Computes all nsubspaces pnp^-1 intersected with the space of type A quiver reps
for all permutations p.
"""
function allnSubspaces(weight::Vector{Int})::Vector{Tuple{CoordSubspace, Int}}
    return [(nSubspace(weight, p), coxeterLength(p)) for p in PermGen(sum(weight))]
end