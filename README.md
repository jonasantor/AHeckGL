# Representations of the Affine Hecke Algebra of $GL_n$

This repository contains a Julia-implementation of the algorithms developed in the paper "Canonical bases via pairing monomials" henceforth cited as [1]. 

## Functionality
The code in this repository defines a Julia-struct **AHeckGL** which can store and compute the dimensions of the irreducible representations of the affine Hecke algebra of $GL_n$. In the process the canonical basis of a corresponding quantum group of type A is also computed. Given a weight $w = [w_1 , ..., w_l]$ (= a central character for the affine Hecke algebra), the struct AHeckGL bundles the following data:

- **kostParts**: Stores all type $A$ Kostant partitions of the given weight. This is an indexing set for the irreducible representations of the affine Hecke algebra with the given weight. At the same time, this is an indexing set for the corresponding canonical basis (see [1]).
- **dimOfIrreducibles**: Stores the dimension of all the irreducible representation with the given weight. This is a vector which is sorted according to the order of the Kostant partitions in `kostParts`.
- **p**: Stores the $P$-matrix from [1] whose entries are Laurent polynomials. This encodes the change of basis between the standard and the canonical basis. Alternatively, this encodes the stalks of the corresponding simple perverse sheaves.
- **q**: Stores the $Q$-matrix from [1] whose entries are Laurent polynomials. This encodes the change of basis between the canonical and the monomial basis.
- **psi**: Stores the $\Psi$-matrix from [1] which encodes the bilinear form on the monomial basis. Note that the definition of the psi matrix in [1] differs by a factor of $(1-v^2)^{-n}$ from the implementation in this repository. This does not change the $P$-matrix or the $Q$-matrix.

## Computational limits

The algorithm in [1] is of exponential complexity in $n = w_1 + ... + w_l$, so it is only feasible for small values of n. For practical purposes, the computations should be doable for weights with $n \le 10$, though computations for $n=10$ may take several minutes (up to an hour for $w = [2,2,2,2,2]$).

The canonical basis of a quantum group can be computed faster using the <a href="https://gap-packages.github.io/quagroup/">QuaGroups package</a>. However, this does not compute the dimensions of irreducible representations for the affine Hecke algebra or the change of basis to the monomial basis.

## File structure

The file *Example.ipynb* contains a demonstration of the functionality of the struct AHeckGL. The file *AHeckGL.jl* defines the struct AHeckGL and implements the algorithms from [1]. This builds on the functionality in *TypeAKostantPartitions.jl* which implements and computes Kostant partitions of type A using the language of triangular arrays from "Combinatorics of Fourier transforms for type A quiver representations" by Pramod N. Achar, Maitreyee C. Kulkarni and Jacob P. Matherne.