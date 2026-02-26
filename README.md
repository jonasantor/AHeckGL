# Representations of the Affine Hecke algebra of $GL_n$

This repository contains a Julia-implementation of the algorithms developed in the paper "Canonical bases via pairing monomials" henceforth cited as [1]. 

## Functionality
The code in this repository defines a Julia-struct **AHeckGL** which can store (and compute) the dimensions of the irreducible representations of the affine Hecke algebra of $GL_n$. In the process the canonical basis of a corresponding quantum group of type A is also computed. Given a weight (= central character for the affine Hecke algebra), the struct AHeckGL bundles the following data:

- **kostparts**: Stores all Kostant partitions of the given weight. This is an indexing set for the irreducible representations of the affine Hecke algebra and the corresponding canonical basis of a quantum group with the given weight/central character.
- **dimOfIrreducible**: Stores the dimension of all the irreducible representation with the given weight. This is a vector whose order is according to the order of the Kostant partitions in `kostparts`.
- **p**: Stores the p-matrix from [1] whose entries are Laurent polynomials. This encodes the change of basis between the standard and the canonical basis in a quantum group. Alternatively, this encodes the stalks of the corresponding simple perverse sheaves.
- **q**: Stores the q-matrix from [1] which encodes the change of basis between the canonical and the monomial basis.
- **psi**: Stores the psi-matrix from [1] which encodes the bilinear form on the monomial basis. Note that the definition of the psi matrix in [1] differs by a factor of $(1-x^2)^-n$ from the implementation in this repository.


Note that the canonical basis of a quantum group can be computed faster using the QuaGroups package, but this does not compute the dimensions of irreducible representations for the affine Hecke algebra, change of basis to the monomial basis or bilinear form.

The file *Example.ipynb* contains a demonstration of the functionality of the struct AHeckGL. The file *AHeckGL.jl* defines the struct AHeckGL and implements the algorithms from [1]. This builds on the functionality in *TypeAKostantPartitions* which implements and computes Kostant partitions of type A using the language of triangular arrays from "Combinatorics of Fourier transforms for type A quiver representations" by Pramod N. Achar, Maitreyee C. Kulkarni and Jacob P. Matherne.