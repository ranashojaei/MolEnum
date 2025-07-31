# MolEnum
Enumeration of Chemical Isomers Using Generating Functions and Pólya’s Counting Theorem

This repository contains Julia code developed for the paper:

Shojaei, R. & Gross, T.
Counting Chemical Isomers with Multivariate Generating Functions

The code computes and analyzes the number of stereo and structural isomers for various types of organic molecules. 

Overview of Codes
1. AcyclicIsoEnum.jl
Computes the coefficient a for acyclic organic molecules.
Represents the number of rooted, tree-like molecular isomers.
Includes optional visualization to illustrate isomer counts based on atom types.

2. BenzeneIsoEnum.jl
Computes the base coefficient a (as in the acyclic case).
Uses this to compute the number of isomers for benzene-derived molecules.

3. NaphthaleneIsoEnum.jl
Also begins by computing the coefficient a for acyclic structures.
Then calculates the number of isomers for naphthalene-like molecules, i.e., structures with two fused benzene rings.

Indexing Note
In Julia (as in most programming languages), array indexing starts at 1. To reflect this, all coefficients (calculated from the 
enumeration equations in the paper) are stored with an offset of +1 in each index.
For example: The coefficient a[0, 0, 1] is stored at a[1, 1, 1].
Likewise, a coefficient corresponding to (L, M, N) atoms is stored at a[L+1, M+1, N+1].
Therefore, to calculate the number of isomers for a molecule with: 5 carbon, 2 oxygen, and 13 hydrogen atoms, Enter L=6, M=3, N=14 in the codes.

However, you do not need to adjust indices in the output — the printed results already reflect the correct values
for L, M, and N. The offset is handled internally for storage and access.
