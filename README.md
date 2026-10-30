# DAC-SUB-DELAUNAY: Divide and Conquer Subspace Delaunay

C++17 implementation of the construction described in *Constructing ε-Constrained
Sparsified βˢ-Complexes using Space Partitioning Trees* (R. P. Singh, P. A. Wilsey).

| Paper | Code |
|---|---|
| Algorithm 1: SP-tree via spherical k-means | `buildSPTree`, `KWayTreeNode` |
| Algorithm 2: incremental mesh generation (post-order) | `Mesher::rec` |
| Algorithm 3: merged mesh, sort by circum-radius, validate | `Mesher::finish` |
| Algorithm 4: cross-simplex enumeration in the ε-neighbourhood | `Mesher::crossSimplices` |
| ε-DIST neighbourhood search | `ballQuery` |
| β-criterion (β ≥ 1 / β < 1 exclusion regions) | `exclusionFactor`, `Mesher::valid` |

For each simplex dimension d = 1..s the program builds all d-simplices (d+1 vertices)
with all edges ≤ ε whose β-exclusion ball (centred at the circum-centre inside the
simplex's own affine hull) contains no data point. β = 0 gives the Vietoris–Rips
simplices, β = 1 the Alpha-like (empty circumsphere) simplices.

## Build and run

    make
    ./dac_sub_delaunay --gen sphere --dim 4 --n 500 --s 3 --eps 0.5 --beta 1 --verify
    ./dac_sub_delaunay --input points.csv --s 2 --eps 0.75 --beta 1 --out complex.txt

`--verify` recomputes the mesh by brute force and checks that both agree. `make test` runs
three such checks. Run `./dac_sub_delaunay --help` for all options.

Output lines (`--out`): `dim weight v0 v1 ...` where weight is the circum-radius.
A small H0 summary (union-find on the 1-simplices) is printed; for full persistent
homology feed the filtered complex to a PH backend such as LHF.

## Notes and deviations from the paper

* Cross simplices are enumerated among "exposed" points (points with an ε-neighbour in a
  different child), which is exact for ε-bounded complexes. The paper's hyperplane-based
  β-exposure test and the PRUNEINVALIDSIMPLICES step are not implemented; validation uses
  full SP-tree range queries, so results are identical to brute force (verified).
* For β > 1 the same exclusion radius as β = 1 is used, as in the paper's formulas.
