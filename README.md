# SPARSEMATRIX: A REDUCE sparse matrix package

**[Francis Wright](https://sites.google.com/site/fjwcentaur)**<br/>
Time-stamp: <2026-05-15 18:21:36 franc>

A [*sparse matrix*](https://en.wikipedia.org/wiki/Sparse_matrix) is a matrix in which most of the elements are zero.  By contrast, if most of the elements are non-zero, the matrix is considered *dense*.  Sparse matrices benefit from being stored using different data structures and manipulated using different algorithms from dense matrices.  Whether it is more efficient to regard a matrix (or more likely a set of matrices) as dense or sparse is ill defined and probably depends on context, so it may be determinable only by experiment, but it is reasonable to assume that in a sparse matrix fewer than half the elements are nonzero.

The REDUCE MATRIX package implicitly assumes dense matrices.  SPARSEMATRIX is a re-implementation of the MATRIX package that assumes sparse matrices.  It uses hash tables as the primary data structures to store matrix elements, and the canonical form for a sparse matrix is a LISP list of the form `(<hash> <m> <n>)`, where `<hash>` is a hash table, `<m>` is the number of rows (row dimension), and `<n>` is the number of columns (column dimension).  Only nonzero elements are ever stored in the has table; missing elements are implicitly zero.  By contrast, the canonical form for a dense matrix is a LISP list of rows of the form `(<row_1> <row_2> ... <row_m>)`, where each `<row_i>` is a list of the matrix elements in the *i-th* row.  The algorithms used to manipulate sparse matrices avoid accessing implicitly-zero (i.e. non-stored) matrix elements as much as possible, whereas the algorithms used to manipulate dense matrices always run through all matrix elements, regardless of their values.

Currently, sparse and dense matrices cannot be mixed, but I plan to implement full and automatic interoperability in future.  However, a sparse matrix can be explicitly converted to a dense matrix, and vice versa.

## Critique of SPARSE, an alternative sparse matrix package

There is already a REDUCE package called SPARSE to support sparse matrices, written by Stephen Scowcroft in 1995 at the Konrad-Zuse-Zentrum für Informationstechnik Berlin.  It uses a different data structure to represent a sparse matrix, namely a LISP list of the form `(<vector> (spm <m> <n>))`, where `<vector>` is a vector of rows of the form `[nil <row_1> <row_2> ... <row_m>]` and each `<row_i>` is a list of the matrix elements in the *i-th* row, each represented as a dotted pair, of the form `((nil) (j_1 . val_1) (j_2 . val_2) ...)`.  The `nil` elements are to allow for the fact that REDUCE indexes for vectors and lists start from 0, whereas matrix indices start from 1.  Note that the SPARSEMATRIX and SPARSE packages are completely incompatible: use one of the other!

Adding support for hash tables to REDUCE on Common Lisp inspired me to try to write a better package to support sparse matrices.  I developed the SPARSEMATRIX package entirely using REDUCE on Steel Bank Common Lisp (on Windows), and used other version of REDUCE and other platforms to test portability.

The SPARSE package has a number of issues:
* If sparse and dense matrices are both used in an expression then the result appears always to be a sparse matrix.  I think that the result should usually be a dense matrix [BUT THIS NEEDS CHECKING]!
* The REDUCE operators `map`, `sub`, `cofactor` and `nullspace` are not supported.  The aggregate property (that appropriate operators automatically map over a data structure) is not supported.  Inverses can only be computed for numerical matrices.  Raising a matrix to the power 0 produces the inverse.  Non-positive power of matrices in products fail.
* The `mateigen` operator is implemented as `spmateigen`, rather than overloading `mateigen`.

## Supported sparse matrix operations

In the following description, `m` represents a dense matrix and `s` represents a sparse matrix (or a variable to which such a matrix is assigned).  Similarly, `i` and `j` represent positive integers.

Description | Dense matrix operation | Sparse matrix operation
------------|------------------------|------------------------
Declaration | `matrix m(i,j)` | `sparse_matrix s(i,j)`
Creation | `m := mat(...)` | `sparse_random_matrix(s,i,j)`
Dimensions | `length m` | `length s`
Display (with `on nat`) | 2D table | list of nonzero elements
Element access | `m(i,j)` | `s(i,j)`
Aggregate | `<function> m` | `<function> s`
Substitution | `sub(<equations>, m)` | `sub(<equations>, s)`
Trace | `trace m` | `sparse_trace s`
Transpose | `tp m` | `sparse_tp s`
Determinant | `det m` | `sparse_det s`
Arithmetic | `+ - * / ^` | `+ - * / ^`
Inverse | `m^(-1)` | `s^(-1)`
Rank | `rank m` | `sparse_rank s`
Cofactors | `cofactor(m,i,j)` | `sparse_cofactor(s,i,j)`
Map | `map(1 + ~w, m)` | `map(1 + ~w, s)`

* `densify`: converts a sparse matrix to a dense matrix
* `sparsify`: converts a dense matrix to a sparse matrix

## Planned sparse matrix operations

Overload all standard matrix operators and allow combinations of dense and sparse matrices.

* `MATEIGEN` Operator
* `NULLSPACE` Operator
* Operators in `LINALG` package (maybe)
* Operators in `SPARSE` package (maybe)


<!-- Local Variables: -->
<!-- fill-column: 1000 -->
<!-- eval: (auto-fill-mode -1) -->
<!-- eval: (visual-line-mode 1) -->
<!-- eval: (visual-wrap-prefix-mode 1) -->
<!-- End: -->
