# SPARSEMATRIX: A REDUCE sparse matrix package

**[Francis Wright](https://sites.google.com/site/fjwcentaur)**<br/>
Time-stamp: <2026-05-26 18:13:47 franc>

A [*sparse matrix*](https://en.wikipedia.org/wiki/Sparse_matrix) is a matrix in which most of the elements are zero.  Common examples of sparse matrices are [diagonal](https://en.wikipedia.org/wiki/Diagonal_matrix) and [band](https://en.wikipedia.org/wiki/Band_matrix) matrices.  By contrast, if most of the elements are non-zero, the matrix is considered *dense*.  Sparse matrices benefit from being stored using different data structures and manipulated using different algorithms from dense matrices.  Whether it is more efficient to regard a matrix (or more likely a set of matrices) as dense or sparse is ill defined and probably depends on context, so it may be determinable only by experiment, but it is reasonable to assume that in a sparse matrix fewer than half the elements are nonzero.  (Of course, a dense matrix can be treated as a sparse matrix, and vice versa, which is likely to be less efficient but is useful for testing.)

The REDUCE MATRIX package implicitly assumes dense matrices.  SPARSEMATRIX is a re-implementation of the MATRIX package that assumes sparse matrices.  It uses hash tables as the primary data structures to store matrix elements, and the canonical form for a sparse matrix is a LISP list of the form `(<hash> <m> <n>)`, where `<hash>` is a hash table, `<m>` is the number of rows (row dimension), and `<n>` is the number of columns (column dimension).  Only nonzero elements are ever stored in the has table; missing elements are implicitly zero.  By contrast, the canonical form for a dense matrix is a LISP list of rows of the form `(<row_1> <row_2> ... <row_m>)`, where each `<row_i>` is a list of the matrix elements in the *i-th* row.  The algorithms used to manipulate sparse matrices avoid accessing implicitly-zero (i.e. non-stored) matrix elements as much as possible, whereas the algorithms used to manipulate dense matrices always run through all matrix elements, regardless of their values.

In tests using large sparse matrices (500&times;500 matrices with 1000 nonzero rational number elements), the `SPARSEMATRIX` and `SPARSE` packages are comparable and both very significantly faster than the `MATRIX` package for addition and multiplication.  They are faster for inversion and `SPARSEMATRIX` is faster for determinant (but `SPARSE` crashes).

Currently, sparse and dense matrices cannot generally be mixed, but I plan to implement full and automatic interoperability in future.  However, a sparse matrix can be explicitly converted to a dense matrix, and vice versa.

## Critique of SPARSE, an alternative sparse matrix package

There is already a REDUCE package called SPARSE to support sparse matrices, written by Stephen Scowcroft in 1995 at the Konrad-Zuse-Zentrum für Informationstechnik Berlin.  It uses a different data structure to represent a sparse matrix, namely a LISP list of the form `(<vector> (spm <m> <n>))`, where `<vector>` is a vector of rows of the form `[nil <row_1> <row_2> ... <row_m>]` and each `<row_i>` is a list of the matrix elements in the *i-th* row, each represented as a dotted pair, of the form `((nil) (j_1 . val_1) (j_2 . val_2) ...)`.  The `nil` elements are to allow for the fact that REDUCE indexes for vectors and lists start from 0, whereas matrix indices start from 1.  Note that the SPARSEMATRIX and SPARSE packages are completely incompatible: use one of the other!

Adding support for hash tables to REDUCE on Common Lisp inspired me to try to write a better package using hash tables to support sparse matrices.  I developed the SPARSEMATRIX package entirely using REDUCE on Steel Bank Common Lisp (on Windows), and used other versions of REDUCE and other platforms to test portability.

The SPARSE package has a number of issues:
* Computing the determinant of a sparse 50*50 integer matrix with 100 nonzero elements fails with heap overflow.
* Inverses can only be computed for numerical matrices.  Raising a matrix to the power 0 produces the inverse.  Non-positive powers of matrices in products fail.
* The REDUCE operators `map`, `sub`, `cofactor` and `nullspace` are not supported.  The aggregate property (that appropriate operators automatically map over a data structure) is not supported.
* If sparse and dense matrices are both used in an expression then the result appears always to be a sparse matrix.  I think that the result should usually be a dense matrix [BUT THIS NEEDS CHECKING]!
* The `mateigen` operator is implemented as `spmateigen`, rather than overloading `mateigen`.

## Supported SPARSEMATRIX operations

In the following description, `m` represents a dense matrix and `s` represents a sparse matrix (or a variable to which such a matrix is assigned).  Similarly, `i` and `j` represent positive integers.

Description | Dense matrix operation | Sparse matrix operation
------------|------------------------|------------------------
Declaration | `matrix m(i,j)` | `sparse_matrix s(i,j)`
Creation | `m := mat(...)` | `sparse_random_matrix s(i,j)`
Dimensions | `length m` | `length s`
Display (with `on nat`) | 2D table | generally list of nonzero elements (but 2D table for small matrices)
Element access | `m(i,j)` | `s(i,j)`
Map | `map(1 + ~w, m)` | `map(1 + ~w, s)`
Aggregate | `<function> m` | `<function> s`
Substitution | `sub(<equations>, m)` | `sub(<equations>, s)`
Trace | `trace m` | `sparse_trace s`
Transpose | `tp m` | `sparse_tp s`
Determinant | `det m` | `sparse_det s`
Cofactors | `cofactor(m,i,j)` | `sparse_cofactor(s,i,j)`
Inverse | `m^(-1)` | `s^(-1)`
Arithmetic | `+ - * / ^` | `+ - * / ^`
Rank | `rank m` | `sparse_rank s`
Nullspace | `nullspace m` | `sparse_nullspace s`

* `densify`: converts a sparse matrix to a dense matrix
* `sparsify`: converts a dense matrix to a sparse matrix

Small sparse matrices (with no more than 10 columns) are displayed the same as dense matrices (mainly to facilitate testing).

See `sparsematrix.rlg` for examples of using the above `SPARSEMATRIX` versions of the `MATRIX` operators with small dense matrices that the files `speed*.tst` for some timed examples using large sparse matrices.


## Support for `LINALG` operators

* `sparse_matrix_augment`, cf. `matrix_augment`
* `sparse_select_columns` (synonym `sparse_augment_columns`), cf.  `augment_columns`
* `sparse_remove_columns`, cf. `remove_columns`

The `SPARSEMATRIX` versions of these operators are more general than those in the `LINALG` package.  The input matrices can use either sparse or dense representation, but the output is always uses sparse representation.  Arguments specifying column indices can be sequences of integers, integer lists or integer intervals.  Column indices can be negative, meaning count from the right, and intervals can be descending.  The operator `sparse_select_columns` allows columns to be duplicated.

Currently, support is in the file `sparselinalg.red`, which needs to be input separately from the main `SPARSEMATRIX` package.  See `sparselinalg.rlg` for examples of using the above `SPARSEMATRIX` versions of the `LINALG` operators.

## Planned SPARSEMATRIX support

Overload all standard matrix operators and allow combinations of dense and sparse matrices.

* `MATEIGEN` operator (maybe)
* More operators from `LINALG` package (maybe)
* More operators from `SPARSE` package (maybe)
* Operators from `NORMFORM` package (maybe)

<!-- Local Variables: -->
<!-- fill-column: 1000 -->
<!-- eval: (auto-fill-mode -1) -->
<!-- eval: (visual-line-mode 1) -->
<!-- eval: (visual-wrap-prefix-mode 1) -->
<!-- End: -->
