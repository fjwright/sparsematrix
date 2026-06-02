# SPARSEMATRIX: A REDUCE sparse matrix package

**[Francis Wright](https://sites.google.com/site/fjwcentaur)**<br/>
Time-stamp: <2026-06-02 17:30:13 franc>

A [*sparse matrix*](https://en.wikipedia.org/wiki/Sparse_matrix) is a matrix in which most of the elements are zero.  Common examples of sparse matrices are [diagonal](https://en.wikipedia.org/wiki/Diagonal_matrix) and [band](https://en.wikipedia.org/wiki/Band_matrix) matrices.  By contrast, if most of the elements are non-zero, the matrix is considered *dense*.  Sparse matrices benefit from being stored using different data structures and manipulated using different algorithms from dense matrices.  Whether it is more efficient to regard a matrix (or more likely a set of matrices) as dense or sparse is ill defined and probably depends on the context, so it may be determinable only by experiment, but it is reasonable to assume that in a sparse matrix no more than half the elements are nonzero.  A common borderline case is triangular matrices.  (Of course, a dense matrix can be treated as a sparse matrix, and vice versa, which is likely to be less efficient but is useful for testing.)

The REDUCE MATRIX package implicitly assumes dense matrices.  SPARSEMATRIX is a re-implementation of the MATRIX package that assumes sparse matrices.  It uses hash tables to store matrix elements, and the canonical form for a sparse matrix is a LISP list of the form `(<hash> <m> <n>)`, where `<hash>` is a hash table, `<m>` is the number of rows (row dimension), and `<n>` is the number of columns (column dimension).  Only nonzero elements should be stored in the hash table and all missing elements are implicitly zero.  By contrast, the canonical form for a dense matrix is a LISP list of rows of the form `(<row_1> <row_2> ... <row_m>)`, where each `<row_i>` is a list of the matrix elements in the *i-th* row.  The algorithms used to manipulate sparse matrices avoid accessing implicitly-zero (i.e. non-stored) matrix elements as much as possible, whereas the algorithms used to manipulate dense matrices always run through all matrix elements, regardless of their values.

In tests using large sparse matrices (500&times;500 matrices with 1000 nonzero rational number elements), the `SPARSEMATRIX` and `SPARSE` packages are comparable and both very significantly faster than the `MATRIX` package for addition and multiplication.  They are faster for inversion and `SPARSEMATRIX` is faster for determinant (but `SPARSE` crashes).

Currently, sparse and dense matrices cannot be mixed within algebraic expressions, but I plan to implement full and automatic interoperability.  However, a sparse matrix can be explicitly converted to a dense matrix, and vice versa.

## Critique of SPARSE, an alternative sparse matrix package

There is already a REDUCE package called SPARSE to support sparse matrices, written by Stephen Scowcroft in 1995 at the Konrad-Zuse-Zentrum für Informationstechnik Berlin.  It uses a different data structure to represent a sparse matrix, namely a LISP list of the form `(<vector> (spm <m> <n>))`, where `<vector>` is a vector of rows of the form `[nil <row_1> <row_2> ... <row_m>]` and each `<row_i>` is a list of the matrix elements in the *i-th* row, each represented as a dotted pair, of the form `((nil) (j_1 . val_1) (j_2 . val_2) ...)`.  The `nil` elements are to allow for the fact that REDUCE indexes vectors and lists from 0, whereas matrix indices conventionally start from 1.  Note that the SPARSEMATRIX and SPARSE packages are completely incompatible: use one of the other!

Adding support for hash tables to REDUCE on Common Lisp inspired me to try to write a better package using hash tables to support sparse matrices.  I developed the SPARSEMATRIX package entirely using REDUCE on Steel Bank Common Lisp (on Windows), and used other versions of REDUCE and other platforms only to test portability.

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
Element access | `m(i,j)` | `s(i,j)`
Arithmetic | `+ - * / ^` | `+ - * / ^`
Inverse | `m^(-1)` | `s^(-1)`
Transpose | `tp m` | `sparse_tp s`
Trace | `trace m` | `sparse_trace s`
Determinant | `det m` | `sparse_det s`
Cofactors | `cofactor(m,i,j)` | `sparse_cofactor(s,i,j)`
Rank | `rank m` | `sparse_rank s`
Nullspace | `nullspace m` | `sparse_nullspace s`
Substitution | `sub(<equations>, m)` | `sub(<equations>, s)`
Explicit mapping | `map(1 + ~w, m)` | `map(1 + ~w, s)`
Implicit mapping | `<function> m` | `<function> s`

* `densify` converts a sparse matrix to a dense matrix, e.g. `m := densify s`
* `sparsify` converts a dense matrix to a sparse matrix, e.g. `s := sparsify m`

If the switch `sparse_matrix_dense_print` is `on`, which it is by default, then small sparse matrices are displayed the same as dense matrices (mainly to facilitate testing), where _small_ means having no more than the number of columns specified by the value of the (shared) variable `sparse_matrix_dense_print_colmax`, which is 10 by default.

See `sparsematrix.rlg` for examples of using the above `SPARSEMATRIX` versions of the `MATRIX` operators with small (dense) matrices and the files `speed*.tst` for some timed examples using large sparse matrices.

## Predicates

The first three of these predicates mirror predicates in the `LINALG` package.

* `sparse_matrix_p`, cf. `LINALG` `matrixp`
* `sparse_square_matrix_p`, cf. `LINALG` `squarep`
* `sparse_symmetric_matrix_p`, cf. `LINALG` `symmetricp`
* `sparse_skew_symmetric_matrix_p`
* `sparse_hermitian_matrix_p`
* `sparse_skew_hermitian_matrix_p`
* `sparse_diagonal_matrix_p`
* `sparse_upper_triangular_matrix_p`
* `sparse_lower_triangular_matrix_p`
* `sparse_identity_matrix_p`
* `sparse_orthogonal_matrix_p`
* `sparse_unitary_matrix_p`

These predicates all take a single argument that can be anything.  They return `true` if the argument evaluates to a sparse matrix with the property implied by the name of the predicate, and `false` otherwise.

Currently, support is in the file `sparsepredicates.red`, which needs to be input separately from the main `SPARSEMATRIX` package.  See `sparsepredicates.rlg` for examples of using the above predicates.


## Support for `LINALG` operators

Some potentially useful facilities for working with sparse matrices modelled loosely on `LINALG`, the REDUCE Linear Algebra Package, by Matt Rebbeck.

### Matrix construction:

* `sparse_identity_matrix`, cf. `LINALG` `make_identity`
* `sparse_band_matrix`, cf. `LINALG` `band_matrix` (with reversed arguments)

### Whole matrix manipulation:

* `sparse_matrix_augment`, cf. `LINALG` `matrix_augment`
* `sparse_block_diagonal_matrix`, cf. `LINALG` `diagonal`

### Column manipulation:

* `sparse_select_columns` (synonym `sparse_augment_columns`), cf. `LINALG`  `augment_columns`
* `sparse_remove_columns`, cf. `LINALG` `remove_columns`
* `sparse_get_columns`, cf. `LINALG` `get_columns`

These operators are all more general than those in the `LINALG` package.  The input matrices can use either sparse or dense representation, but the output always uses sparse representation.  Arguments specifying column indices can be a sequence of integers, integer lists or integer intervals, which can be freely intermixed.  Column indices can be negative, meaning count from the right, and intervals can be descending.  The operator `sparse_select_columns` allows columns to be duplicated.

Currently, support is in the file `sparselinalg.red`, which needs to be input separately from the main `SPARSEMATRIX` package.  See `sparselinalg.rlg` for examples of using the above `SPARSEMATRIX` versions of the `LINALG` operators.


## TO DO

* Overload all standard matrix operators and allow combinations of dense and sparse matrices in algebraic expressions.
* Support for special matrices -- triangular, symmetric, etc. -- via access functions.
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
