# SPARSEMATRIX: An alternative REDUCE sparse matrix package

**[Francis Wright](https://sites.google.com/site/fjwcentaur)**<br/>
Time-stamp: <2026-07-21 16:44:47 franc>

**This package is currently under development and hence experimental.**

It currently runs correctly on Common Lisp (only tested on SBCL so far!) and PSL, but not on CSL revision 7366 or earlier due to a bug in the hash-table support, which I hope will be fixed soon.

* [Introduction](#introduction)
* [Supported matrix operations](#supported-matrix-operations)
* [Predicates](#predicates)
* [Support for `LINALG` operators](#support-for-linalg-operators)
* [SPARSE: the original REDUCE sparse matrix package](#sparse-the-original-reduce-sparse-matrix-package)
* [TO DO](#to-do)


## Introduction

A [*sparse matrix*](https://en.wikipedia.org/wiki/Sparse_matrix) is a matrix in which most of the elements are zero.  Common examples of sparse matrices are [diagonal](https://en.wikipedia.org/wiki/Diagonal_matrix) and [band](https://en.wikipedia.org/wiki/Band_matrix) matrices.  By contrast, if most of the elements are nonzero, the matrix is considered to be *dense*.  Sparse matrices benefit from being stored using different data structures and manipulated using different algorithms from dense matrices.  Whether it is more efficient to regard a matrix (or more likely a set of matrices) as dense or sparse is ill defined and probably depends on the context, so it may be determinable only by experiment, but it is reasonable to assume that in a sparse matrix no more than about half the elements are nonzero.  A common borderline case is triangular matrices.  (Of course, a dense matrix can be treated as a sparse matrix, and vice versa, which is likely to be non-optimal but is useful for testing.)

The standard REDUCE MATRIX package implicitly assumes dense matrices.  The SPARSEMATRIX package is a re-implementation of the MATRIX package to support sparse matrices.  It uses hash-tables to store matrix elements, and the canonical form for a sparse matrix (henceforth referred to as *sparse representation*) is a LISP list of the form `(<hash> <m> <n>)`, where `<hash>` is a hash-table, `<m>` is the number of rows (row dimension), and `<n>` is the number of columns (column dimension).  Only nonzero elements should be stored in the hash-table and all missing elements are implicitly zero.  By contrast, the canonical form for a dense matrix (henceforth referred to as *dense representation*) is a LISP list of rows of the form `(<row_1> <row_2> ... <row_m>)`, where each `<row_i>` is a list of the matrix elements in the *i-th* row.  The algorithms used in the SPARSEMATRIX package try to avoid accessing implicitly-zero (i.e. non-stored) matrix elements, whereas the algorithms used in the MATRIX package always run through all matrix elements.  Results obtained using the `MATRIX` and `SPARSEMATRIX` packages should be identical (apart from memory use and time), and the test files compare the results of using them both.

In tests (on Steel Bank Common Lisp) using large sparse matrices (500&times;500 matrices with 1000 nonzero rational number elements &ndash; 0.4% density), the `SPARSEMATRIX` package is very much faster than the `MATRIX` package for addition and multiplication, and faster for inversion and determinant computation.

I wrote this package specifically to experiment with hash-table support in REDUCE on Common Lisp, and then ported it to PSL and CSL.  The code naturally relies heavily on lexical scoping, which Common Lisp uses but Standard Lisp (hence PSL/CSL) does not, so I emulate it in PSL/CSL by declaring variables that should be lexically scoped to be `fluid`.  I based the hash-table access on the functions provided by Common Lisp, which are different in PSL/CSL (and different between PSL and CSL), so I have implemented any missing Common Lisp functions that I need using what is available in PSL/CSL.


## Supported matrix operations

In the following description, the letters `d` and `s` represent variables to which matrices in respectively dense or sparse representation will be, or have been, assigned.  Similarly, the letters `i` and `j` represent positive integers.

Description | Dense matrix operation | Sparse matrix operation
------------|------------------------|------------------------
Declaration | `matrix d(i,j), ...` | `sparse_matrix s(i,j), ...`
Creation | `d := mat(...)` | `s := sparse_random_matrix(i,j,<types>)`

In the following description, the letter `m` represents any valid *matrix expression* involving matrices using *dense and/or sparse* representation, and the letter `s` represents any valid *matrix expression* involving matrices using only *sparse* representation.  (In other words, all standard REDUCE matrix operations are generic and accept any mixture of matrix types, but there are also special cases of some operators that accept only matrices in sparse representation, which are provided primarily to facilitate testing.)

Description | Generic matrix operation
------------|-------------------------
Element access | `m(i,j)`
Arithmetic | `+ - * / ^`
Inverse | `m^(-1)`
Dimensions | `length m`
Determinant | `det m` (or `sparse_det s`)
Eigenvectors | `mateigen(m,id)` (or `sparse_mateigen(s,id)`)
Transpose | `tp m` (or `sparse_tp s`)
Trace | `trace m` (or `sparse_trace s`)
Cofactors | `cofactor(m,i,j)` (or `sparse_cofactor(s,i,j)`)
Nullspace | `nullspace m` (or `sparse_nullspace s`)
Rank | `rank m` (or `sparse_rank s`)
Substitution | `sub(<equations>, m)`
Explicit mapping | `map(<expression>, m)`
Implicit mapping | `<function> m`
Density | `matrix_density m`

### Implicit matrix type conversion

Matrices in dense and sparse representation can be freely mixed in algebraic expressions, i.e. on either side of the binary operators `+ - * /`.  The type of an algebraic expression involving matrices in both representations currently defaults to dense, but can be changed by the operator `sparse_matrix_auto_convert_type`.  If an argument is supplied then it must be one of the symbols `dense`, `sparse` or `none`, and the operator always returns the previous type.  Hence, `sparse_matrix_auto_convert_type();` just displays the current auto-conversion type, without changing it.

### Explicit matrix type conversion

* `densify` converts a matrix from sparse to dense representation, e.g. `d := densify s`
* `sparsify` converts a matrix from dense to sparse representation, e.g. `s := sparsify d`

Note that using explicit matrix type conversion in an expression disables implicit matrix type conversion for that expression.

### Matrix output

If the switch `sparse_matrix_dense_print` is `on`, which it is by default, then small matrices in sparse representation are displayed the same as matrices in dense representation (mainly to facilitate testing), where _small_ means having no more than the number of columns specified by the value of the (shared) variable `sparse_matrix_dense_print_colmax`, which is 10 by default.

See `sparsematrix.rlg` for examples of using the above `SPARSEMATRIX` versions of the `MATRIX` operators with small (dense) matrices and the files `speed*.tst` for some timed examples using large sparse matrices.


## Predicates

The first three of these predicates mirror predicates in the `LINALG` package.

* `sparse_matrix_p` (cf. `LINALG` `matrixp`)
* `sparse_square_matrix_p` (cf. `LINALG` `squarep`)
* `sparse_symmetric_matrix_p` (cf. `LINALG` `symmetricp`)
* `sparse_skew_symmetric_matrix_p`
* `sparse_hermitian_matrix_p`
* `sparse_skew_hermitian_matrix_p`
* `sparse_diagonal_matrix_p`
* `sparse_upper_triangular_matrix_p`
* `sparse_lower_triangular_matrix_p`
* `sparse_identity_matrix_p`
* `sparse_orthogonal_matrix_p`
* `sparse_unitary_matrix_p`

These predicates all take a single argument that can be anything.  They return `true` if the argument evaluates to a sparse matrix with the property implied by the name of the predicate, and `false` otherwise.  Currently, these predicates are *exact* and do not allow for numerical error in floating-point matrices!  See `sparsepredicates.rlg` for examples of using the above predicates.


## Support for `LINALG` operators

The following potentially useful functions for working with sparse matrices are modelled loosely on functions in `LINALG`, the REDUCE Linear Algebra Package, by Matt Rebbeck.

### Matrix construction:

* `sparse_random_matrix` (cf. `LINALG` `random_matrix`)
* `sparse_identity_matrix` (cf. `LINALG` `make_identity`)
* `sparse_band_matrix` (cf. `LINALG` `band_matrix`, but with reversed arguments)

### Whole matrix manipulation:

* `sparse_matrix_augment` (cf. `LINALG` `matrix_augment`)
* `sparse_block_diagonal_matrix` (cf. `LINALG` `diagonal`)

### Column manipulation:

* `sparse_select_columns` (synonym `sparse_augment_columns`, cf. `LINALG`  `augment_columns`)
* `sparse_remove_columns` (cf. `LINALG` `remove_columns`)
* `sparse_get_columns` (cf. `LINALG` `get_columns`)

These operators are all more general than those in the `LINALG` package.  The input matrices can use either sparse or dense representation, but the output always uses sparse representation.  Arguments specifying column indices can be a sequence of integers, integer lists or integer intervals, which can be freely intermixed.  Column indices can be negative, meaning count from the right, and intervals can be descending.  The operator `sparse_select_columns` allows columns to be duplicated.  See `sparselinalg.rlg` for examples of using the above `SPARSEMATRIX` versions of the `LINALG` operators.


## SPARSE: the original REDUCE sparse matrix package

This was written by Stephen Scowcroft in 1995 at the Konrad-Zuse-Zentrum für Informationstechnik Berlin.  It uses a different data structure to represent a sparse matrix, namely a LISP list of the form `(<vector> (spm <m> <n>))`, where `<vector>` is a vector of rows of the form `[nil <row_1> <row_2> ... <row_m>]` and each `<row_i>` is a list of the matrix elements in the *i-th* row, each represented as a dotted pair, of the form `((nil) (j_1 . val_1) (j_2 . val_2) ...)`.  The `nil` elements are to allow for the fact that REDUCE indexes vectors and lists from 0, whereas matrix indices conventionally start from 1.  **Note that the SPARSEMATRIX and SPARSE packages are completely incompatible: use one of the other!**

In tests using large sparse matrices (500&times;500 matrices with 1000 nonzero rational number elements &ndash; 0.4% density), the `SPARSE` and `SPARSEMATRIX` packages are broadly comparable in speed.  However, the SPARSE package has a number of issues:
* Computing the determinant of a sparse 50*50 integer matrix with 100 nonzero elements fails with heap overflow.
* Inverses can only be computed for numerical matrices.  Raising a matrix to the power 0 produces the inverse.  Non-positive powers of matrices in products fail.
* The REDUCE operators `map`, `sub`, `cofactor` and `nullspace` are not supported.  The aggregate property (that appropriate operators automatically map over a data structure) is not supported.
* If matrices in sparse and dense representation are both used in an expression then the result appears always to use sparse representation, whereas I think that the result should usually be a dense matrix.
* The `mateigen` operator is implemented only as `spmateigen`, rather than overloading `mateigen`.


## TO DO

* Check predicates apply substitutions.
* Better argument checking.
* More efficient maphash-and on Common Lisp.
* Overall efficiency review -- avoid repeated simplifications, etc.
* Check timings on PSL and CSL.
* Support for special matrices -- triangular, symmetric, etc. -- via access functions.
* More operators from the `LINALG` package (maybe), e.g. row analogues of column manipulations, sub_matrix (i.e. both row and column manipulation).
* More operators from the `SPARSE` package (maybe).
* Operators from the `NORMFORM` package (maybe).

<!-- Local Variables: -->
<!-- fill-column: 10000 -->
<!-- eval: (auto-fill-mode -1) -->
<!-- eval: (visual-line-mode 1) -->
<!-- eval: (visual-wrap-prefix-mode 1) -->
<!-- End: -->
