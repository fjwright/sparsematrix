module sparselinalg; % Construction and manipulation of sparse matrices

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-06-26 12:07:36 franc>
% Created: May 2026

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%  * Redistributions of source code must retain the relevant copyright
%    notice, this list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the relevant
%    copyright notice, this list of conditions and the following
%    disclaimer in the documentation and/or other materials provided
%    with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% $Id$

#if (not (memq 'common!-lisp lispsystem!*))
fluid '(hash!* n!* dim!* alist!* columns!* hashes!*);
#endif

% Some potentially useful facilities for working with sparse matrices
% modelled loosely on LINALG, the REDUCE Linear Algebra Package, by
% Matt Rebbeck.

symbolic procedure safe!-cdr x;
   % This function is not currently provided in mainstream REDUCE.
   % (The file "alg/general.red" containing this definition is not
   % included in the build.)  It is provided in REDUCE on Common Lisp
   % and flagged lose.  (It is aliased to cdr, because Common Lisp cdr
   % is safe.)
   if atom x then nil else cdr x;

symbolic procedure sparse!-reval!&flatten u;
   % Reval and flatten a list that may include REDUCE lists.
   for each el in revlis u join
      if eqcar(el, 'list) then cdr el else {el};

%                         %%%%%%%%%%%%%%%%%%%
%                         MATRIX CONSTRUCTION
%                         %%%%%%%%%%%%%%%%%%%

% sparse_identity_matrix
% cf. LINALG make_identity

put('sparse_identity_matrix, 'psopfn, 'sparse_identity_matrix);
put('sparse_identity_matrix, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_identity_matrix u; % (dim)
   % DIM is a positive integer.  Return a DIM*DIM sparse identity
   % matrix.
   begin scalar dim := reval carx(u, "sparse_identity_matrix"), hash;
      if not(fixp dim and dim > 0) then typerr(dim, "matrix dimension");
      hash := mk!-sparse!-matrix!-hash();
      for i := 1 : dim do puthash(i.i, hash, 1);
      return {'sparse!-mat, hash, dim, dim}
   end;

% sparse_band_matrix
% cf. LINALG band_matrix (with reversed arguments)

put('sparse_band_matrix, 'psopfn, 'sparse_band_matrix);
put('sparse_band_matrix, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_band_matrix u; % (dim, scalars)
   % DIM is a positive integer specifying the dimensions of the square
   % matrix generated.  SCALARS is a sequence of scalar expressions,
   % or lists thereof.  Return a DIM*DIM sparse matrix with the n
   % scalars in the order specified in each row i in columns from j =
   % i - fix((n-1)/2) to j + n.  Normally, n will be odd.
   u and cdr u and
   begin scalar dim := reval car u, n, j0, hash;
      if not(fixp dim and dim > 0) then typerr(dim, "matrix dimension");
      for each el in (u := sparse!-reval!&flatten cdr u) do
         if getrtype el then typerr(el, "scalar");
      n := length u;  j0 := quotient(n-1, 2);
      hash := mk!-sparse!-matrix!-hash();
      for i := 1 : dim do
         begin scalar j := i - j0;
            for each x in u do <<
               if 1 <= j and j <= dim and
                  x neq 0 then puthash(i.j, hash, x);
               j := j + 1;
            >>;
         end;
      return {'sparse!-mat, hash, dim, dim}
   end;


%                      %%%%%%%%%%%%%%%%%%%%%%%%%
%                      WHOLE MATRIX MANIPULATION
%                      %%%%%%%%%%%%%%%%%%%%%%%%%

% sparse_matrix_augment
% cf. LINALG matrix_augment

put('sparse_matrix_augment, 'psopfn, 'sparse_matrix_augment);
put('sparse_matrix_augment, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_matrix_augment u;
   % Accept any number of arguments that are either sparse or dense
   % matrices, or lists thereof, where all matrices must have the same
   % number of rows.  The result is the matrices adjoined horizontally
   % in the order specified, as a sparse matrix.
   u and begin scalar m;   % common row dimension
      % Build a list of sparse matrix canonical forms:
      u := for each el in sparse!-reval!&flatten u collect sparse!-matsm
         begin scalar newm, mtrx :=
            if eqcar(el, 'sparse!-mat) then <<
               newm := caddr el;        % row dim of sparse matrix
               el
            >> else if eqcar(el, 'mat) then <<
               newm := length cdr el;   % row dim of dense matrix
               sparsify el
            >> else typerr(el, "matrix");
            if m then
               (if newm neq m then rederr
                  "matrices must have the same row dimensions")
            else m := newm;
            return mtrx;
         end;
      return sparse!-matsm!*1 sparse!-matrix!-augment(u, m);
   end;

symbolic procedure sparse!-matrix!-augment(u, m);
   % Adjoin horizontally a list U of sparse matrix canonical forms all
   % with row dim M, and return a sparse matrix canonical form.
   begin scalar hash!* := copyhash caar u,
         n!* := caddar u;               % col dim of result so far
      for each el in cdr u do <<
         maphash(function
            (lambda(key, value);
            puthash(car key . (n!* + cdr key), hash!*, value)),
            car el);
         n!* := n!* + caddr el;
      >>;
      return {hash!*, m, n!*};
   end;

% sparse_block_diagonal_matrix
% cf. LINALG diagonal

put('sparse_block_diagonal_matrix, 'psopfn, 'sparse_block_diagonal_matrix);
put('sparse_block_diagonal_matrix, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_block_diagonal_matrix u;
   % Accept any number of arguments that are either sparse or dense
   % matrices or scalars, or lists thereof, where all matrices must be
   % square.  Scalars are treated as 1*1 matrices.  The result is the
   % matrices adjoined into a block diagonal matrix in the order
   % specified, as a square sparse matrix.
   u and <<
      % Build a list of sparse matrix canonical forms or SQs:
      u := for each el in sparse!-reval!&flatten u collect
         begin scalar rtype := getrtype el;
            if null rtype then return simp el;
            if rtype eq 'sparse!-matrix then
               (if caddr el = caddr el then % square
                  return sparse!-matsm el)
            else if rtype eq 'matrix then
               (if length cdr el = length cadr el then % square
                  return sparse!-matsm sparsify el);
            typerr(el, "square matrix");
         end;
      sparse!-matsm!*1 sparse!-block!-diagonal!-matrix u
   >>;

symbolic procedure sparse!-block!-diagonal!-matrix u;
   % Adjoin diagonally a list U of square sparse matrix canonical
   % forms or SQs, and return a square sparse matrix canonical form.
   begin scalar hash!* := mk!-sparse!-matrix!-hash();
      integer dim!*;                  % row & col dim of result so far
      for each el in u do
         if hash!-table!-p car el then << % sparse matrix canonical form
            maphash(function
               (lambda(key, value);
               puthash((dim!* + car key) . (dim!* + cdr key), hash!*, value)),
               car el);
            dim!* := dim!* + caddr el;
         >> else <<                     % scalar = 1*1 matrix
            dim!* := dim!* + 1;
            puthash(dim!* . dim!*, hash!*, el);
         >>;
      return {hash!*, dim!*, dim!*};
   end;


%                         %%%%%%%%%%%%%%%%%%%
%                         COLUMN MANIPULATION
%                         %%%%%%%%%%%%%%%%%%%

% sparse_select_columns / sparse_augment_columns
% cf. LINALG augment_columns

put('sparse_select_columns, 'psopfn, 'sparse_select_columns);
put('sparse_augment_columns, 'psopfn, 'sparse_select_columns);
put('sparse_select_columns, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_select_columns u; % (mtrx, columns)
   % MTRX should be an algebraic sparse or dense matrix.  COLUMNS
   % should be a sequence of integers, integer lists or integer
   % intervals representing column indices, i.e. `col_1, col_2, ...,
   % col_n', or `{col_1, col_2, ..., col_n}' or `col_1 .. col_n',
   % which means all column indices in the interval from col_1 to
   % col_n inclusive (where `..' is the REDUCE interval operator).
   % Negative indices are allowed, and count from the right.  In an
   % interval `a .. b', if b < a then the interval is expanded as `a,
   % a-1, a-2, ..., b'.  Return an algebraic sparse matrix copy of
   % MTRX containing only the specified columns in the order
   % specified; duplicate column indices are respected.
   u and sparse!-process!-mtrx!&cols(u, function sparse!-select!-columns);

symbolic procedure sparse!-index!-check(idx, n);
   if fixp idx and not zerop idx and abs idx <= n then
      if minusp idx then n + 1 + idx else idx
   else typerr(idx, "matrix column index");

symbolic procedure sparse!-process!-mtrx!&cols(u, fn);
   begin scalar mtrx, n, columns;
      u := revlis u;
      % Process matrix:
      mtrx := car u;
      if eqcar(mtrx, 'mat) then
         mtrx := sparsify mtrx
      else if not eqcar(mtrx, 'sparse!-mat) then
         typerr(mtrx, "matrix");
      mtrx := sparse!-matsm mtrx;       % sparse matrix canonical form
      n := caddr mtrx;                  % column dimension
      % Process column indices:
      columns := for each el in cdr u join
         if fixp el then {sparse!-index!-check(el, n)}
         else if eqcar(el, 'list) then
            for each col in cdr el collect
               sparse!-index!-check(col, n)
         else if eqcar(el, '!*interval!*) then
         begin scalar
            a := sparse!-index!-check(cadr el, n),
            b := sparse!-index!-check(caddr el, n),
            s := if a <= b then +1 else -1;
            return for col := a step s until b collect col;
         end;
      return apply2(fn, mtrx, columns);
   end;

symbolic procedure sparse!-select!-columns(mtrx, columns);
   % MTRX is a sparse matrix canonical form.  COLUMNS is a list of
   % column indices.  Return an algebraic copy of MTRX containing only
   % the specified columns, with the columns correctly re-indexed.
   begin scalar hash!* := mk!-sparse!-matrix!-hash(), alist!*;
      % Convert columns to an alist with elements of the form
      % (old_col_ind new_col_ind_1 new_col_ind_2 ...):
      integer newcol;                   % initialised to 0
      for each oldcol in columns do
      begin scalar el;
         newcol := newcol + 1;
         if el := assoc(oldcol, alist!*) then
            nconc(el, {newcol})
         else alist!* := {oldcol, newcol} . alist!*;
      end;
      maphash(function
         (lambda(key, value);
         for each newcol in safe!-cdr assoc(cdr key, alist!*) do
            puthash(car key . newcol, hash!*, value)),
         car mtrx);
      return sparse!-matsm!*1 {hash!*, cadr mtrx, newcol};
   end;

% sparse_remove_columns
% cf. LINALG remove_columns

put('sparse_remove_columns, 'psopfn, 'sparse_remove_columns);
put('sparse_remove_columns, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_remove_columns u; % (mtrx, columns)
   % MTRX should be an algebraic sparse or dense matrix.  COLUMNS
   % should be a sequence of integers, integer lists or integer
   % intervals representing column indices, i.e. col_1, col_2, ...,
   % col_n, or {col_1, col_2, ..., col_n} or col_1 .. col_n, which
   % means all column indices in the interval from col_1 to col_n
   % inclusive (where `..' is the REDUCE interval operator).  Negative
   % indices are allowed, and count from the right.  In an interval `a
   % .. b', if b < a then the interval is expanded as `a, a-1, a-2,
   % ..., b'.  Return an algebraic sparse matrix copy of MTRX without
   % the specified columns.
   u and sparse!-process!-mtrx!&cols(u, function sparse!-remove!-columns);

symbolic procedure sparse!-remove!-columns(mtrx, columns!*);
   % MTRX is a sparse matrix canonical form.  COLUMNS is a list of
   % column indices.  Return an algebraic copy of MTRX without the
   % specified columns, with the columns correctly re-indexed.
   begin scalar hash!* := mk!-sparse!-matrix!-hash();
      integer newcol;           % initialised to 0
      % Convert columns to a selection list:
      columns!* := for col := 1 : caddr mtrx join
         if not member(col, columns!*) then {col};
      % Convert columns to an alist of elements of the form
      % (old_col_ind . new_col_ind):
      columns!* := for each oldcol in columns!* collect
         (oldcol . (newcol := newcol + 1));
      maphash(function
         (lambda(key, value);
          begin scalar el;
             if (el := assoc(cdr key, columns!*)) then
                puthash(car key . cdr el, hash!*, value);
          end),
         car mtrx);
      return sparse!-matsm!*1 {hash!*, cadr mtrx, newcol};
   end;

% sparse_get_columns
% cf. LINALG get_columns

put('sparse_get_columns, 'psopfn, 'sparse_get_columns);

symbolic procedure sparse_get_columns u; % (mtrx, columns)
   % Like sparse_select_columns, but returns a list of sparse column
   % matrices.
   u and sparse!-process!-mtrx!&cols(u, function sparse!-get!-columns);

symbolic procedure sparse!-get!-columns(mtrx, columns);
   % MTRX is a sparse matrix canonical form.  COLUMNS is a list of
   % column indices.  Return an algebraic list of the specified
   % columns of MTRX, in the order specified, as algebraic sparse
   % column matrices.
   begin scalar hashes!*, alist!*, m := cadr mtrx;
      % Convert columns to an alist with elements of the form
      % (old_col_ind new_col_ind_1 new_col_ind_2 ...):
      integer newcol;                   % initialised to 0
      hashes!* := for each oldcol in columns collect
      begin scalar el;
         newcol := newcol + 1;
         if el := assoc(oldcol, alist!*) then
            nconc(el, {newcol})
         else alist!* := {oldcol, newcol} . alist!*;
         return mk!-sparse!-matrix!-hash();
      end;
      maphash(function
         (lambda(key, value);
         for each newcol in safe!-cdr assoc(cdr key, alist!*) do
            puthash(car key . 1, nth(hashes!*, newcol), value)),
         car mtrx);
      return 'list . for each hash in hashes!* collect
         sparse!-matsm!*1 {hash, m, 1};
   end;

endmodule;

end;
