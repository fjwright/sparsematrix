module sparselinalg;    % Useful linalg operations for sparse matrices

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-27 12:45:44 franc>
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

% Some potentially useful facilities for working with sparse matrices
% modelled loosely on LINALG, the REDUCE Linear Algebra Package, by
% Matt Rebbeck.

% %%%%%%%%%%%%%%%%%%%%%
% sparse_matrix_augment
% %%%%%%%%%%%%%%%%%%%%%
% cf. LINALG matrix_augment

put('sparse_matrix_augment, 'psopfn, 'sparse_matrix_augment);
put('sparse_matrix_augment, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_matrix_augment u;
   % sparse_matrix_augment accepts any number of arguments that are
   % either sparse or dense matrices, or lists thereof, where all
   % matrices must have the same number of rows.  The result is the
   % matrices adjoined horizontally in the order specified, as a
   % sparse matrix.
   begin scalar m;                      % common row dimension
      if null u then return;
      % Flatten any lists among the arguments before processing:
      u := for each el in revlis u join
         if eqcar(el, 'list) then cdr el else {el};
      % Build a list of sparse matrix canonical forms:
      u := for each el in u collect sparse!-matsm
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
   begin scalar hash := copyhash caar u,
         n := caddar u;                 % col dim of result so far
      for each el in cdr u do <<
         maphash(car el,
            (lambda(key, value);
            puthash(car key . (n + cdr key), hash, value)));
         n := n + caddr el;
      >>;
      return {hash, m, n};
   end;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparse_select_columns / sparse_augment_columns
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cf. LINALG augment_columns

put('sparse_select_columns, 'psopfn, 'sparse_select_columns);
put('sparse_augment_columns, 'psopfn, 'sparse_select_columns);
put('sparse_select_columns, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_select_columns u; % (mtrx, columns)
   % MTRX should be a sparse or dense matrix.  COLUMNS should be a
   % sequence of integers, integer lists or integer intervals
   % representing column indices, i.e. `col_1, col_2, ..., col_n', or
   % `{col_1, col_2, ..., col_n}' or `col_1 .. col_n', which means all
   % column indices in the interval from col_1 to col_n inclusive
   % (where `..' is the REDUCE interval operator).  Negative indices
   % are allowed, and count from the right.  In an interval `a .. b',
   % if b < a then the interval is expanded as `a, a-1, a-2, ..., b'.
   % Return a sparse matrix copy of MTRX containing only the specified
   % columns in the order specified; duplicate column indices are
   % respected.
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
      return sparse!-matsm!*1 apply2(fn, mtrx, columns);
   end;

symbolic procedure sparse!-select!-columns(mtrx, columns);
   % MTRX is a sparse matrix canonical form.
   % COLUMNS is a list of column indices.
   % Return a copy of MTRX containing only the specified columns,
   % with the columns correctly re-indexed.
   begin scalar hash := mk!-sparse!-matrix!-hash(), alist,
      % Convert columns to an alist with elements of the form
      % (old_col_ind new_col_ind_1 new_col_ind_2 ...):
      integer newcol;                   % initialised to 0
      for each oldcol in columns do
      begin scalar el;
         newcol := newcol + 1;
         if el := assoc(oldcol, alist) then
            nconc(el, {newcol})
         else alist := {oldcol, newcol} . alist;
      end;
      maphash(car mtrx,
         (lambda(key, value);
         for each newcol in cdr assoc(cdr key, alist) do
            puthash(car key . newcol, hash, value)));
      return {hash, cadr mtrx, newcol};
   end;


% %%%%%%%%%%%%%%%%%%%%%
% sparse_remove_columns
% %%%%%%%%%%%%%%%%%%%%%
% cf. LINALG remove_columns

put('sparse_remove_columns, 'psopfn, 'sparse_remove_columns);
put('sparse_remove_columns, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_remove_columns u; % (mtrx, columns)
   % MTRX should be a sparse or dense matrix.  COLUMNS should be a
   % sequence of integers, integer lists or integer intervals
   % representing column indices, i.e. col_1, col_2, ..., col_n, or
   % {col_1, col_2, ..., col_n} or col_1 .. col_n, which means all
   % column indices in the interval from col_1 to col_n inclusive
   % (where `..' is the REDUCE interval operator).  Negative indices
   % are allowed, and count from the right.  In an interval `a .. b',
   % if b < a then the interval is expanded as `a, a-1, a-2, ..., b'.
   % Return a sparse matrix copy of MTRX without the specified
   % columns.
   u and sparse!-process!-mtrx!&cols(u, function sparse!-remove!-columns);

symbolic procedure sparse!-remove!-columns(mtrx, columns);
   % MTRX is a sparse matrix canonical form.
   % COLUMNS is a list of column indices.
   % Return a copy of MTRX without the specified columns,
   % with the columns correctly re-indexed.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      integer newcol;           % initialised to 0
      % Convert columns to a selection list:
      columns := for col := 1 : caddr mtrx join
         if not member(col, columns) then {col};
      % Convert columns to an alist of elements of the form
      % (old_col_ind . new_col_ind):
      columns := for each oldcol in columns collect
         (oldcol . (newcol := newcol + 1));
      maphash(car mtrx,
         (lambda(key, value);
          begin scalar el;
             if (el := assoc(cdr key, columns)) then
                puthash(car key . cdr el, hash, value);
          end));
      return {hash, cadr mtrx, newcol};
   end;

endmodule;

end;
