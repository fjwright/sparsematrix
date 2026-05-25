module sparselinalg;    % Useful linalg operations for sparse matrices

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-25 18:04:19 franc>
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


% %%%%%%%%%%%%%%%%%%%%%
% sparse_select_columns
% %%%%%%%%%%%%%%%%%%%%%
% cf. LINALG augment_columns

put('sparse_select_columns, 'psopfn, 'sparse_select_columns);
put('sparse_select_columns, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_select_columns u; % (mtrx, columns)
   % MTRX should be a sparse or dense matrix.

   % COLUMNS should be a sequence, list or range of column indices,
   % i.e. col_1, col_2, ..., col_n, or {col_1, col_2, ..., col_n} or
   % col_1 .. col_n, which means all column indices in the interval
   % from col_1 to col_n inclusive (where `..' is the REDUCE interval
   % operator).  Return a sparse matrix copy of MTRX containing only
   % the specified columns.  (Column indices out of range are
   % ignored.)
   begin scalar mtrx, columns, range;
      u := revlis u;
      % Process matrix:
      mtrx := car u;
      if eqcar(mtrx, 'mat) then
         mtrx := sparsify mtrx
      else if not eqcar(mtrx, 'sparse!-mat) then
         typerr(el, "matrix");
      mtrx := sparse!-matsm mtrx;
      % Process columns:
      u := cdr u;
      if eqcar(car u, 'list) then
         columns := cdar u              % (col_1, col_2, ..., col_n)
      else if eqcar(car u, '!*interval!*) then
         range := cdar u                % (start finish)
      else columns := u;                % (col_1, col_2, ..., col_n)
      for each el in columns or range do
         if not fixp el then typerr(el, "matrix column index");
      if range then
         columns := for col := car range : cadr range collect col;
       return sparse!-matsm!*1 sparse!-select!-columns(mtrx, columns);
   end;

symbolic procedure sparse!-select!-columns(mtrx, columns);
   % MTRX is a sparse matrix canonical form.
   % COLUMNS is a list of column indices.
   % Return a copy of MTRX containing only the specified columns,
   % with the columns correctly re-indexed.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      integer j;
      % Convert columns to an alist of elements of the form
      % (old_col_ind . new_col_ind):
      columns := for each col in columns collect
         (col . (j := j + 1));
      maphash(car mtrx,
         (lambda(key, value);
          begin scalar el;
             if (el := assoc(cdr key, columns)) then
                puthash(car key . cdr el, hash, value);
          end));
      return {hash, cadr mtrx, length columns};
   end;

endmodule;

end;
