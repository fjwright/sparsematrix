module sparserank;                % Sparse matrix rank, cofactor, etc.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-06-23 15:58:16 franc>
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
fluid '(hash!* i!* j!*);
#endif

% %%%%
% Rank
% %%%%

put('rank, 'psopfn, 'generic!-rank!-eval); % updates "matrix/rank.red"

symbolic procedure generic!-rank!-eval u;
   % Return the rank of a generic, i.e. dense or sparse, matrix
   % expression U.
   generic!-matfn(function rank!-eval,
      function sparse!-rank!-eval, {u}, getrtype car u);

% The rank code is based on "matrix/rank.red" by Eberhard Schruefer.

put('sparse_rank, 'psopfn, 'sparse!-rank!-eval);

symbolic procedure sparse!-rank!-eval u;
   if cdr u then rerror(sparse!-matrix, 17, "Wrong number of arguments")
   else if getrtype (u := car u) eq 'sparse!-matrix
   then sparse!-rank!-matrix sparse!-matsm u
   else typerr(u, "sparse matrix");

symbolic procedure sparse!-rank!-matrix u;
   % U = (<hash> <n_rows> <n_cols>) is a sparse matrix canonical form.
   % Return the row rank, i.e. the number of linearly independent rows.
   begin scalar x, y, z := 1; integer n;
      scalar hash := car u, el,
         n_rows := cadr u, n_cols := caddr u;
      %% for each row
      for i := 1 : n_rows do <<
         y := 1;
         %% y := LCM of denominators of elements of row.
         for j := 1 : n_cols do
            if (el := gethash(i.j, hash)) then
               y := lcm(y, denr el);
         x := nil;
         %% x := sum over nonzero elements of the row of (index *
         %% element * LCM of denominators).
         for j := 1 : n_cols do
            if (el := gethash(i.j, hash)) then
               x := list j .* multf(numr el,quotf(y,denr el)) .+ x;
         if y := c!:extmult(x,z)
         then <<z := y; n := n + 1>>
      >>;
      return n
   end;

% %%%%%%%%%%%%%%
% Cofactors, etc
% %%%%%%%%%%%%%%

symbolic operator sparse_submatrix;

put('sparse_submatrix, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_submatrix(u, i, j);
   % Return the submatrix of sparse matrix u excluding row i and
   % column j.  Sparse matrices are represented as tagged algebraic
   % forms.
   if not eqcar(u, 'sparse!-mat) then typerr(u, "sparse matrix")
   else if not fixp i or i <= 0 then typerr(i, "positive integer")
   else if not fixp j or j <= 0 then typerr(j, "positive integer")
   else if i > caddr u then
      rerror(sparse!-matrix, 23, {"Sparse matrix row number",i,"out of range"})
   else if j > cadddr u then
      rerror(sparse!-matrix, 24, {"Sparse matrix column number",j,"out of range"})
   else 'sparse!-mat . sparse!-submatrix(cdr u, i, j);

symbolic procedure sparse!-submatrix(u, i!*, j!*);
   % Return the submatrix of sparse matrix u excluding row i and
   % column j.  Sparse matrices are represented as (<hash> <m> <n>).
   begin scalar hash!* := mk!-sparse!-matrix!-hash();
      maphash(function
         (lambda(key, value);
          begin scalar ii := car key, jj := cdr key;
             if ii < i!* then <<
                if jj < j!* then
                   puthash(key, hash!*, value)
                else if jj > j!* then
                   puthash(ii.(jj-1), hash!*, value)
             >> else if ii > i!* then <<
                if jj < j!* then
                   puthash((ii-1).jj, hash!*, value)
                else if jj > j!* then
                   puthash((ii-1).(jj-1), hash!*, value)
             >>;
          end),
         car u);
      return {hash!*, cadr u - 1, caddr u - 1}
   end;

% The following cofactor code is based partly on "matrix/cofactor.red"
% by Alan Barnes.

put('cofactor, 'simpfn, 'generic!-simpcofactor); % updates "matrix/cofactor.red"
flag('(cofactor), 'immediate);
flag('(cofactor), 'matfn);

symbolic procedure generic!-simpcofactor u;
   % Return the cofactor of a generic, i.e. dense or sparse, matrix
   % expression U.
   generic!-matfn(function simpcofactor, function sparse!-simpcofactor,
      {u}, getrtype car u);

put ('sparse_cofactor, 'simpfn, 'sparse!-simpcofactor);
flag('(sparse_cofactor), 'immediate);

symbolic procedure sparse!-simpcofactor u;
   sparse!-cofactorq(reval car u,
      ieval cadr u, ieval carx(cddr u,'sparse_cofactor));

symbolic procedure sparse!-cofactorq(u, i, j);
   % Return the cofactor of the element in row I and column J of the
   % sparse matrix tagged algebraic form U as a standard quotient.
   % (sparse_submatrix checks its arguments and sparse!-detq checks
   % its argument is square.)
   <<
      u := sparse!-detq sparse!-matsm sparse_submatrix(u, i, j);
      if oddp(i + j) then negsq u else u
   >>;

% The following procedure is not currently used:

symbolic procedure sparse!-generateident n;
   % Return sparse matrix canonical form of identity matrix of order N.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      for i := 1 : n do puthash(i.i, hash, 1 ./ 1);
      return {hash, n, n}
   end;

endmodule;

end;
