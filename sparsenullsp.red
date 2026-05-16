module sparsenullsp;       % Compute the nullspace of a sparse matrix.

% Author: Herbert Melenk <melenk@sc.zib-berlin.de>.
% Revised for sparse matrices represented as hash tables by FJW.
% Time-stamp: <2026-05-16 17:14:11 franc>

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

% This file is a reworking of "matrix/nullsp.red" to use hash tables
% to represent sparse matrices.

% Algorithm: Rational Gaussian elimination with standard quotients.

put('sparse_nullspace, 'psopfn, 'sparse!-nullspace!-eval);

% TO DO: Make sparse_nullspace accept an optional second argument that
% determines whether it returns a list of sparse column vectors or a
% list of lists (as at present).

symbolic procedure sparse!-nullspace!-eval u;
   % Interface for the nullspace calculation.
   % U is a sparse matrix in tagged algebraic form.  Return the basis
   % of the kernel of U as a list of lists.
   begin scalar v := reval car u;
      if not eqcar(v, 'sparse!-mat) then
         typerr ("sparse!-matrix", u);
      v := sparse!-nullspace!-alg cdr v;
      return 'list . for each vect in v collect
         'list . vect;
   end;

symbolic procedure sparse!-nullspace!-alg u;
   % U is a sparse matrix, encoded as (<hash> <m> <n>) where the
   % matrix elements are algebraic expressions.
   % Result is the basis of the kernel of U as a list of lists.
   begin scalar mp, vars, rvars, res, oldorder;
      scalar hash := car u;
      integer n := caddr u;             % # columns
      vars := for i := 1 : n collect gensym();
      rvars := reverse vars;
      oldorder := setkorder rvars;
      % Build a list of SQs, one for each row, where each SQ is a sum
      % of terms of the form v*<matrix element> and each v is a new
      % variable (gensym) that is distinct for each column.
      mp := for i := 1 : cadr u collect % for each row
         begin scalar el,
               v := vars,
               r := nil ./ 1;           % SQ zero
               for j := 1 : n do <<     % for each row element
                  el := gethash(i.j, hash);
                  if el then r := addsq(r, simp {'times, car v, el});
                  v := cdr v;
               >>;
               return r;
         end;
      res := nullspace!-elim(mp, rvars);
      setkorder oldorder;
      return reverse for each q in res collect
         for each x in vars collect
            cdr atsoc(x, q);
   end;

% symbolic procedure nullspace!-elim(m,vars) is defined in
% "matrix/nullsp.red".

endmodule;

end;
