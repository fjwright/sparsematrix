module sparsepredicates;                % Sparse matrix predicates

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-28 15:42:31 franc>
% Created: April 2026

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

% TO DO:
% sparse_matrix_p, cf. LINALG matrixp
% sparse_square_matrix_p, cf. LINALG squarep
% sparse_symmetric_matrix_p, cf. LINALG symmetricp
% sparse_antisymmetric_matrix_p
% sparse_hermitean_matrix_p
% sparse_antihermitean_matrix_p

% Currently, the symmetry-related predicates test for EXACT symmetry
% and do not allow for numerical error in floating-point matrices!

symbolic operator sparse_matrix_p, sparse_square_matrix_p,
   sparse_symmetric_matrix_p, sparse_antisymmetric_matrix_p;

flag('(sparse_matrix_p sparse_matrix_p
   sparse_symmetric_matrix_p sparse_antisymmetric_matrix_p),
   'boolean);

symbolic procedure sparse_matrix_p u;
   % Return t if u is a sparse matrix (algebraic form); nil otherwise.
   eqcar(u, 'sparse!-mat);

symbolic procedure sparse_square_matrix_p u;
   % Return t if u is a sparse matrix (algebraic form) and it is
   % square; nil otherwise.
   eqcar(u, 'sparse!-mat) and caddr u = cadddr u;

symbolic procedure sparse_symmetric_matrix_p u;
   % Return t if u is a sparse matrix (algebraic form) and it is
   % symmetric; nil otherwise.
   eqcar(u, 'sparse!-mat) and cadr(u := cdr u) = caddr u and
      begin scalar hash := car u, result := t;
         maphash(hash,
            (lambda(key,val1);
            result and                  % efficiency hack!
            begin scalar i := car key, j := cdr key, val2;
               if i < j and (val2 := gethash(j.i, hash)) and
                  val1 neq val2 then
                     result := nil;
            end));
         return result;
      end;

symbolic procedure sparse_antisymmetric_matrix_p u;
   % Return t if u is a sparse matrix (algebraic form) and it is
   % anti-symmetric; nil otherwise.
   eqcar(u, 'sparse!-mat) and cadr(u := cdr u) = caddr u and
      begin scalar hash := car u, result := t;
         maphash(hash,
            (lambda(key,val1);
            result and                  % efficiency hack!
            begin scalar i := car key, j := cdr key, val2;
               if (i = j and val1 neq 0) or
                  (i < j and (val2 := gethash(j.i, hash)) and
                     val1 neq reval {'minus, val2}) then
                        result := nil;
            end));
         return result;
      end;

endmodule;

end;
