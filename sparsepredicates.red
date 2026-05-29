module sparsepredicates;                % Sparse matrix predicates

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-29 17:39:17 franc>
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

% Useful for testing; this should be included in "matrix/matrix.red"!

flag('(conj repart impart), 'matmapfn);

% sparse_matrix_p, cf. LINALG matrixp
% sparse_square_matrix_p, cf. LINALG squarep
% sparse_symmetric_matrix_p, cf. LINALG symmetricp
% sparse_skew_symmetric_matrix_p
% sparse_hermitian_matrix_p
% sparse_skew_hermitian_matrix_p
% sparse_identity_matrix_p
% sparse_orthogonal_matrix_p
% sparse_unitary_matrix_p

% Currently, these predicates are EXACT and do not allow for numerical
% error in floating-point matrices!

symbolic operator sparse_matrix_p, sparse_square_matrix_p,
   sparse_symmetric_matrix_p, sparse_skew_symmetric_matrix_p,
   sparse_hermitian_matrix_p, sparse_skew_hermitian_matrix_p,
   sparse_identity_matrix_p,
   sparse_orthogonal_matrix_p, sparse_unitary_matrix_p;

flag('(sparse_matrix_p sparse_matrix_p
   sparse_symmetric_matrix_p sparse_skew_symmetric_matrix_p
      sparse_hermitian_matrix_p sparse_skew_hermitian_matrix_p
         sparse_identity_matrix_p
            sparse_orthogonal_matrix_p sparse_unitary_matrix_p),
   'boolean);

symbolic procedure sparse_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form), nil otherwise.
   eqcar(u, 'sparse!-mat);

symbolic procedure sparse_square_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is square,
   % nil otherwise.
   eqcar(u, 'sparse!-mat) and caddr u = cadddr u;

% The following functions would benefit from a version of maphash that
% can stop early!  This could be provided via the Common Lisp macro
% with-hash-table-iterator.

symbolic procedure sparse_symmetric_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) symmetric, nil otherwise.
   eqcar(u, 'sparse!-mat) and cadr(u := cdr u) = caddr u and
      begin scalar hash := car u, result := t;
         maphash(hash,
            (lambda(key,val1);
            result and                  % efficiency hack!
            begin scalar i := car key, j := cdr key, val2;
               if i < j then
                  result := (val2 := gethash(j.i, hash)) and val1 = val2;
            end));
         return result;
      end;

symbolic procedure sparse_skew_symmetric_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) skew-symmetric, nil otherwise.
   eqcar(u, 'sparse!-mat) and cadr(u := cdr u) = caddr u and
      begin scalar hash := car u, result := t;
         maphash(hash,
            (lambda(key,val1);
            result and                  % efficiency hack!
            begin scalar i := car key, j := cdr key, val2;
               if i = j then
                  result := val1 = 0
               else if i < j then
                  result := (val2 := gethash(j.i, hash)) and
                  reval {'plus, val1, val2} = 0;
            end));
         return result;
      end;

symbolic procedure sparse_hermitian_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) Hermitian, nil otherwise.
   eqcar(u, 'sparse!-mat) and cadr(u := cdr u) = caddr u and
      begin scalar hash := car u, result := t;
         maphash(hash,
            (lambda(key,val1);
            result and                  % efficiency hack!
            begin scalar i := car key, j := cdr key, val2;
               if i = j then
                  result := reval {'impart, val1} = 0
               else if i < j then
                  result := (val2 := gethash(j.i, hash)) and
                  reval {'difference, val1, {'conj, val2}} = 0;
            end));
         return result;
      end;

symbolic procedure sparse_skew_hermitian_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) skew-Hermitian, nil otherwise.
   eqcar(u, 'sparse!-mat) and cadr(u := cdr u) = caddr u and
      begin scalar hash := car u, result := t;
         maphash(hash,
            (lambda(key,val1);
            result and                  % efficiency hack!
            begin scalar i := car key, j := cdr key, val2;
               if i = j then
                  result := reval {'repart, val1} = 0
               else if i < j then
                  result := (val2 := gethash(j.i, hash)) and
                  reval {'plus, val1, {'conj, val2}} = 0;
            end));
         return result;
      end;

symbolic procedure sparse_identity_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) an identity matrix, nil otherwise.
   eqcar(u, 'sparse!-mat) and caddr u = cadddr u and
      % More efficient to process algebraic form here!
      sparse!-identity!-p sparse!-matsm u;

symbolic procedure sparse!-identity!-p u;
   % U is a sparse matrix in canonical form with SQ elements that may
   % not yet be fully simplified.  Return t if it is an identity
   % matrix, nil otherwise.
   begin scalar result := t;
      maphash(car u,
         (lambda(key,val);
         result and                     % efficiency hack!
            (result := if car key = cdr key then subs2!* val = (1 ./ 1)
            else numr subs2!* val eq nil)));
      !*sub2 := nil;                    % since all substitutions done
      return result;
   end;

symbolic procedure sparse_orthogonal_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) orthogonal, nil otherwise.  A matrix A is orthogonal if
   % A*A^T = A^T*A = I, where ^T denote transpose.
   eqcar(u, 'sparse!-mat) and caddr u = cadddr u and
   begin scalar v;
      u := sparse!-matsm u; % canonical form, SQ elements
      v := sparse!-tp1 u;   % transpose as canonical form, SQ elements
      u := sparse!-multm(u,v);  % A*A^T as canonical form, SQ elements
      return sparse!-identity!-p u;;
   end;

symbolic inline procedure sparse!-conjsq u;
   % See simpconj in "poly/compopr.red".
   multsq(cmpx_conjsf numr u, invsq cmpx_conjsf denr u);

symbolic procedure sparse_unitary_matrix_p u;
   % Return t if U is a sparse matrix (algebraic form) that is (square
   % and) unitary, nil otherwise.  A matrix A is unitary if A*A^H =
   % A^H*A = I, where ^H denotes conjugate transpose.
   eqcar(u, 'sparse!-mat) and caddr u = cadddr u and
   begin scalar v, result := t, hash;
      u := sparse!-matsm u; % canonical form, SQ elements
      v := sparse!-tp1 u;   % transpose as canonical form, SQ elements
      hash := car v;
      maphash(hash,         % conjugate as canonical form, SQ elements
         (lambda(key,val);
         puthash(key, hash, sparse!-conjsq val)));
      u := sparse!-multm(u,v);  % A*A^H as canonical form, SQ elements
      return sparse!-identity!-p u;;
   end;

endmodule;

end;
