module sparseechelon;    % Reduce a sparse matrix to row echelon form.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-24 15:36:14 franc>
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

put('sparse_echelon, 'rtypefn, 'getrtypecar); % declares algebraic operator

symbolic procedure sparse_echelon u;
   % Return the sparse matrix in row echelon form.
   % U is a tagged algebraic form.
   % Return a sparse matrix canonical form
   begin scalar hash, m, n;
      u := sparse!-matsm u;
      hash := car u;
      m := cadr u;
      n := caddr u;
      % Reduce hash (destructively) to row echelon form:
      sparse!-echelon(hash, m, n, nil);
      mathprint densify sparse!-matsm!*1 {hash, m, n};
      return sparse!-matsm!*1 {hash, m, n};
   end;

% Using reduction to row echelon form (Gaussian elimination).
% No support for Bareiss algorithm at present!

% The following row reduction code is based on
% https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode

% sparse!-echelon is faster than sparse!-echelon1 and so
% sparse!-echelon is used for determinant computation.

symbolic procedure sparse!-echelon(hash, m, n, det);
   % HASH contains the elements of a sparse M*N matrix (A).
   % The elements are assumed to be standard quotients.
   % On return the elements in HASH are in row echelon form.
   % Return non-nil if odd # row swaps, nil otherwise.
   % If DET is non-nil then return 'singular as soon as a singular
   % determinant is detected.
   begin scalar
      h := 1,                           % initial pivot row
      k := 1,                           % initial pivot column
      neg;                              % true if odd # row swaps
      while h <= m and k <= n do
      begin scalar i_piv := h, pivot;
         % Find the first (nonzero) pivot below row h in column k:
         while i_piv <= m and null (pivot := gethash(i_piv.k, hash)) do
            i_piv := i_piv + 1;
         if i_piv > m then <<
            if det then return 'singular;
            % No pivot in this column, pass to next column
            k := k + 1
         >> else <<
            if i_piv > h then <<
               % Swap rows h and i_piv:
               for j := k : n do sparse!-el!-swap(hash, h.j, i_piv.j);
               neg := not neg;
            >>;
            % Do for all rows below pivot:
            for i := h + 1 : m do
               begin scalar f := gethash(i.k, hash); % A[i, k]
                  if null f then return; % row already in echelon form
                  f := negsq quotsq(f, pivot); % - A[i, k] / A[h, k]
                  % Fill lower part of pivot column with zeros:
                  remhash(i.k, hash);   % A[i, k] := 0
                  % Do for all remaining elements in this row:
                  for j := k + 1 : n do
                     % A[i, j] := A[i, j] - A[h, j] * f
                     begin scalar change := gethash(h.j, hash);
                        if change then <<
                           change := multsq(change, f);
                           sparse!-add!-to!-el(hash, i.j, change);
                        >>;
                     end;
               end;
            % Increase pivot row and column:
            h := h + 1;
            k := k + 1;
         >>;
      end;
      return neg;
   end;

symbolic procedure sparse!-el!-swap(hash, i1_j1, i2_j2);
   % Swap elements with keys I1_J1 and I2_J2 in hash table HASH.
   begin scalar
      val1 := gethash(i1_j1, hash),
      val2 := gethash(i2_j2, hash);
      if val1 then <<
         puthash(i2_j2, hash, val1);
         if val2 then
            puthash(i1_j1, hash, val2)
         else
            remhash(i1_j1, hash);
      >> else if val2 then <<
         puthash(i1_j1, hash, val2);
         remhash(i2_j2, hash);
      >>;
   end;

symbolic procedure sparse!-add!-to!-el(hash, i_j, value);
   % Add VALUE to element with key I_J in hash table HASH.
   % Do not save a zero element.  Assume values are SQs.
   begin scalar old_val := gethash(i_j, hash);
      if old_val then value := addsq(old_val, value);
      puthash!-nzsq(i_j, hash, value);
   end;

put('sparse_canonical, 'rtypefn, 'getrtypecar); % declares algebraic operator

symbolic procedure sparse_canonical u;
   % Return the sparse matrix in row canonical form.
   % U is a tagged algebraic form.
   % Return a sparse matrix canonical form
   begin scalar hash, m, n, neg;
      u := sparse!-matsm u;
      hash := car u;
      m := cadr u;
      n := caddr u;
      % Reduce hash (destructively) to row canonical form:
      neg := sparse!-echelon(hash, m, n, t);
      mathprint densify sparse!-matsm!*1 {hash, m, n};
      if neg eq 'singular then rederr("Singular leading submatrix");
      sparse!-canonical(hash, m, n);
      mathprint densify sparse!-matsm!*1 {hash, m, n};
      return sparse!-matsm!*1 {hash, m, n};
   end;

symbolic procedure sparse!-canonical(hash, m, n);
   % HASH contains the elements of a sparse M*N (augmented) matrix (A)
   % in row echelon form.  The elements are assumed to be standard
   % quotients.  On return the elements in HASH are in row canonical
   % form.  Assume that the leading M*M submatrix is non-singular.
   for i := m step -1 until 1 do <<
      % Re-scale row(i) so that A[i,i] = 1:
      begin scalar f := invsq gethash(i.i, hash); % f = 1/A[i,i]
         % row(i) := row(i) / A[i,i]
         puthash(i.i, hash, 1 ./ 1); % A[i,i] := 1
         for j := i+1 : n do
            begin scalar el := gethash(i.j, hash);
               if el then
                  puthash(i.j, hash, multsq(el, f));
            end;
      end;
      % Zero col(i) above the (now unit) pivot A[i,i]:
      for ii := i-1 step -1 until 1 do
         begin scalar f := gethash(ii.i, hash); % f = A[ii,i]
            if f then <<
               % row(ii) := row(ii) - A[ii,i]*row(i) where A[i,i] = 1
               remhash(ii.i, hash);     % A[ii,i] := 0
               f := negsq f;            % f = - A[ii,i]
               for j := i+1 : n do
                  % A[ii,j] := A[ii,j] - A[ii,i]*A[i,j] (if A[i,j] neq 0)
                  begin scalar A_i_j := gethash(i.j, hash), A_ii_j;
                     if A_i_j then <<
                        A_i_j := multsq(f, A_i_j);
                        A_ii_j := gethash(ii.j, hash);
                        puthash!-nzsq(ii.j, hash,
                           if A_ii_j then addsq(A_ii_j, A_i_j) else A_i_j);
                     >>;
                  end;
            >>;
         end;
   >>;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse and linear solve (inverse times matrix)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Functions called by sparse!-matsm1 to support matrix inverse
% arithmetic:

put('sparse!-mat, 'inversefn, 'sparse!-matinverse);

symbolic inline procedure sparse!-matinverse u;
   % Return the inverse of sparse matrix U.
   sparse!-lnrsolve(u, nil);

put('sparse!-mat, 'lnrsolvefn, 'sparse!-lnrsolve);

symbolic procedure sparse!-lnrsolve(u, v);
   % U is a sparse matrix.  Return U^(-1)*V if V is a sparse matrix or
   % U^(-1) if V is nil, in which case it defaults to the identity
   % matrix.  Use reduction of the augmented matrix to row canonical
   % form.  Assume U is m*m and V is m*n, so the product is m*n, and
   % all matrices are represented as sparse matrix canonical forms.
   begin scalar hash := copyhash car u,
         m := cadr u, n, sing, newhash;
      n := if v then <<                 % augment U with V
         maphash(car v,
            (lambda(key, value);
            puthash(car key . (cdr key + m), hash, value)));
         caddr v
      >> else <<                        % augment U with a unit matrix
         for i := 1 : m do puthash(i . (i + m), hash, 1 ./ 1);
         m
      >>;
      % Reduce hash (destructively) to row canonical form:
      sing := sparse!-echelon(hash, m, m + n, t);
      if sing eq 'singular then
         rerror(sparse!-matrix, 13, "Singular sparse matrix");
      sparse!-canonical(hash, m, m + n);
      % Extract the product or inverse matrix:
      newhash := mk!-sparse!-matrix!-hash();
      maphash(hash,
         (lambda(key, value);
         if cdr key > m then
            puthash(car key . (cdr key - m), newhash, value)));
      return {newhash, m, n};
   end;

endmodule;

end;
