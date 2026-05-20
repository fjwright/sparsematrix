module sparsedet;          % Determinant and trace of a sparse matrix.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-20 18:00:53 franc>
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

% This file is a reworking of "matrix/det.red" to use hash tables to
% represent sparse matrices.

% %%%%%%%%%%%
% Determinant
% %%%%%%%%%%%

put('sparse_det, 'simpfn, 'simpsparse!-det);
flag('(sparse_det), 'immediate);

% Using reduction to row echelon form (Gaussian elimination).
% No support for Bareiss algorithm at present!

symbolic procedure simpsparse!-det u;
   % Return the determinant of a sparse matrix, cf. det.
   sparse!-detq sparse!-matsm carx(u, 'sparse_det);

symbolic procedure sparse!-detq u;
   % Top level determinant function.
   % U is a sparse matrix canonical form (<hash> <m> <n>).
   begin scalar m := cadr u, hash, neg, d := 1 ./ 1;
      if caddr u neq m then rederr "Non square sparse matrix";
      if m = 1 then return gethash(1 . 1, car u) or 0;
      hash := car u;
      neg := sparse!-echelon(hash, m, m);
      for i := 1 : m do
         d := multsq(d, gethash(i.i, hash) or (nil ./ 1));
      return if neg then negsq d else d;
   end;

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
      sparse!-echelon(hash, m, n);
      return sparse!-matsm!*1 {hash, m, n};
   end;

% The following row reduction code is based on
% https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode

symbolic procedure sparse!-echelon(hash, m, n);
   % HASH contains the elements of a sparse M*N matrix.
   % The elements are assumed to be standard quotients.
   % On return the elements in HASH are in row echelon form.
   % Return non-nil if odd # row swaps, nil otherwise.
   begin scalar
      h := 1,                           % initial pivot row
      k := 1,                           % initial pivot column
      neg;                              % true if odd # row swaps
      while h <= m and k <= n do
      begin scalar i_piv := h, pivot;
         % Find the first (nonzero) pivot below row h in column k:
         while i_piv <= m and null (pivot := gethash(i_piv.k, hash)) do
            i_piv := i_piv + 1;
         if i_piv > m then
            % No pivot in this column, pass to next column
            k := k + 1
         else <<
            if i_piv > h then <<
               % Swap rows h and i_piv:
               for j := k : n do sparse!-el!-swap(hash, h.j, i_piv.j);
               neg := not neg;
            >>;
            % Do for all rows below pivot:
            for i := h + 1 : m do
               begin scalar f := negsq quotsq(gethash(i.k, hash), pivot);
                  % Fill with zeros the lower part of pivot column:
                  remhash(i.k, hash);
                  % Do for all remaining elements in this row:
                  for j := k + 1 : n do
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
      % mathprint densify sparse!-matsm!*1 {hash, m, n};
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
   % Assume all values are SQs.
   begin scalar old_val := gethash(i_j, hash);
      puthash(i_j, hash,
         if old_val then addsq(old_val, value) else value);
   end;

% %%%%%
% Trace
% %%%%%

put('sparse_trace, 'simpfn, 'simpsparse!-trace);

symbolic procedure simpsparse!-trace u;
   % Return the trace of a sparse matrix, cf. trace.
   begin scalar m, hash, el, z;
      u := sparse!-matsm carx(u, 'sparse_trace); % (<hash> <m> <n>)
      if (m := cadr u) neq caddr u then
         rederr "Non square sparse matrix";
      hash := car u;
      % The matrix elements are standard quotients.
      z := nil ./ 1;                    % zero standard quotient
      for i := 1 : m do
         if (el := gethash(i.i, hash)) then z := addsq(el, z);
      return z
   end;

endmodule;

end;
