module sparselinalg;    % Useful linalg operations for sparse matrices

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-23 18:19:27 franc>
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

% Some useful facilities modelled on LINALG, the REDUCE Linear Algebra
% Package, by Matt Rebbeck.

% sparse_matrix_augment(mat_1, mat_2, ..., mat_n) or
% sparse_matrix_augment({mat_1, mat_2, ..., mat_n}) where mat_i are
% either dense or sparse matrices, all having the same number of rows.
% The result is the matrices adjoined horizontally in the order specified.

put('sparse_matrix_augment, 'psopfn, 'sparse_matrix_augment);

symbolic procedure sparse_matrix_augment u;
   % Adjoin horizontally multiple compatible dense or sparse matrices.
   % Any lists among the arguments are flattened before processing.
   begin scalar m;
      u := for each el in revlis u join
         if eqcar(el, 'list) then cdr el else {el};
      m := cadr lengthreval {car u};
      for each el in cdr u do if cadr lengthreval {el} neq m then
         rederr
            "sparse_matrix_augment arguments must have the same row dimensions";
      % Build a list of sparse matrix canonical forms:
      u := for each el in u collect
         sparse!-matsm if eqcar(el, 'mat) then sparsify el else el;
      return sparse!-matsm!*1 sparse!-matrix!-augment(u, m);
   end;

symbolic procedure sparse!-matrix!-augment(u, m);
   % Adjoin horizontally a list U of sparse matrix canonical forms all
   % with row dim M, and return a sparse matrix canonical form.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      integer n;               % col dim of result so far, initially 0
      for each el in u do <<
         maphash(car el,
            (lambda(key, value);
            puthash(car key . (n + cdr key), hash, value)));
         n := n + caddr el;
      >>;
      return {hash, m, n};
   end;

endmodule;

end;
