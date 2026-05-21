module sparsedet;          % Determinant and trace of a sparse matrix.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-21 18:00:44 franc>
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
   % Return determinant as a SQ.
   begin scalar m := cadr u, hash, neg, d := 1 ./ 1;
      if caddr u neq m then rederr "Non square sparse matrix";
      if m = 1 then return gethash(1 . 1, car u) or (nil ./ 1);
      hash := car u;
      % Reduce hash (destructively) to row echelon form and return
      % 'singular if the matrix is singular; otherwise return non-nil
      % if the sign of the determinant has been changed (by an odd
      % number of row swaps):
      neg := sparse!-echelon(hash, m, m, t);
      if neg eq 'singular then return (nil ./ 1);
      for i := 1 : m do d := multsq(d, gethash(i.i, hash));
      return if neg then negsq d else d;
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
