module sparsedet;          % Determinant and trace of a sparse matrix.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-06 15:50:30 franc>
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

% Using expansion by minors.
% No support for Bareiss algorithm at present!

symbolic procedure simpsparse!-det u;
   % Return the determinant of a sparse matrix, cf. det.
   sparse!-detq sparse!-matsm carx(u, 'sparse_det);

symbolic procedure sparse!-detq u;
   % Top level determinant function.
   % U is a sparse matrix canonical form (<hash> <m> <n>).
   begin scalar len := cadr u;          % number of rows <m>
      if caddr u neq len then rederr "Non square sparse matrix";
      return if len = 1 then
         gethash(1 . 1, car u) or 0
      else sparse!-detq1(car u, len, 1, for j := 1:len collect j)
   end;

symbolic procedure sparse!-detq1(hash, len, i, jlist);
   % HASH contains elements of a sparse square matrix of order LEN.
   % The elements are assumed to be standard quotients.
   % Return the determinant of the matrix.
   % Algorithm is recursive expansion by minors of first row.
   % I is the current "first" row index, initially 1, finally LEN.
   % JLIST is a list of column indices for the current sub-matrix,
   % initially {1,...,LEN}.
   begin scalar i1, result, neg;
      % Base case: last (or single) row, single element.
      if i = len then return gethash(i . car jlist, hash);
      i1 := i + 1;
      result := nil ./ 1;               % zero standard quotient
      for each j in jlist do
      begin scalar el := gethash(i.j, hash);
         if el then <<
            if neg then el := negsq el;
            result := addsq(result, multsq(el,
               % determinant of matrix excluding row i and column j:
               sparse!-detq1(hash, len, i1, delete(j,jlist))));
         >>;
         neg := not neg;
      end;
      return result
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
