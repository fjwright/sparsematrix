module sparserank;                      % Rank of a sparse matrix

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-08 17:02:52 franc>
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

% Based on "matrix/rank.red" by Eberhard Schruefer.

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

endmodule;

end;
