module sparsematrix;   % Header for sparse matrices using hash tables.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-04-04 17:11:10 franc>
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

% This file is a reworking of "matrix/matrix.red" to use hash tables
% to represent sparse matrices.

fluid '(!*sub2 subfg!*);

global '(nxtsym!*);

% The representation of a sparse matrix is
%   (sparsemat <hash> <m> <n>),
% where <hash> is a hash table, <m> is the maximum row index and <n>
% is the maximum column index.

% The rtype of a sparse matrix is sparsematrix.

symbolic procedure sparsematrix u;
   % Declare list U as sparse matrices (represented as hash tables).
   % cf. matrix.
   begin scalar x, y;
        for each j in u do
           if atom j then if null (x := gettype j)
                            then put(j,'rtype,'sparsematrix)
                           else if x eq 'sparsematrix
                            then <<lprim {x,j,"redefined"};
                                   put(j,'rtype,'sparsematrix)>>
                           else typerr({x,j},"sparse matrix")
            else if not idp car j then errpri2(j,'hold)
            else if not (x := gettype car j) or x eq 'sparsematrix
             then <<if length j neq 3 then typerr(j,'sparsematrix);
                    x := reval_without_mod cadr j;
                    if not fixp x or x<=0 then typerr(x,"positive integer");
                    y := reval_without_mod caddr j;
                    if not fixp y or y<=0 then typerr(y,"positive integer");
                    put(car j, 'rtype, 'sparsematrix);
                    put(car j, 'avalue, {'sparsematrix,
                       {'sparsemat, mkhash(10, 1), x, y}})>>
            else typerr({x,car j},"sparse matrix")
   end;

rlistat '(sparsematrix);

put('sparsemat, 'rtypefn, 'quotesparsematrix);

symbolic procedure quotesparsematrix u; 'sparsematrix;

% Dimensions access / length interface

symbolic inline procedure sparsematdims u;
   % Return dimensions (m n) of sparse matrix u.
   % Assume u = (sparsemat hashtable m n).
   cddr u;

put('sparsematrix, 'lengthfn, 'sparsematlength);

symbolic procedure sparsematlength u;
   % Return dimensions {m,n} of sparse matrix u.
   % cf. matlength.
   if not eqcar(u,'sparsemat) then
      rerror(sparsematrix,2,{"Sparse matrix",u,"not set"})
   else 'list . sparsematdims u;

% Element access

symbolic procedure accesssparsematelem u;
   % Access an element of a sparse matrix u = (id i j ...).
   % Return (hash i j).
   begin scalar x, i, j, dims;
      if length u neq 3 then typerr(u,"sparse matrix element");
      x := get(car u,'avalue);
      if null x or not(car x eq 'sparsematrix) then
         typerr(car u,"sparse matrix")
      else if not eqcar(x := cadr x,'sparsemat) then
         rerror(sparsematrix,1,{"Sparse matrix",car u,"not set"});
      i := reval_without_mod cadr u;
      if not fixp i or i<=0 then typerr(i,"positive integer");
      dims := sparsematdims x;          % dims = (m n)
      if i > car dims then
         rerror(sparsematrix,23,{"Sparse matrix row number",i,"out of range"});
      j := reval_without_mod caddr u;
      if not fixp j or j<=0 then typerr(j,"positive integer");
      if j > cadr dims then
         rerror(sparsematrix,24,{"Sparse matrix column number",j,"out of range"});
      return {cadr x, i, j}
   end;

put('sparsematrix, 'getelemfn, 'getsparsematelem);

symbolic procedure getsparsematelem u;
   % Return an element of a sparse matrix u = (id i j).
   % cf. getmatelem.
   (gethash(cdr x, car x) or 0) where x = accesssparsematelem u;

put('sparsematrix, 'setelemfn, 'setsparsematelem);

symbolic procedure setsparsematelem(u,v);
   % Assign v to an element of a sparse matrix u = (id i j)
   % and return v, cf. setmatelem.
   puthash(cdr x, car x, v) where x = accesssparsematelem u;

endmodule;

end;
