module sparsematrix;   % Header for sparse matrices using hash tables.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-04-10 16:47:36 franc>
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

% The representation of a sparse matrix is
%   (sparse-mat <hash> <m> <n> . <name>),
% where <hash> is a hash table, <m> is the maximum row index, <n> is
% the maximum column index, and <name> is either the name of the
% sparse matrix (an identifier) to be used by the print routine or nil
% if it has no name.

% The rtype of a sparse matrix is sparse-matrix.

symbolic inline procedure mk!-sparse!-matrix!-hash;
   mkhash(10, 1);

symbolic procedure sparsematrix u;
   % Declare list U as sparse matrices (represented as hash tables).
   % cf. matrix.
   begin scalar x, y;
        for each j in u do
           if atom j then if null (x := gettype j)
                            then put(j,'rtype,'sparse!-matrix)
                           else if x eq 'sparse!-matrix
                            then <<lprim {x,j,"redefined"};
                                   put(j,'rtype,'sparse!-matrix)>>
                           else typerr({x,j},"sparse matrix")
            else if not idp car j then errpri2(j,'hold)
            else if not (x := gettype car j) or x eq 'sparse!-matrix
             then <<if length j neq 3 then typerr(j,'sparse!-matrix);
                    x := reval_without_mod cadr j;
                    if not fixp x or x<=0 then typerr(x,"positive integer");
                    y := reval_without_mod caddr j;
                    if not fixp y or y<=0 then typerr(y,"positive integer");
                    put(car j, 'rtype, 'sparse!-matrix);
                    put(car j, 'avalue, {'sparse!-matrix,
                       {'sparse!-mat, mk!-sparse!-matrix!-hash(), x, y}})>>
            else typerr({x,car j},"sparse matrix")
   end;

rlistat '(sparsematrix);

put('sparse!-mat, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure quotesparse!-matrix u; 'sparse!-matrix;

flag('(sparse!-mat), 'sparse!-matflg);

flag('(sparse!-mat), 'noncommuting);

put('sparse!-matrix, 'fn, 'sparse!-matflg);

put('sparse!-matrix, 'tag, 'sparse!-mat);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions access / length interface
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symbolic inline procedure sparse!-matdims u;
   % Return dimensions (m n) of sparse matrix u.
   % Assume u = (sparse-mat <hash> <m> <n> . <name>).
   {caddr u, cadddr u};

put('sparse!-matrix, 'lengthfn, 'sparse!-matlength);

symbolic procedure sparse!-matlength u;
   % Return dimensions {m,n} of sparse matrix u.
   % cf. matlength.
   if not eqcar(u,'sparse!-mat) then
      rerror(sparse!-matrix,2,{"Sparse matrix",u,"not set"})
   else 'list . sparse!-matdims u;

% %%%%%%%%%%%%%%
% Element access
% %%%%%%%%%%%%%%

symbolic procedure access!-sparse!-matelem u;
   % Access an element of a sparse matrix u = (id i j ...).
   % Return (hash i j).
   begin scalar x, i, j, dims;
      if length u neq 3 then typerr(u,"sparse matrix element");
      x := get(car u,'avalue);
      if null x or not(car x eq 'sparse!-matrix) then
         typerr(car u,"sparse matrix")
      else if not eqcar(x := cadr x,'sparse!-mat) then
         rerror(sparse!-matrix,1,{"Sparse matrix",car u,"not set"});
      i := reval_without_mod cadr u;
      if not fixp i or i<=0 then typerr(i,"positive integer");
      dims := sparse!-matdims x;          % dims = (m n)
      if i > car dims then
         rerror(sparse!-matrix,23,{"Sparse matrix row number",i,"out of range"});
      j := reval_without_mod caddr u;
      if not fixp j or j<=0 then typerr(j,"positive integer");
      if j > cadr dims then
         rerror(sparse!-matrix,24,{"Sparse matrix column number",j,"out of range"});
      return {cadr x, i, j}
   end;

put('sparse!-matrix, 'getelemfn, 'get!-sparse!-matelem);

symbolic procedure get!-sparse!-matelem u;
   % Return an element of a sparse matrix u = (id i j).
   % cf. getmatelem.
   (gethash(cdr x, car x) or 0) where x = access!-sparse!-matelem u;

put('sparse!-matrix, 'setelemfn, 'set!-sparse!-matelem);

symbolic procedure set!-sparse!-matelem(u,v);
   % Assign v to an element of a sparse matrix u = (id i j)
   % and return v, cf. setmatelem.
   puthash(cdr x, car x, v) where x = access!-sparse!-matelem u;

% %%%%%%%%
% Printing
% %%%%%%%%

put('sparse!-mat, 'prifn, 'sparse!-matpri);

symbolic procedure sparse!-matpri u;
   % Print a sparse matrix u = (sparse!-mat <hash> <m> <n> . <name>)
   % If no (null) name then display name as "?".
   begin scalar alist, name;
      alist := hashcontents cadr u;
      if null alist then return write "Empty matrix";
      % Each element has the form ((i j) value).
      % Sort matrix elements by row index and then by column index:
      alist := sort(alist,
         lambda(x,y);
      caar x < caar y or
         (caar x = caar y and cadar x < cadar y));
      name := cddddr u;
      name := if name then concat2(id2string name, "(") else "?(";
      for each el in alist do eval formwrite(
         'write . {name, caar el, ",", cadar el, ") := ", mkquote cdr el},
         nil, 'algebraic)
   end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random sparse matrices (for testing)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

operator randomsparsematrix;

symbolic procedure randomsparsematrix(s, m, n);
   % s must be an identifier.  Generate an m*n sparse matrix s
   % containing 10 random positive integers.
   begin scalar i, j;
      if not idp s then rederr({s, "invalid as identifier"});
      sparsematrix s(m,n);
      for count := 1:10 do <<
         i := random(m) + 1;
         j := random(n) + 1;
         s(i,j) := random(1000);
      >>;
   end;

endmodule;

in "sparsematsm.red";
in "sparsedet.red";

end;
