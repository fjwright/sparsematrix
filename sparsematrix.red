module sparsematrix;   % Header for sparse matrices using hash tables.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-06-24 15:29:11 franc>
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
% where <hash> is a hash table, <m> is the maximum row index (row
% dimension), <n> is the maximum column index (column dimension), and
% <name> is either the name of the sparse matrix (an identifier) to be
% used by the print routine or nil if it has no name.

% Matrix elements are stored in the hash table under the key
%   (<i> . <j>),
% where <i> is the row index and <j> is the column index.

% The rtype of a sparse matrix is sparse-matrix.

#if (not (memq 'common!-lisp lispsystem!*))
fluid '(fn!* newhash!* f!* u!*);
#endif

% %%%%%%%%%%%%%%%%%
% Utility functions
% %%%%%%%%%%%%%%%%%

% 1000 hash table entries accommodates a 500*500 sparse matrix with
% nonzero diagonal and 500 other nonzero elements.  Also, the REDUCE
% simplifier uses hash-tables with 1000 elements initially (see
% "alg/simp.red")

symbolic inline procedure mk!-sparse!-matrix!-hash;
#if (memq 'psl lispsystem!*)
   mkhash(10, 'equal); % OK in both PSL and CSL, but PSL displays contents!
#else
   mkhash(1000, 3);                     % PSL only accepts 0 & 3 as arg 2!
#endif

% Proposed new Standard Lisp functions, implemented in
% "sl-on-cl.lisp".  The versions below provide a fallback if they are
% not available.  Note that the Standard Lisp function hashcontents
% returns a list of pairs of the form (key . value).

#if (not (getd 'hash!-table!-p))
% Provided in Common Lisp and CSL but not PSL.
symbolic inline procedure hash!-table!-p u;
   % This implementation is not reliable, but currently I only need to
   % distinguish a sparse matrix canonical form from a standard
   % quotient by applying this to the car.
   atom u and not (idp u or numberp u);
#endif

#if (not (getd 'maphash))
% Provided in Common Lisp and CSL but not PSL.
symbolic procedure maphash(fn!*, hash);
   % Iterate over all entries in the hash-table HASH and return nil.
   % For each entry, the function FN is called with two arguments --
   % the key and the value of that entry.
   % mapc(hashcontents hash,
   %    function(lambda el; apply2(fn!*, car el, cdr el)));
   % The above definition doesn't work well in PSL, so...
   for each el in hashcontents hash do
      apply2(fn, car el, cdr el);
#endif

#if (not (getd 'hash!-table!-count))
% Provided in Common Lisp but not CSL or PSL.
symbolic inline procedure hash!-table!-count hash;
   % Return the number of entries in the hash-table HASH.
   length hashcontents hash;
#endif

#if (not (getd 'copyhash))
% Provided in SL-on-CL only.
symbolic procedure copyhash hash;
   % Copy each element of hash table HASH to a new hash table and
   % return the latter.
   begin scalar newhash!* := mk!-sparse!-matrix!-hash();
      maphash(
         function(lambda(key, value); puthash(key, newhash!*, value)),
         hash);
      return newhash!*;
   end;
#endif

symbolic procedure maphash!-new(fn!*, hash);
   % Iterate over all entries in the hash-table HASH and return a new
   % hash-table.  For each entry in HASH, the function FN is called
   % with two arguments -- oldkey, oldval -- and should return a pair
   % (newkey . newval).  Oldkey is the key used to look up an entry
   % with value oldval in hash-table HASH, and newkey is the key used
   % to save newval in the new hash-table.
   begin scalar newhash!* := mk!-sparse!-matrix!-hash();
      maphash(function
         (lambda(oldkey, oldval);
         puthash(car new, newhash!*, cdr new) where
            new = apply2(fn!*, oldkey, oldval)),
         hash);
      return newhash!*;
   end;

symbolic procedure maphash!-new!-values(fn!*, hash);
   % Iterate over all entries in the hash-table HASH and return a new
   % hash-table.  For each entry in HASH, the function FN is called
   % with the entry value as its single argument and its return value
   % is used as the entry value in the new hash-table.
   begin scalar newhash!* := mk!-sparse!-matrix!-hash();
      maphash(function
         (lambda(key, val);
         puthash(key, newhash!*, apply1(fn!*, val))),
         hash);
      return newhash!*;
   end;

symbolic inline procedure puthash!-nzsq(key, hash, value);
   % Avoid putting a zero SQ entry into a sparse matrix hash table.
   % VALUE is a SQ; if it is nonzero then insert it into hash table
   % HASH, otherwise remove the entry in hash table HASH.
   if numr value then puthash(key, hash, value)
   else remhash(key, hash);

% This should probably be in the main REDUCE source:

symbolic operator mat2list;

symbolic procedure mat2list m;
   % Convert matrix M to a list of lists
   begin scalar mm := reval m;
      if not eqcar(mm, 'mat) then typerr(m, "matrix");
      return 'list . for each row in cdr mm collect 'list . row;
   end;

% Oddp is defined in SL-on-CL and in "rtools/general.red" in the
% development system but not in the last snapshot release, revision
% 7327), build date 08-Mar-2026, so ...
#if (null (getd 'oddp))
symbolic inline procedure oddp n;  not evenp n;
#endif


% %%%%%%%%%%%
% Declaration
% %%%%%%%%%%%

symbolic procedure sparse_matrix u;
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

rlistat '(sparse_matrix);

put('sparse!-mat, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure quotesparse!-matrix u; 'sparse!-matrix;

% flag('(sparse!-mat), 'sparse!-matflg);

flag('(sparse!-mat), 'noncommuting);

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
   if not eqcar(u, 'sparse!-mat) then
      rerror(sparse!-matrix, 2, {"Sparse matrix",u,"not set"})
   else 'list . sparse!-matdims u;


% %%%%%%%%%%%%%%
% Element access
% %%%%%%%%%%%%%%

symbolic procedure access!-sparse!-matelem u;
   % Access an element of a sparse matrix, where U = (id i j).
   % Return ((i . j) . hash).
   begin scalar x, i, j, dims;
      if length u neq 3 then typerr(u, "sparse matrix element");
      x := get(car u, 'avalue);
      if null x or not(car x eq 'sparse!-matrix) then
         typerr(car u, "sparse matrix")
      else if not eqcar(x := cadr x, 'sparse!-mat) then
         rerror(sparse!-matrix, 1, {"Sparse matrix",car u,"not set"});
      i := reval_without_mod cadr u;
      if not fixp i or i <= 0 then typerr(i, "positive integer");
      dims := sparse!-matdims x;        % dims = (m n)
      if i > car dims then
         rerror(sparse!-matrix, 23, {"Sparse matrix row number",i,"out of range"});
      j := reval_without_mod caddr u;
      if not fixp j or j <= 0 then typerr(j, "positive integer");
      if j > cadr dims then
         rerror(sparse!-matrix, 24, {"Sparse matrix column number",j,"out of range"});
      return ((i . j) . cadr x)
   end;

put('sparse!-matrix, 'getelemfn, 'get!-sparse!-matelem);

symbolic procedure get!-sparse!-matelem u;
   % Return an element of a sparse matrix u = (id i j).
   % cf. getmatelem.
   (gethash(car x, cdr x) or 0)
      where x = access!-sparse!-matelem u;

put('sparse!-matrix, 'setelemfn, 'set!-sparse!-matelem);

symbolic procedure set!-sparse!-matelem(u,v);
   % Assign v to an element of a sparse matrix u = (id i j)
   % and return v, cf. setmatelem.
   (if zerop v then << remhash(car x, cdr x); 0 >>
   else puthash(car x, cdr x, v))
      where x = access!-sparse!-matelem u;


% %%%%%%%
% Mapping
% %%%%%%%

% Explicitly map an operator over the elements of a sparse matrix:

put('sparse!-mat, 'mapfn, 'map!-sparse!-mat);

symbolic procedure map!-sparse!-mat(f!*,o);
   {'sparse!-mat,
      maphash!-new(function
         (lambda(key, value); (key . apply1(f!*, value))),
         cadr o),
         caddr o, cadddr o};

% Automatically map an operator over the elements of a sparse matrix:

put('sparse!-matrix, 'aggregatefn, 'sparse!-matrixmap);
put('sparse!-matrix, 'fn, 'matflg);

flag('(sparse_det sparse_trace sparse_cofactor), 'matfn);

symbolic procedure sparse!-matrixmap(u!*, v);
   % U = (<function> <sparse matrix> <other args>).
   % Apply <function> to each element of <sparse matrix>, cf. matrixmap.
   % The sparse matrix is input and output in tagged algebraic form.
   if flagp(car u!*, 'matmapfn)
   then sparse!-matsm!*1
      ({maphash!-new(function
         (lambda(key, value);
         (key . simp!*(car u!* . mk!*sq value . cddr u!*))),
         car sparse!-matsm sm),
         caddr sm, cadddr sm} where sm = cadr u!*)
   else if flagp(car u!*, 'matfn) then reval2(u!*, v)
   else typerr(car u!*, "sparse matrix operator");


% %%%%%%%%
% Printing
% %%%%%%%%

share sparse_matrix_dense_print_colmax;
sparse_matrix_dense_print_colmax := 10;
switch sparse_matrix_dense_print = on;

% The following code is used by assgnpri.

% Needed for special printing of assignments of sparse matrices:
flag('(sparse!-matrix), 'sprifn);
put('sparse!-mat, 'assgnpri, 'sparse!-assgnpri);

symbolic procedure sparse!-assgnpri uvw;
   % Called by assgnpri to print a sparse matrix
   % or an assignment of the form
   %   <variable> := <sparse matrix>
   % UVW = (u v w), where
   % U = (sparse!-mat <hash> <m> <n> . <name>)
   % V = (<variable>) or null if not an assignment
   % W = only
   begin scalar u := car uvw, v := cadr uvw;
      % Display as a dense matrix if feasible, mainly for testing
      % with small sparse matrices:
      if !*sparse_matrix_dense_print and
         caddr u <= sparse_matrix_dense_print_colmax then
            return assgnpri(densify u, v, 'only);
      if v then
         u := 'sparse!-mat . cadr u . caddr u . cadddr u . car v;
      sparse!-matpri u;
   end;

% Needed to print sparse matrices inside other structures, such as
% lists:
put('sparse!-mat, 'prifn, 'sparse!-matpri);

symbolic procedure sparse!-matpri u;
   % Display as a dense matrix if feasible, mainly for testing
   % with small sparse matrices:
   if !*sparse_matrix_dense_print and
      caddr u <= sparse_matrix_dense_print_colmax
   then matpri densify u else
      % Print a sparse matrix u = (sparse!-mat <hash> <m> <n> . <name>)
      % If no (null) name then display name as "?".
      begin scalar alist := hashcontents car (u := cdr u),
            msg := {cadr u, "#times;", caddr u,
            "sparse matrix #mdash;"};
         if null alist then return
            lprim append(msg, {"no nonzero elements"});
         lprim append(msg, {length alist, "nonzero elements:"});
         % Each alist element has the form ((i . j) . value).
         % Sort by row index and then by column index:
         alist := sort(alist, function
            (lambda(x,y);
            caar x < caar y or
               (caar x = caar y and cdar x < cdar y)));
         for each el in alist do
            assgnpri(cdr el, {{cdddr u or '!?, caar el, cdar el}}, 'only);
      end;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate sparse random matrices (for testing)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symbolic procedure sparse_random_matrix u;
   % U must evaluate to a list of elements of the form (s m n).
   % S must be an identifier.  Generate an M*N sparse matrix S
   % containing (M+N)/2 random positive integers.
   for each v in u do
   begin scalar m, n, i, j, hash;
      % Essentially, sparse_matrix s(m, n):
      if (m := gettype car v) and not (m eq 'sparse!-matrix)
      then typerr({m, car v}, "sparse matrix");
      if length v neq 3 then typerr(v, 'sparse!-matrix);
      m := reval_without_mod cadr v;
      if not fixp m or m <= 0 then typerr(m, "positive integer");
      n := reval_without_mod caddr v;
      if not fixp n or n <= 0 then typerr(n, "positive integer");
      put(car v, 'rtype, 'sparse!-matrix);
      hash := mk!-sparse!-matrix!-hash();
      put(car v, 'avalue, {'sparse!-matrix,
         {'sparse!-mat, hash, m, n}});
      % Now assign some elements of S:
      for count := 1 : fix((m+n)/2) do <<
         i := random(m) + 1;
         j := random(n) + 1;
         puthash(i.j, hash, random(1000));
      >>;
   end;

rlistat '(sparse_random_matrix);

% %%%%%%%%%%%%%%%%%%%
% Density of a matrix
% %%%%%%%%%%%%%%%%%%%

symbolic operator matrix_density;

symbolic procedure matrix_density u;
   % U must evaluate to a dense or sparse matrix.  Return its density,
   % namely the proportion of nonzero elements, as a percentage
   % truncated to the nearest integer.
   begin scalar type := getrtype u;
      integer nz;                    % count of nonzero elements
      return
         if type eq 'matrix then <<
            u := matsm u;
            for each row in u do for each el in row do
               if numr el then nz := nz + 1;
            quotient(nz * 100, length u * length car u)
         >> else if type eq 'sparse!-matrix then <<
            u := sparse!-matsm u;
            nz := hash!-table!-count car u;
            quotient(nz * 100, cadr u * caddr u)
         >>;
   end;

endmodule;

end;
