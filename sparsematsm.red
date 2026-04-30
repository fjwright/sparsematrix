module sparsematsm;               % Simplification of sparse matrices.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-04-30 17:57:32 franc>
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

% This file is a reworking of "matrix/matsm.red" to use hash tables to
% represent sparse matrices.

% The canonical form of an <m>*<n> matrix is
%   ((el_11 el_12 ... el_1n)
%    (el_21 el_22 ... el_2n)
%    ...
%    (el_m1 el_m2 ... el_mn))
% where el_ij is the ij matrix element in SQ form.

% The canonical form of an <m>*<n> sparse matrix is
%   (<hash> <m> <n>)
% where <hash> is a hash table and the ij matrix element is stored in
% SQ form in the hash table with key (i j).

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation and simplification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

put('sparse!-matrix, 'evfn, 'sparse!-matsm!*);

symbolic procedure sparse!-matsm!*(u,v);
   % Sparse matrix expression simplification function.
   sparse!-matsm!*1 sparse!-matsm u;

symbolic procedure sparse!-matsm!*1 u;
   % Assume u evaluates to a sparse matrix canonical form
   %   (<hash> <m> <n> . <name>)
   % where the matrix elements are in SQ form.
   % Return a COPY in which each element is converted to an ALGEBRAIC
   % EXPRESSION in the form
   %   (sparse!-mat <hash> <m> <n> . <name>).
   % *** TEMPORARY HACK TO CHECK SIMPLER FACILITIES! ***
   begin scalar hash := mk!-sparse!-matrix!-hash();
      for each el in hashcontents car u do
         puthash(car el, hash, !*q2a cdr el);
      return 'sparse!-mat . hash . cdr u;
   end;

symbolic procedure sparse!-matsm u;
   % Return a sparse matrix canonical form
   %   (<hash> <m> <n> . <name>)
   % where <name> is either an identifier or nil, and the hash table
   % elements are STANDARD QUOTIENT FORMS.

   % If U evaluates to a variable assigned a sparse matrix then return
   %   (<hash> <m> <n> . u).
   % If U evaluates to the form
   %   (sparse!-mat <hash> <m> <n>)
   % then return
   %   (<hash> <m> <n>).
   % If U evaluates to an operator expression of the form (op s) then
   % simplify s to the form (<hash> <m> <n>) and apply op it.
   % *** TEMPORARY HACK TO CHECK SIMPLER FACILITIES! ***
   begin scalar x, name;
      if idp u and (x := get(u, 'avalue))
         and eqcar(x, 'sparse!-matrix)
            and eqcar(x := cadr x, 'sparse!-mat) then <<
               x := cdr x;
               name := u;
            >>
      else if eqcar(u, 'sparse!-mat) then
         x := cdr u
      else return apply(car u, {sparse!-matsm(cadr u)});
      % Convert matrix elements to standard quotients in a NEW hash
      % table:
      return begin scalar hash := mk!-sparse!-matrix!-hash();
         for each el in hashcontents car x do
            puthash(car el, hash, simp cdr el);
         return hash . cadr x . caddr x . name
      end;
   end;

% %%%%%%%%
% Addition
% %%%%%%%%

% This procedure would benefit from an efficient hash table copy
% function.

symbolic procedure sparse!-addm(u,v);
   % Return the sum of two sparse matrix canonical forms U and V as a
   % sparse matrix canonical form.  Return U + 0 as U and 0 + V as V.
   % U & V have the form (<hash> <m> <n>).
   if v = '(((nil . 1))) then u
   else if u = '(((nil . 1))) then v
   else if not(cadr u = cadr v and caddr u = caddr v) then
      rerror(sparse!-matrix,8,"Sparse matrix mismatch")
   else
   begin scalar hash := mk!-sparse!-matrix!-hash(), val;
      % Each element of hashcontents list has the form
      % ((i j) . value).
      for each el in hashcontents car u do
         puthash(car el, hash, cdr el);
      for each el in hashcontents car v do
         puthash(car el, hash,
            if val := gethash(car el, hash) then
               addsq(val, cdr el)
            else cdr el);
      return {hash, cadr u, caddr u}
   end;

% %%%%%%%%%
% Transpose
% %%%%%%%%%

symbolic procedure sparsetp u; sparse!-tp1 sparse!-matsm u;

flag('(sparsetp), 'sparse!-matflg);
put('sparsetp, 'rtypefn, 'getrtypecar); % declares algebraic operator

symbolic procedure sparse!-tp1 u;
   % Return the transpose of the sparse matrix internal form U =
   % (<hash> <m> <n>) as a new sparse matrix internal form.
   begin scalar
      alist := hashcontents(car u),
      % Each alist element has the form ((i j) value).
      newhash := mk!-sparse!-matrix!-hash();
      for each el in alist do           % write transposed element
         puthash({cadar el,caar el}, newhash, cdr el);
      return {newhash, caddr u, cadr u}
   end;

% %%%%%%%%%%
% Conversion
% %%%%%%%%%%

% NB: s := sparsify m; does not display correctly!

symbolic operator sparsify;

symbolic procedure sparsify u;
   % Convert dense matrix algebraic form U to a sparse matrix
   % algebraic form.
   sparse!-matsm!*1 sparsify!-matrix matsm u;

symbolic procedure sparsify!-matrix u;
   % Convert dense matrix canonical form U to a sparse matrix
   % canonical form.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      integer i, j;
      for each row in u do <<
         i := i + 1;  j := 0;
         for each el in row do <<
            j := j + 1;
            if el neq '(((nil . 1))) then puthash({i,j}, hash, el);
         >>;
      >>;
      return {hash, i, j}
   end;

symbolic operator densify;

symbolic procedure densify u;
   % Convert sparse matrix algebraic form U to a dense matrix
   % algebraic form.
   matsm!*1 densify!-matrix sparse!-matsm u;

symbolic procedure densify!-matrix u;
   % Convert sparse matrix canonical form U to a dense matrix
   % canonical form.
   begin scalar hash := car u, m := cadr u, n := caddr u, el;
      return
         for i := 1 : m collect
            for j := 1 : n collect
               if el := gethash({i,j}, hash) then el else '(((nil . 1)))
   end;

endmodule;

end;
