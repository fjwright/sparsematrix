module sparsematsm;               % Simplification of sparse matrices.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-12 17:50:21 franc>
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

load_package matrix;                    % needed for densify

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
   % U is an arbitrary sparse matrix expression in algebraic form.
   % Return a sparse matrix expression in tagged algebraic form.
   sparse!-matsm!*1 sparse!-matsm u;

symbolic procedure sparse!-matsm!*1 u;
   % Assume u evaluates to a sparse matrix canonical form
   %   (<hash> <m> <n> . <name>)
   % where the matrix elements are in SQ form.
   % Return a COPY in which each element is converted to an ALGEBRAIC
   % EXPRESSION in the form
   %   (sparse!-mat <hash> <m> <n> . <name>).
   % *** TEMPORARY HACK TO CHECK SIMPLER FACILITIES! ***
   'sparse!-mat . map!-sparse!-matrix(u, function !*q2a, t);

symbolic procedure sparse!-matsm u;
   % Simplify an arbitrary sparse matrix expression U in algebraic
   % form and return a sparse matrix canonical form
   %   (<hash> <m> <n> . <name>)
   % where <name> is U if U is an identifier or nil, and the matrix
   % elements are STANDARD QUOTIENT FORMS.

   % nssimp returns a sparse matrix expression as a list representing
   % a sum, where each summand is a list of the form (c m1 m2 ...)
   % representing a product, c is a scalar in standard quotient form,
   % and the mi are tagged algebraic sparse matrix forms.
   begin scalar x, y, name := idp u and u;
      for each j in nssimp(u, 'sparse!-matrix) do <<
         % multiply out each product...
         y := sparse!-multsm(car j, sparse!-matsm1(cdr j, name));
         % and add them...
         x := if null x then y else sparse!-addm(x,y)
      >>;
      return x
   end;

symbolic procedure sparse!-matsm1(u, name);
   % U is a sparse matrix symbol product, i.e. a general sparse matrix
   % expression; NAME is an identifier (for a single sparse matrix) or
   % nil.

   % Return nil or a sparse matrix canonical form
   %   (<hash> <m> <n> . <name>)
   % where <name> is either an identifier or nil, and the hash table
   % elements are STANDARD QUOTIENT FORMS.
   begin scalar x,y,z; integer n;
   a:
      if null u then return z
      else if eqcar(car u, '!*div) then go to d % inverse
      else if atom car u then go to er          % not set
      else if caar u eq 'sparse!-mat then go to c1 % tagged alg form
      else if flagp(caar u, 'matmapfn) and cdar u  % map applied
         and getrtype cadar u eq 'sparse!-matrix
      then x := sparse!-matsm sparse!-matrixmap(car u,nil)
      else if (x := get(caar u, 'psopfn)) % psopfn function call
      then << x := lispapply(x, list cdar u);
         if eqcar(x,'sparse!-mat) then x := sparse!-matsm x >>
      else <<x := lispapply(caar u,cdar u); % other function call
         if eqcar(x,'sparse!-mat) then x := sparse!-matsm x>>;
   b: % Multiplication (scalar or matrix):
      z := if null z then x
      else if null cdr z and null cdar z then sparse!-multsm(caar z,x)
      else sparse!-multm(x,z);
   c: % Loop through elements of matrix expression u:
      u := cdr u;
      go to a;
      % End of main loop, followed by special cases.
   c1:                       % tagged algebraic form
      %   car u = (sparse!-mat <hash> <m> <n>)
      % Return a sparse matrix canonical form
      %   (<hash> <m> <n> . <name>):
      x := map!-sparse!-matrix(cdar u, function xsimp, name);
      go to b;
   d: % Inverse
      y := sparse!-matsm cadar u;
      if (n := cadr y) neq caddr y
      then rerror(sparse!-matrix,4,"Non square sparse matrix")
      else if (z and n neq cadr z)
      then rerror(sparse!-matrix,5,"Sparse matrix mismatch")
      else if cddar u then go to h
      else if n = 1 then go to e; % 1*1 matrix
      x := subfg!*;
      subfg!* := nil;
      if null z then z := apply1(get('sparse!-mat,'inversefn),y)
      else if null(x := get('sparse!-mat,'lnrsolvefn))
      then z := multm(apply1(get('sparse!-mat,'inversefn),y),z)
      else z := apply2(get('sparse!-mat,'lnrsolvefn),y,z);
      subfg!* := x;
      % Make sure there are no power substitutions:
      z := map!-sparse!-matrix(z, lambda value;
                               <<!*sub2 := t; subs2 value>>);
      go to c;
   e: % y is 1*1 matrix, cf. y = ((el))
      y := gethash(1 . 1, car y); % nil or value of single element as SQ
      if null(y and numr y) then
         % y is 1*1 zero matrix, cf. y = ((nil ./ 1))
         % cf. mat(())^-1 or  mat((0))^-1 or M/mat(()) or M/mat((0))
         rerror(sparse!-matrix,6,"Zero divisor");
      y := revpr y;                     % now y is a scalar
      z := if null z then list list y else sparse!-multsm(y,z);
      go to c;
   h:
      if null z then z := sparse!-generateident n;
      go to c;
   er:
      rerror(sparse!-matrix,7,list("Sparse matrix",car u,"not set"))
   end;

% %%%%%%%%
% Addition
% %%%%%%%%

% This procedure would benefit from an efficient hash table copy
% function.

symbolic procedure sparse!-addm(u,v);
   % Return the sum of two sparse matrix canonical forms U and V as a
   % sparse matrix canonical form.
   % U & V have the form (<hash> <m> <n>).
   if not(cadr u = cadr v and caddr u = caddr v) then
      rerror(sparse!-matrix,8,"Sparse matrix mismatch")
   else begin scalar hash := mk!-sparse!-matrix!-hash();
      % Copy each nonzero element of sparse matrix U to hash:
      maphash(car u,
         (lambda(key, u_val);
         puthash(key, hash, u_val)));
      % Add or in each nonzero element of sparse matrix V to hash:
      maphash(car v,
         (lambda(key, v_val);
          begin scalar u_val;
             puthash(key, hash,
                if (u_val := gethash(key, hash)) then
                   addsq(u_val, v_val)
                else v_val)
          end));
      return {hash, cadr u, caddr u}
   end;

% %%%%%%%%%
% Transpose
% %%%%%%%%%

% This code currently works but seems a bit convoluted, in that it
% appears to end up calling sparse!-matsm multiple time!

symbolic procedure sparse_tp u; sparse!-tp1 sparse!-matsm u;

put('sparse_tp, 'rtypefn, 'getrtypecar); % declares algebraic operator
% flag('(sparse_tp), 'sparse!-matflg);

symbolic procedure sparse!-tp1 u;
   % Return the transpose of the sparse matrix canonical form U =
   % (<hash> <m> <n>) as a new sparse matrix canonical form.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      maphash(car u,
         (lambda(key, value);
         puthash(cdr key . car key, hash, value)));
      return {hash, caddr u, cadr u}
   end;

% %%%%%%%%%%%%%%
% Multiplication
% %%%%%%%%%%%%%%

symbolic procedure sparse!-multm(u,v);
   % Return the product of two sparse matrix canonical forms U and V
   % as a new sparse matrix canonical form.  Assume U and V are
   % conformable, i.e. caddr u = cadr v.
   begin scalar hash := mk!-sparse!-matrix!-hash();
      maphash(car u,
         (lambda(u_key, u_value);
          begin scalar i := car u_key, k := cdr u_key;
             maphash(car v,
                (lambda(v_key, v_value);
                if car v_key = k then
                   % The product of this pair of matrix elements is a
                   % summand of the scalar product forming the i,j element
                   % of the product matrix.
                   begin scalar j := cdr v_key,
                         scalprod := gethash(i.j, hash),
                         prod := multsq(u_value, v_value);
                      puthash(i.j, hash,
                         if scalprod then addsq(scalprod, prod) else prod);
                   end));
          end));
      return {hash, cadr u, caddr v}
   end;

symbolic procedure sparse!-multsm(u,v);
   % Return the product of standard quotient U and sparse matrix
   % canonical form V as a new sparse matrix canonical form.
   if u = (1 ./ 1) then v else
      map!-sparse!-matrix(v,
         % Ordering of multsq arguments to preserve the ordering of
         % noncom scalars in matrix elements!
         (lambda value; multsq(value, u)));

% %%%%%%%%%%%%
% Substitution
% %%%%%%%%%%%%

put('sparse!-matrix, 'subfn, 'sparse!-matsub);

symbolic procedure sparse!-matsub(u,v);
   % V is a tagged algebraic sparse matrix form;
   % U is a substitution equation represented as a dotted pair.
   % Return a new tagged algebraic sparse matrix form with
   % substitution U applied to every element, cf. matsub.
   'sparse!-mat .
      map!-sparse!-matrix(cdr v, (lambda value; subeval1(u, value)));

% %%%%%%%%%%
% Conversion
% %%%%%%%%%%

symbolic operator sparsify;

put('sparsify, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparsify u;
   % Convert dense matrix algebraic form U to a sparse matrix
   % algebraic form.
   sparsify!-matrix matsm u;

symbolic procedure sparsify!-matrix u;
   % Convert dense matrix canonical form U to a sparse matrix
   % algebraic form (cf. matsm!*1).
   begin scalar hash := mk!-sparse!-matrix!-hash();
      integer i, j;
      for each row in u do <<
         i := i + 1;  j := 0;
         for each el in row do <<
            j := j + 1;
            if numr el then             % nonzero standard quotient
               puthash(i.j, hash, !*q2a el);
         >>;
      >>;
      return {'sparse!-mat, hash, i, j}
   end;

symbolic operator densify;

put('densify, 'rtypefn, 'quotematrix);

symbolic procedure densify u;
   % Convert sparse matrix algebraic form U to a dense matrix
   % algebraic form.
   densify!-matrix sparse!-matsm u;

symbolic procedure densify!-matrix u;
   % Convert sparse matrix canonical form U to a dense matrix
   % algebraic form (cf. matsm!*1).
   begin scalar hash := car u, m := cadr u, n := caddr u, el;
      return 'mat .
         for i := 1 : m collect
            for j := 1 : n collect
               if el := gethash(i.j, hash) then !*q2a el else 0
   end;

endmodule;

end;
