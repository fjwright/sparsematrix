module sparserandmat;                   % cf. LINALG random_matrix

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-07-22 11:59:08 franc>
% Created: June 2026

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

% A sparse-matrix analogue of the LINALG random_matrix operator with
% more flexibility but without using switches.

% sparse_random_matrix(m, n, types) returns a random sparse-matrix,
% where n is optional and types can be zero or more of the following.

% Allowed types -- at most one of each of the following
% mutually-exclusive type classes, where the actual input is shown in
% upper case:

% Element type:
%   either numeric: default, or specified by
%     range: LIMIT=N where n is a positive integer or LO .. HI where
%     lo, hi are integers and lo < hi
%     RATIONAL
%     COMPLEX
%   or SYMBOL

% Density:
%   DENSITY=N, where n is a positive integer, rational or float

% Matrix type (only if square) -- at most one of the following
% mutually-exclusive type:
%   DIAGONAL, UPPER, LOWER,
%   BAND or BAND=N, where n is a positive integer that defaults to 3,
%   SYMMETRIC, ANTI/SKEW_SYMMETRIC,
%   HERMITIAN, ANTI/SKEW_HERMITIAN

%   INVERTIBLE

% Number type:

%   By default, a random integer r is generated in the range lo <= r <
%   hi, where lo = -limit, hi = limit, and limit = 1000.

%   To obtain positive matrix elements use lo = 1; it doesn't make
%   sense to use lo = 0 in a sparse matrix because 0 elements will be
%   essentially ignored!

%   rational; default is integer
%   complex; default is real

%   A random rational number is a quotient of a random integer num in
%   the range lo <= num < hi and a random integer den in the range 0 <
%   den < hi.  A random complex number has real and imaginary parts
%   that are random integer or rational numbers.  With ON ROUNDED, all
%   numerical elements will appear as floats.

% Matrix density:

%   An integer is interpreted as a percentage, whereas a rational or
%   float is interpreted as a fraction, so 100, 1/1 and 1.0 all mean
%   the same.  The default density assigns nonzero values to a number
%   of elements equal to the mean matrix dimension, except that for
%   specified matrix types (not including invertible) the default
%   density is 100% (on the grounds that such a matrix is sparse by
%   definition).

% Square matrix types -- invalid if the matrix is not square:

%   invertible: not allowed for an anti-symmetric matrix, otherwise
%   assigns the appropriate unit element to any diagonal element that
%   would otherwise be zero.

%   band=n: creates a band matrix with n nonzero elements in each row
%   centred (approximately) about the main diagonal.

%   upper, lower: create upper and lower triangular matrices.
%   hermitian, anti/skew_hermitian also set complex element type.

put('sparse_random_matrix, 'rtypefn, 'quotesparse!-matrix);
put('sparse_random_matrix, 'formfn, 'form!-sparse!-random!-matrix);

symbolic procedure form!-sparse!-random!-matrix(u, vars, mode);
   % Allow symmetric as an argument (even though it is a REDUCE
   % keyword).
   begin scalar args := cadr u .
      for each arg in cddr u collect
         if eqcar(arg, 'symmetric) then ''symmetric else arg;
      return form1('sparse!-random!-matrix . args, vars, mode);
   end;

put('sparse!-random!-matrix, 'psopfn, 'sparse!-random!-matrix);

fluid '(rational!* complex!* symbol!*
   diagonal upper lower symm anti_symm herm anti_herm);
global '(sparse!-random!-matrix!-types);
sparse!-random!-matrix!-types := {'(diagonal), '(upper), '(lower),
   '(symmetric . symm), '(hermitian . herm),
   '(anti_symmetric . anti_symm), '(anti_hermitian . anti_herm),
   '(skew_symmetric . anti_symm), '(skew_hermitian . anti_herm)};

symbolic procedure sparse!-random!-matrix u; % (m n types)
   % M must evaluate to a positive integer.  N is optional and if
   % specified must evaluate to a positive integer; it defaults to the
   % value of M.  TYPES is an optional sequence of type specifiers.
   % Return an M*N sparse matrix containing by default (M+N)/2 random
   % positive integers.
   if null u then
      rederr "Wrong number of arguments to sparse_random_matrix"
   else
   begin scalar m, n, hash, element_type, matrix_type, tp,
         lo := -1000, hi := 1000, density, rational!*, complex!*, symbol!*,
         diagonal, band, upper, lower,
         symm, anti_symm, herm, anti_herm, invertible,
         maxcount, bandspread;
      m := reval_without_mod car u;
      if not fixp m or m <= 0 then typerr(m, "positive integer");
      if null(u := cdr u) then n := m else <<
         n := reval_without_mod car u;
         if fixp n and n > 0 then u := cdr u else n := m;
      >>;
      hash := mk!-sparse!-matrix!-hash();

      % Process types:
      for each type in u do
         if eqcar(type, '!*interval!*) then << % INTERVAL element type
            if element_type then
               rederr "Element type already set";
            if not (fixp(lo := cadr type) and fixp(hi := caddr type)
               and lo < hi) then
                  typerr(type, "sparse random matrix element interval");
            element_type := 'numeric;
         >> else if eqcar(type, 'equal) then << % equational type form
            if (tp := cadr type) eq 'limit then << % LIMIT element type
               if element_type then
                  rederr "Element type already set";
               if not (fixp(hi := caddr type) and hi > 0) then
                  typerr(type, "sparse random matrix element limit");
               lo := -hi; element_type := 'numeric;
            >> else if tp eq 'density then << % DENSITY
               density := caddr type;
               if fixp density and      % percentage
                  0 < density and density <= 100 then
                     density := {'quotient, density, 100}
               else if eqcar(density, '!:dn!:) and % non-negative float
                  cddr density < 0 then
                     density := {'quotient, cadr density, 10^(-cddr density)}
               else if not eqcar(density, 'quotient) % fraction
               then typerr(type, "sparse random matrix density");
            >> else if tp eq 'band then << % BAND matrix type
               if matrix_type then
                  rederr "Matrix type already set";
               if not (fixp(band := caddr type) and band > 0) then
                  typerr(type, "sparse random matrix band type");
               matrix_type := t;
            >> else typerr(type, "sparse random matrix type");
         >> else if idp type then <<   % keyword type form
            if type eq 'rational then << % RATIONAL element type
               if element_type eq 'symbolic then
                  rederr "Element type already set";
               rational!* := t;
               element_type := 'numeric;
            >> else if type eq 'complex then << % COMPLEX element type
               if element_type eq 'symbolic then
                  rederr "Element type already set";
               complex!* := t;
               element_type := 'numeric;
            >> else if type eq 'symbol then << % SYMBOL element type
               if element_type eq 'numeric then
                  rederr "Element type already set";
               symbol!* := t;
               element_type := 'symbolic;
            >> else if type eq 'band then << % default BAND matrix type
               if matrix_type then rederr "Matrix type already set";
               band := 3;
               matrix_type := t;
            >> else if tp := assoc(type, sparse!-random!-matrix!-types) then <<
               if m neq n then rederr "Matrix must be square";
               if matrix_type then rederr "Matrix type already set";
               set(if cdr tp then cdr tp else type, t);
               if herm or anti_herm then complex!* := t;
               matrix_type := t;
            >> else if type eq 'invertible then <<
               if m neq n then rederr "Matrix must be square";
               if invertible then rederr "Matrix type already set";
               invertible := t;
            >> else typerr(type, "sparse random matrix type");
         >> else typerr(type, "sparse random matrix type");
      if anti_symm and invertible then
         rederr {"anti/skew_symmetric and invertible together",
            "is an invalid sparse random matrix type combination"};
      if matrix_type and not density then density := '(quotient 1 1);

      maxcount := if density then
         numr simp {'fix, {'times, density, m, n}} or 0
      else (m+n)/2;                     % integer division

      % Assign random values to random elements:
      if band then <<
         bandspread := (band-1)/2;
         for i := 1 : m do
            for j := i - bandspread : i + bandspread do
               if 1 <= j and j <= n then
                  puthash(i.j, hash, srm!-nonzero!-value(lo, hi));
      >> else for count := 1 : maxcount do
         begin scalar i, j, val := srm!-nonzero!-value(lo, hi);
            i := random(m) + 1;
            j := if diagonal then i else random(n) + 1;
            if (upper and j < i) or (lower and j > i) then return;
            if anti_symm and i = j then return;
            puthash(i.j, hash, val);
            if symm then puthash(j.i, hash, val)
            else if anti_symm then puthash(j.i, hash, -val)
            else if herm then
               if i = j then puthash(j.i, hash, {'repart, val})
               else puthash(j.i, hash, {'conj, val})
            else if anti_herm then
               if i = j then puthash(j.i, hash, {'times, 'i, {'impart, val}})
               else puthash(j.i, hash, {'minus, {'conj, val}});
         end;
      if invertible then                % => square and not anti-symm
         for i := 1 : m do
            if not gethash(i.i, hash) then
               if anti_herm then puthash(i.i, hash, 'i)
               else puthash(i.i, hash, 1);
      return {'sparse!-mat, hash, m, n}
   end;

symbolic procedure srm!-integer!-value(lo, hi);
   % Return a random integer r such that lo <= r < hi.
   lo + random(hi - lo);

symbolic procedure srm!-posint!-value limit;
   % Return a random positive integer less than limit.
   1 + random(limit - 1);

symbolic procedure srm!-real!-value(lo, hi);
   % Return a random rational or integer value.
   if rational!* then
      {'quotient, srm!-integer!-value(lo, hi), srm!-posint!-value hi}
   else
      srm!-integer!-value(lo, hi);

symbolic procedure srm!-value(lo, hi);
   % Return a random numerical value.
   if complex!* then
      {'plus, srm!-real!-value(lo, hi),
         {'times, 'i, srm!-real!-value(lo, hi)}}
   else srm!-real!-value(lo, hi);

symbolic procedure srm!-nonzero!-value(lo, hi);
   % Return a symbol or nonzero random numerical value.
   if symbol!* then gensym()
   else begin scalar val;
      % Filter out 0 values:
      repeat val := srm!-value(lo, hi)
         until numr simp val;
      return val;
   end;

endmodule;

end;

% TO DO:

% More testing!
% Dense type that allows random zeros?
