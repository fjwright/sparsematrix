module sparse_random_matrix;
% cf. LINALG sparse_random_matrix

% A sparse analogue of the LINALG random_matrix operator but with more
% flexibility without using switches.

% sparse_random_matrix(m, n, types), where n is optional and types can
% be zero or more of the following.

% Allowed types -- at most one of each of the following
% mutually-exclusive type classes, where the actual input is shown in
% upper case:

% Element type:
%   either numeric: default, or specified by
%     range: LIMIT=N where n is a positive integer or LO .. HI where
%     lo, hi are integers and lo < hi
%     RATIONAL: keyword
%     COMPLEX: keyword
%   or SYMBOL: keyword

% Density:
%   DENSITY=N, where n is a positive integer, rational or float

% Matrix type (only if square) -- at most one of the following
% mutually-exclusive type:
%   DIAGONAL,
%   BAND=N or BAND, where n is a positive integer that defaults to 3,
%   UPPER, LOWER,
%   SYMMETRIC, ANTI/SKEW_SYMMETRIC,
%   HERMITIAN, ANTI/SKEW_HERMITIAN

%   INVERTIBLE

% Number type:

%   By default, a random integer r is generated in the range lo <= r <
%   hi, where lo = -lim, hi = lim, and lim = 1000.

%   To obtain positive matrix elements use lo = 1; it doesn't makes
%   sense to use lo = 0 in a sparse matrix because 0 elements will be
%   essentially ignored!

%   rational; default is integer
%   complex; default is real

%   A random rational number is a quotient of a random integer num in
%   the range lo <= num < hi and a random integer den in the range 0 <
%   den < hi.  A random complex number has real and imaginary parts
%   that are random integer or rational numbers.

% Matrix density:

%   An integer is interpreted as a percentage, whereas a rational or
%   float is interpreted as a fraction, so 100, 1/1 and 1.0 all mean
%   the same.  The default density assigns values to a number of
%   elements equal to the mean matrix dimension.

% Square matrix types -- ignored if the matrix is not square:

%   invertible: assigns 1 to any diagonal elements that would
%   otherwise be zero.

%   band(n): creates a band matrix with n nonzero elements in each row
%   centred about the main diagonal.

%   upper, lower: create upper and lower triangular matrices.

remprop('sparse_random_matrix, 'stat);  % TEMPORARY!

put('sparse_random_matrix, 'rtypefn, 'quotesparse!-matrix);
put('sparse_random_matrix, 'formfn, 'form_sparse_random_matrix);

symbolic procedure form_sparse_random_matrix(u, vars, mode);
   % Allow symmetric as an argument (even though it is a REDUCE
   % keyword).
   begin scalar args := cadr u .
      for each arg in cddr u collect
         if eqcar(arg, 'symmetric) then ''symmetric else arg;
      return form1('eval_sparse_random_matrix . args, vars, mode);
   end;

put('eval_sparse_random_matrix, 'psopfn, 'eval_sparse_random_matrix);

fluid '(diagonal upper lower symm anti_symm herm anti_herm);
global '(sparse_random_matrix_types);
sparse_random_matrix_types := {'(diagonal), '(upper), '(lower),
   '(symmetric . symm), '(hermitian . herm),
   '(anti_symmetric . anti_symm), '(anti_hermitian . anti_herm),
   '(skew_symmetric . anti_symm), '(skew_hermitian . anti_herm)};

symbolic procedure eval_sparse_random_matrix u; % (m n types)
   % M must evaluate to a positive integer.  N is optional and if
   % specified must evaluate to a positive integer; it defaults to the
   % value of M.  TYPES is an optional sequence of type specifiers.
   % Return an M*N sparse matrix containing by default (M+N)/2 random
   % positive integers.
   if null u then
      rederr "Wrong number of arguments to sparse_random_matrix"
   else
   begin scalar m, n, hash, element_type, matrix_type, tp,
         lo := -1000, hi := 1000, density, rational, complex, symbol,
         diagonal, band, upper, lower,
         symm, anti_symm, herm, anti_herm, invertible,
         maxcount, bandspread, realvalue, value;
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
                  typerr(type, "sparse random matrix element range");
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
               rational := t;
               element_type := 'numeric;
            >> else if type eq 'complex then << % COMPLEX element type
               if element_type eq 'symbolic then
                  rederr "Element type already set";
               complex := t;
               element_type := 'numeric;
            >> else if type eq 'symbol then << % SYMBOL element type
               if element_type eq 'numeric then
                  rederr "Element type already set";
               symbol := t;
               element_type := 'symbolic;
            >> else if type eq 'band then << % default BAND matrix type
               if matrix_type then rederr "Matrix type already set";
               band := 3;
               matrix_type := t;
            >> else if tp := assoc(type, sparse_random_matrix_types) then <<
               if m neq n then rederr "Matrix must be square";
               if matrix_type then rederr "Matrix type already set";
               set(if cdr tp then cdr tp else type, t);
               matrix_type := t;
            >> else if type eq 'invertible then <<
               if m neq n then rederr "Matrix must be square";
               if invertible then rederr "Matrix type already set";
               invertible := t;
            >> else typerr(type, "sparse random matrix type");
         >> else typerr(type, "sparse random matrix type");

      maxcount := if density then
         numr simp {'fix, {'times, density, m, n}} or 0
      else (m+n)/2;                     % integer division

      % Set element value function:
      realvalue := if rational then
         (lambda(); {'quotient, num!-value(lo, hi), den!-value hi})
      else
         (lambda(); num!-value(lo, hi));
      value := if complex or herm or anti_herm then
         (lambda (); {'plus, apply(realvalue, nil),
            {'times, 'i, apply(realvalue, nil)}})
      else realvalue;
      % Assign random values to random elements:
      if band then <<
         bandspread := (band-1)/2;
         for i := 1 : m do
            for j := i - bandspread : i + bandspread do
               if 1 <= j and j <= n then
               begin scalar val;
                  % Filter out 0 values:
                  repeat val := apply(value, nil)
                     until numr simp val;
                  puthash(i.j, hash, val);
               end
      >> else for count := 1 : maxcount do
         begin scalar i, j, val;
            i := random(m) + 1;
            j := if diagonal then i else random(n) + 1;
            if (upper and j < i) or (lower and j > i) then return;
            if anti_symm and i = j then return;
            % Filter out 0 values:
            repeat val := apply(value, nil)
               until numr simp val;
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
      if m = n and not anti_symm and invertible then
         for i := 1 : m do
            if not gethash(i.i, hash) then
               if anti_herm then puthash(i.i, hash, 'i)
               else puthash(i.i, hash, 1);
      return {'sparse!-mat, hash, m, n}
   end;

symbolic procedure num!-value(lo, hi);
   % Return a random integer r such that lo <= r < hi.
   lo + random(hi - lo);

symbolic procedure den!-value limit;
   % Return a random positive integer less than limit.
   1 + random(limit - 1);

endmodule;

end;

% TO DO:

% Fully implement symbol type.
% Don't apply density to special matrix types, which are already sparse by definition?
% Dense type that allows random zeros?
