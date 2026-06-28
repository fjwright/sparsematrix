% sparse_random_matrix
% cf. LINALG sparse_random_matrix

% A sparse analogue of the LINALG random_matrix operator but with more
% flexibility without using switches.

% sparse_random_matrix(m, n, types), where types can be one or more
% of:

% Number type:

%   By default, a random integer r is generated in the range lo <= r <
%   hi, where lo = -limit, hi = limit, and limit = 1000.

%   The first positive integer resets limit to that value.
%   Alternatively, the first interval with integer endpoints of the
%   form lo .. hi, with lo < hi, resets lo and hi to those values.  To
%   obtain positive matrix elements use lo = 1; it doesn't makes sense
%   to use lo = 0 in a sparse matrix because 0 elements will be
%   essentially ignored!

%   rational; default is integer
%   complex; default is real

%   A random rational number is a quotient of a random integer num in
%   the range lo <= num < hi and a random integer den in the range 0 <
%   den < hi.  A random complex number has real and imaginary parts
%   that are random integer or rational numbers.

% Matrix density:

%   density = <positive rational number>; an integer is interpreted as
%   a percentage.  The default density assigns values to a number of
%   elements equal to the mean matrix dimension.

% Square matrix types -- ignored if the matrix is not square:

%   invertible; assigns 1 to any diagonal elements that would
%   otherwise be zero.

%   At most one of the following mutually-exclusive types:

%   diagonal, band(number), upper, lower,
%   symmetric, anti/skew_symmetric,
%   hermitian, anti/skew_hermitian

% Most repeated or invalid types are silently ignored.

remprop('sparse_random_matrix, 'stat);  % TEMPORARY!

put('sparse_random_matrix, 'psopfn, 'sparse_random_matrix);
put('sparse_random_matrix, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_random_matrix u; % (m n types)
   % M and N must evaluate to positive integers.  TYPES is an optional
   % sequence of type identifiers.  Return an M*N sparse matrix
   % containing (M+N)/2 random positive integers.
   if length u < 2 then
      rederr "Wrong number of arguments to sparse_random_matrix"
   else
   begin scalar m, n, hash, types,
         lo := -1000, hi := 1000, density, maxcount,
         !*diagonal, !*upper, !*lower,
         !*symmetric, !*anti_symmetric,
         !*hermitian, !*anti_hermitian,
         realvalue, value;
      m := reval_without_mod car u;
      if not fixp m or m <= 0 then typerr(m, "positive integer");
      n := reval_without_mod cadr u;
      if not fixp n or n <= 0 then typerr(n, "positive integer");
      hash := mk!-sparse!-matrix!-hash();
      types := cddr u;                % list of types
      % Set number range:
      begin scalar tps := types, tp, lo1, hi1;
         while tps do
            if fixp (tp := car tps) and tp > 0
            then << lo := -tp;  hi := tp;  tps := nil >>
            else if eqcar(tp, '!*interval!*) and fixp(lo1 := cadr tp)
               and fixp(hi1 := caddr tp) and lo1 < hi1
            then << lo := lo1;  hi := hi1;  tps := nil >>
            else tps := cdr tps;
      end;
      % Set density:
      begin scalar tps := types, tp;
         while tps do
            if eqcar(tp := car tps, 'equal) and cadr tp eq 'density
            then <<
               density := caddr tp;
               if fixp density then        % percentage
                  density := {'quotient, density, 100}
               else if not eqcar(density, 'quotient) then
                  typerr(density, "sparse random matrix density");
               tps := nil
            >> else tps := cdr tps;
      end;
      maxcount := if density then
         numr simp {'fix, {'times, density, m, n}} or 0
      else (m+n)/2;                     % integer division
      % Set matrix type:
      if m = n then
         if 'diagonal memq types then !*diagonal := t
         else if 'upper memq types then !*upper := t
         else if 'lower memq types then !*lower := t
         else if 'symmetric memq types then !*symmetric := t
         else if 'anti_symmetric memq types
            or 'skew_symmetric memq types then
               !*anti_symmetric := t
         else if 'hermitian memq types then !*hermitian := t
         else if 'anti_hermitian memq types
            or 'skew_hermitian memq types then
               !*anti_hermitian := t;
      % Set element value function:
      realvalue := if 'rational memq types then
         (lambda(); {'quotient, num!-value(lo, hi), den!-value hi})
      else
         (lambda(); num!-value(lo, hi));
      value := if 'complex memq types
         or !*hermitian or !*anti_hermitian then
            (lambda (); {'plus, apply(realvalue, nil),
               {'times, 'i, apply(realvalue, nil)}})
      else realvalue;
      % Assign random values to random elements:
      for count := 1 : maxcount do
         begin scalar i, j, val;
            i := random(m) + 1;
            j := if !*diagonal then i else random(n) + 1;
            if (!*upper and j < i) or (!*lower and j > i) then return;
            if !*anti_symmetric and i = j then return;
            % Filter out 0 values:
            repeat val := apply(value, nil)
               until numr simp val;
            puthash(i.j, hash, val);
            if !*symmetric then puthash(j.i, hash, val)
            else if !*anti_symmetric then puthash(j.i, hash, -val)
            else if !*hermitian then
               if i = j then puthash(j.i, hash, {'repart, val})
               else puthash(j.i, hash, {'conj, val})
            else if !*anti_hermitian then
               if i = j then puthash(j.i, hash, {'times, 'i, {'impart, val}})
               else puthash(j.i, hash, {'minus, {'conj, val}});
         end;
      if m = n and not !*anti_symmetric and 'invertible memq types then
         for i := 1 : m do
            if not gethash(i.i, hash) then
               if !*anti_hermitian then puthash(i.i, hash, 'i)
               else puthash(i.i, hash, 1);
      return {'sparse!-mat, hash, m, n}
   end;

symbolic procedure num!-value(lo, hi);
   % Return a random integer r such that lo <= r < hi.
   lo + random(hi - lo);

symbolic procedure den!-value limit;
   % Return a random positive integer less than limit.
   1 + random(limit - 1);

end;
