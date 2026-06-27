% sparse_random_matrix
% cf. LINALG sparse_random_matrix

% A sparse analogue of the LINALG random_matrix operator but with more
% flexibility without using switches.

% sparse_random_matrix(m, n, types), ... where types can be one or
% more of:

% Number type:

%   rational; default is integer
%   complex; default is real

%   The first positive integer specifies a limit, default 1000,
%   meaning that random integers r are generated in the range 0 <= r <
%   limit.  A random rational number is a quotient of such random
%   integers.  A random complex number has real and imaginary parts
%   that are such random integer or rational numbers.

% Matrix type:

%   diagonal, band(number), upper, lower,
%   symmetric, anti_symmetric/skew_symmetric,
%   hermitian, anti_hermitian/skew_hermitian
%   invertible

% Repeated or invalid types are silently ignored.

remprop('sparse_random_matrix, 'stat);  % TEMPORARY!

put('sparse_random_matrix, 'psopfn, 'sparse_random_matrix);
put('sparse_random_matrix, 'rtypefn, 'quotesparse!-matrix);

symbolic procedure sparse_random_matrix u; % (m n types)
   % M and N must evaluate to positive integers.  TYPES is an optional
   % sequence of type identifiers.  Return an M*N sparse matrix
   % containing (M+N)/2 random positive integers.
   if length u < 2 then
      rederr "Wrong number of arguments to sparse_random_matrix"
   else begin scalar m, n, hash, types, i, j, limit := 1000,
         realvalue, value;
      m := reval_without_mod car u;
      if not fixp m or m <= 0 then typerr(m, "positive integer");
      n := reval_without_mod cadr u;
      if not fixp n or n <= 0 then typerr(n, "positive integer");
      hash := mk!-sparse!-matrix!-hash();
      types := cddr u;                % list of types

      begin scalar tps := types, tp;
         while tps and not numberp car tps do tps := cdr tps;
         if tps and (tp := car tps) > 0 then limit := tp;
      end;

      realvalue := if 'rational memq types then
         (lambda(); {'quotient, random limit, den!-value limit})
      else
         (lambda(); random limit);
      value := if 'complex memq types then
         (lambda (); {'plus, apply(realvalue, nil),
            {'times, 'i, apply(realvalue, nil)}})
      else realvalue;
      print integervalue; print realvalue; print value;

      % Now assign some elements:
      for count := 1 : fix((m+n)/2) do <<
         i := random(m) + 1;
         j := random(n) + 1;
         puthash(i.j, hash, apply(value, nil));
      >>;
      return {'sparse!-mat, hash, m, n}
   end;

symbolic procedure den!-value limit;
   random(limit-1) + 1;

end;
