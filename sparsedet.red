module sparsedet;          % Determinant and trace of a sparse matrix.


% load_package matrix;              % uses some code in "matrix/det.red"

% %%%%%%%%%%%
% Determinant
% %%%%%%%%%%%

put('sparse_det, 'simpfn, 'simpsparse!-det);
flag('(sparse_det), 'immediate);

% Using expansion by minors.
% No support for Bareiss algorithm at present!

symbolic procedure simpsparse!-det u;
   % Return determinant of sparse matrix s, where u = (s) and s is a
   % variable assigned a sparse matrix, cf. det.
   sparse!-detq sparse!-matsm carx(u, 'sparse_det);

%% symbolic procedure sparse!-detq u;
%%    % Top level determinant function.
%%    % u is a sparse matrix canonical form (<hash> <m> <n>).
%%    begin scalar len := cadr u;          % Number of rows <m>.
%%       if caddr u neq len then rederr "Non square matrix";
%%       if len = 1 then
%%          return gethash({1,1}, car u) or 0;
%%       matrix_clrhash();
%%       u := sparse!-detq1(car u, len, 0, 1);
%%       matrix_clrhash();
%%       return u
%%    end;

%% symbolic procedure sparse!-detq1(hash, len, ignnum, i);
%%    % HASH contains elements of a sparse square matrix of (initial) order LEN.
%%    % The elements are assumed to be standard quotients.
%%    % Return the determinant of the matrix.
%%    % Algorithm is expansion by minors of first row.
%%    % IGNNUM is a packed set of column indices to avoid.
%%    % I is the current "first" row index, initially 1.
%%    % The row dimension remains LEN.
%%    begin scalar n2, sign, z;
%%       n2 := 1;
%%       if i = len then return
%%       begin scalar j := 1;
%%          while twomem(n2,ignnum)
%%             do << n2 := 2*n2; j := j+1 >>;
%%          return gethash({i,j}, hash)  % Last row, single element.
%%       end;
%%       if z := matrix_gethash ignnum then return cdr z;
%%       i := i + 1;
%%       z := nil ./ 1;                    % zero standard quotient
%%       for j := 1 : len do
%%          begin scalar x :=  gethash({i,j}, hash);
%%             if not twomem(n2,ignnum) then <<
%%                if x then <<
%%                   if sign then x := negsq x;
%%                   z := addsq(multsq(x, sparse!-detq1(hash, len, n2+ignnum, i)),
%%                      z)
%%                >>;
%%                sign := not sign
%%             >>;
%%             n2 := 2*n2
%%          end;
%%       % z is a standard quotient and hence never NIL, and that makes
%%       % this use of hash tables safe!
%%       matrix_puthash(ignnum,z);
%%       return z
%%    end;

symbolic procedure sparse!-detq u;
   % Top level determinant function.
   % U is a sparse matrix canonical form (<hash> <m> <n>).
   begin scalar len := cadr u;          % number of rows <m>
      if caddr u neq len then rederr "Non square sparse matrix";
      return if len = 1 then
         gethash({1,1}, car u) or 0
      else sparse!-detq1(car u, len, 1, for j := 1:len collect j)
   end;

symbolic procedure sparse!-detq1(hash, len, i, jlist);
   % HASH contains elements of a sparse square matrix of order LEN.
   % The elements are assumed to be standard quotients.
   % Return the determinant of the matrix.
   % Algorithm is recursive expansion by minors of first row.
   % I is the current "first" row index, initially 1, finally LEN.
   % JLIST is a list of column indices for the current sub-matrix,
   % initially {1,...,LEN}.
   begin scalar i1, result, neg;
      % Base case: last (or single) row, single element.
      if i = len then return gethash({i, car jlist}, hash);
      i1 := i + 1;
      result := nil ./ 1;               % zero standard quotient
      for each j in jlist do
      begin scalar el := gethash({i,j}, hash);
         if el then <<
            if neg then el := negsq el;
            result := addsq(result, multsq(el,
               % determinant of matrix excluding row i and column j:
               sparse!-detq1(hash, len, i1, delete(j,jlist))));
         >>;
         neg := not neg;
      end;
      return result
   end;

% %%%%%
% Trace
% %%%%%

put('sparse_trace, 'simpfn, 'simpsparse!-trace);

symbolic procedure simpsparse!-trace u;
   % Return trace of sparse matrix s, where u = (s) and s is a
   % variable assigned a sparse matrix, cf. trace.
   begin scalar m, hash, el, z;
      u := sparse!-matsm!*(carx(u,'trace),nil); % (sparse-mat <hash> <m> <n>)
      if (m := caddr u) neq cadddr u then rederr "Non square matrix";
      hash := cadr u;
      z := nil ./ 1;                    % zero standard quotient
      % Assume elements of sparse matrix are algebraic forms, so
      % simplify them to standard quotients.
      for i := 1 : m do
         if (el := gethash({i,i}, hash)) then z := addsq(simp el, z);
      return z
   end;

endmodule;

end;
