module sparsedet;          % Determinant and trace of a sparse matrix.

symbolic procedure simpsparsetrace u;
   % Return trace of sparse matrix s, where u = (s) and s is a
   % variable assigned a sparse matrix, cf. trace.
   begin scalar m, hash, el, z;
      u := sparsematsm!*(carx(u,'trace),nil); % (sparsemat <hash> <m> <n>)
      if (m := caddr u) neq cadddr u then rederr "Non square matrix";
      hash := cadr u;
      z := nil ./ 1;
      % Assume elements of sparse matrix s are prefix forms, so
      % simplify them to standard quotient.
      for i := 1 : m do
         if (el := gethash({i,i}, hash)) then z := addsq(simp el, z);
      return z
   end;

put('sparsetrace, 'simpfn, 'simpsparsetrace);

endmodule;

end;
