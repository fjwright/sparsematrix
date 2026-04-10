module sparsedet;          % Determinant and trace of a sparse matrix.

symbolic procedure simpsparse!-trace u;
   % Return trace of sparse matrix s, where u = (s) and s is a
   % variable assigned a sparse matrix, cf. trace.
   begin scalar m, hash, el, z;
      u := sparse!-matsm!*(carx(u,'trace),nil); % (sparse-mat <hash> <m> <n>)
      if (m := caddr u) neq cadddr u then rederr "Non square matrix";
      hash := cadr u;
      z := nil ./ 1;
      % Assume elements of sparse matrix s are prefix forms, so
      % simplify them to standard quotient.
      for i := 1 : m do
         if (el := gethash({i,i}, hash)) then z := addsq(simp el, z);
      return z
   end;

put('sparsetrace, 'simpfn, 'simpsparse!-trace);

endmodule;

end;
