% Input this file after "sparsematsm.red" or replace the line
% put('sparse!-matrix, 'evfn, 'sparse!-matsm!*);
% with the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

put('matrix, 'evfn, 'generic!-matsm!*); % updates "matrix/matrix.red"
put('sparse!-matrix, 'evfn, 'generic!-matsm!*);

global '(sparse_matrix_auto_convert);
share sparse_matrix_auto_convert;
sparse_matrix_auto_convert := 'matrix; % or sparse_matrix or anything else

symbolic procedure generic!-matsm!*(u, v);
   % Generic matrix expression simplification function.
   % U is an arbitrary matrix expression in algebraic form.
   % Return a matrix expression in tagged algebraic form converted to
   % dense or sparse representation as appropriate.
   begin scalar type := getrtype u, result;
      put('matrix, 'evfn, 'matsm!*);
      put('sparse!-matrix, 'evfn, 'sparse!-matsm!*);
      sparse_matrix_auto_convert := reval sparse_matrix_auto_convert;
                                        % in case set in algebraic mode!
      result := if not
         (sparse_matrix_auto_convert memq '(matrix sparse_matrix)) or
         sparse!-check!-rtype(u, type) then
            if type eq 'matrix then
               matsm!*(u, v)
            else if type eq 'sparse!-matrix then
               sparse!-matsm!*(u, v)
            else typerr(u, "matrix")
      else
         if sparse_matrix_auto_convert eq 'matrix then <<
            % Convert all sparse matrices to dense:
            u := sparse!-densify!-all u;
            matsm!*(u,v)
         >> else if sparse_matrix_auto_convert eq 'sparse_matrix then <<
            % Convert all dense matrices to sparse:
            u := sparse!-sparsify!-all u;
            sparse!-matsm!*(u,v)
         >>;
      put('matrix, 'evfn, 'generic!-matsm!*);
      put('sparse!-matrix, 'evfn, 'generic!-matsm!*);
      return result;
   end;

symbolic procedure sparse!-check!-rtype(u, type);
   % Return t if the type of every matrix in prefix form algebraic
   % expression U is TYPE, nil otherwise.
   if null u then t
   else if atom u then
      (null x or x eq type) where x = getrtype u
   else sparse!-check!-rtype(car u, type) and
      sparse!-check!-rtype(cdr u, type);

symbolic procedure sparse!-densify!-all u;
   % Recursively densify any sparse matrix variables.
   u and if atom u then
      if getrtype u eq 'sparse!-matrix then {'densify, u} else u
   else sparse!-densify!-all car u .
      sparse!-densify!-all cdr u;

symbolic procedure sparse!-sparsify!-all u;
   % Recursively sparsify any demse matrix variables.
   u and if atom u then
      if getrtype u eq 'matrix then {'sparsify, u} else u
   else sparse!-sparsify!-all car u .
      sparse!-sparsify!-all cdr u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;end;

m := mat((1,2),(3,4));
s := sparsify m;

(2s+m)*m;
(2m+s)^2;
