% Input this file after "sparsematsm.red" or replace the line
put('sparse!-matrix, 'evfn, 'sparse!-matsm!*);
% with the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

put('matrix, 'evfn, 'generic!-matsm!*); % updates "matrix/matrix.red"
put('sparse!-matrix, 'evfn, 'generic!-matsm!*);

%% symbolic procedure generic!-matsm!*(u,v);
%%    % Generic matrix expression simplification function.
%%    % U is an arbitrary matrix expression in algebraic form.
%%    % Return a matrix expression in tagged algebraic form converted to
%%    % dense or sparse representation as appropriate.
%%    begin scalar type := getrtype u, result;
%%       if type eq 'matrix then
%%          (if not errorp
%%             (result := errorset!*({'matsm!*, mkquote u, mkquote v}, nil))
%%          then return car result)
%%       else if type eq 'sparse!-matrix then
%%          (if not errorp
%%             (result := errorset!*({'sparse!-matsm!*, mkquote u, mkquote v}, nil))
%%          then return car result)
%%       else typerr(u, "matrix");
%%       % Convert sparse matrices to dense and try again:
%%       u := sparse!-densify!-all u;
%%       return matsm!*(u,v);
%%    end;

symbolic procedure generic!-matsm!*(u,v);
   % Generic matrix expression simplification function.
   % U is an arbitrary matrix expression in algebraic form.
   % Return a matrix expression in tagged algebraic form converted to
   % dense or sparse representation as appropriate.
   begin scalar type := getrtype u, result, msg := !*msg;
      put('matrix, 'evfn, 'matsm!*);
      put('sparse!-matrix, 'evfn, 'sparse!-matsm!*);
      !*msg := nil;                     % suppress warning from nssimp
      if type eq 'matrix then
         (if not errorp
            (result := errorset2({'matsm!*, mkquote u, mkquote v}))
         then return car result)
      else if type eq 'sparse!-matrix then
         (if not errorp
            (result := errorset2({'sparse!-matsm!*, mkquote u, mkquote v}))
         then return car result)
      else typerr(u, "matrix");
      !*msg := msg;
      % Convert sparse matrices to dense and try again:
      u := sparse!-densify!-all u;
      result := matsm!*(u,v);
      put('matrix, 'evfn, 'generic!-matsm!*);
      put('sparse!-matrix, 'evfn, 'generic!-matsm!*);
      return result;
   end;

symbolic procedure sparse!-densify!-all u;
   % Recursively densify any sparse matrix variables.
   u and if atom u then
      if getrtype u eq 'sparse!-matrix then {'densify, u} else u
   else sparse!-densify!-all car u .
      sparse!-densify!-all cdr u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;end;

m := mat((1,2),(3,4));
s := sparsify m;

s+m;

***** Missing sparse-matrix in
mat((1,2),(3,4))$


[2  4]
[    ]
[6  8]


m+s;

% No result!
