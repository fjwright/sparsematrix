on errcont;
sparse_random_matrix(5);                % => square
sparse_random_matrix(5, 7);
sparse_random_matrix(5, 7, 10);         % error
sparse_random_matrix(5, 7, limit=10);
sparse_random_matrix(5, 7, 1 .. 10);
sparse_random_matrix(5, 7, 10 .. 1);    % error
sparse_random_matrix(5, 7, -10 .. -1);
sparse_random_matrix(5, 7, rational);
sparse_random_matrix(3, 5, complex);
sparse_random_matrix(5, 7, symbol);
sparse_random_matrix(5, 7, limit=10, symbol); % error
sparse_random_matrix(5, 7, rational, symbol); % error
sparse_random_matrix(5, 7, limit=10, rational, complex);
sparse_random_matrix(5, 7, limit=10, density); % error
sparse_random_matrix(5, 7, limit=10, density=100);
sparse_random_matrix(5, 7, limit=10, density=1.0);
sparse_random_matrix(5, 7, limit=10, density=0.5);
sparse_random_matrix(5, 7, limit=10, density=1/2);
sparse_random_matrix(5, 7, limit=10, diagonal); % error
sparse_random_matrix(5, limit=10, diagonal);
sparse_random_matrix(5, limit=10, diagonal, density=0.5);
sparse_random_matrix(5, limit=10, band); % => band=3
sparse_random_matrix(5, limit=10, band=3);
sparse_random_matrix(5, limit=10, upper);
sparse_random_matrix(5, limit=10, lower);
sparse_random_matrix(5, limit=10, invertible);
sparse_random_matrix(5, limit=10, symmetric);
sparse_random_matrix(5, symbol, symmetric);
sparse_random_matrix(5, limit=10, density=0.5, symmetric, invertible);
sparse_random_matrix(5, symbol, density=0.5, symmetric, invertible);
sparse_random_matrix(5, limit=10, density=0.5, anti_symmetric);
sparse_random_matrix(5, limit=10, density=0.5, anti_symmetric, invertible); % error
sparse_random_matrix(5, limit=10, density=0.5, hermitian);
sparse_random_matrix(5, limit=10, density=0.5, hermitian, invertible);
sparse_random_matrix(5, limit=10, density=0.5, anti_hermitian);
sparse_random_matrix(5, limit=10, density=0.5, anti_hermitian, invertible);
sparse_random_matrix(5, limit=10, density=0.5, foo); % error
;end;
