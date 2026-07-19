sparse_random_matrix(3);
sparse_random_matrix(3, 5);
sparse_random_matrix(3, 5, 10);         % error
sparse_random_matrix(3, 5, limit 10);   % ERROR NEEDS FIXING!
sparse_random_matrix(3, 5, 1 .. 10);
sparse_random_matrix(3, 5, 10 .. 1);    % error
sparse_random_matrix(3, 5, -10 .. -1);
sparse_random_matrix(3, 5, rational);
sparse_random_matrix(3, 5, complex);
sparse_random_matrix(3, 5, rational, complex);
sparse_random_matrix(3, 5, -10 .. 10, density);     % error
sparse_random_matrix(3, 5, -10 .. 10, density 100);   % WRONG?
sparse_random_matrix(3, 5, -10 .. 10, density 1.0);
sparse_random_matrix(3, 5, -10 .. 10, density 0.5);
sparse_random_matrix(3, 5, -10 .. 10, density(1/2));
sparse_random_matrix(5, 5, -10 .. 10, diagonal);               % ???
sparse_random_matrix(5, 5, -10 .. 10, diagonal, density 1.00); % ???
sparse_random_matrix(5, 5, -10 .. 10, band);                   % error
sparse_random_matrix(5, 5, -10 .. 10, band 3);
sparse_random_matrix(5, 5, -10 .. 10, upper);
sparse_random_matrix(5, 5, -10 .. 10, lower);
sparse_random_matrix(5, 5, -10 .. 10, symmetric);
sparse_random_matrix(5, 5, -10 .. 10, anti_symmetric);
sparse_random_matrix(5, 5, -10 .. 10, hermitian);
sparse_random_matrix(5, 5, -10 .. 10, anti_hermitian);
sparse_random_matrix(5, 5, -10 .. 10, foo); % error
;end;
