sparse_random_matrix(3, 5);

sparse_random_matrix(3, 5, 10);

sparse_random_matrix(3, 5, 1 .. 10);

sparse_random_matrix(3, 5, 10 .. 1);    % error?

sparse_random_matrix(3, 5, -10 .. -1);

sparse_random_matrix(3, 5, rational);

sparse_random_matrix(3, 5, complex);

sparse_random_matrix(3, 5, rational, complex);

sparse_random_matrix(3, 5, 10, density); % error?

sparse_random_matrix(3, 5, 10, density=50);

sparse_random_matrix(3, 5, 10, density = 0.5);

sparse_random_matrix(5, 5, 10, diagonal);

sparse_random_matrix(5, 5, 10, band);   % error?

sparse_random_matrix(5, 5, 10, band=3);

sparse_random_matrix(5, 5, 10, upper);

sparse_random_matrix(5, 5, 10, lower);

sparse_random_matrix(5, 5, 10, symmetric);

sparse_random_matrix(5, 5, 10, anti_symmetric);

sparse_random_matrix(5, 5, 10, hermitian);

sparse_random_matrix(5, 5, 10, anti_hermitian);

sparse_random_matrix(5, 5, 10, foo);   % error?

;end;
