in "sparselinalg.red"$

on errcont;

m := mat((1,2),(3,4),(5,6));
s := sparsify m;

% %%%%%%%%%%%%%%%%%%%%%
% sparse_matrix_augment
% %%%%%%%%%%%%%%%%%%%%%
% cf. LINALG matrix_augment

sparse_matrix_augment();
sparse_matrix_augment(42);
s1 := sparse_matrix_augment(m, 2s);
sparse_matrix_augment({m, 2s});
sparse_matrix_augment({m, 2s}, {2m, 3s});

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sparse_select_columns / sparse_augment_columns
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cf. LINALG augment_columns

sparse_select_columns();
sparse_select_columns(42);
sparse_select_columns(s1, 2);
sparse_select_columns(s1, 0);
sparse_select_columns(s1, 13);
sparse_select_columns(s1, 2, {-2});
sparse_select_columns(s1, {2, 3});
sparse_select_columns(s1, 2 .. -1);
sparse_select_columns(s1, 1 .. 3, 3 .. 1);
sparse_select_columns(s1, 3 .. 4);

% %%%%%%%%%%%%%%%%%%%%%
% sparse_remove_columns
% %%%%%%%%%%%%%%%%%%%%%
% cf. LINALG remove_columns

sparse_remove_columns();
sparse_remove_columns(42);
sparse_remove_columns(s1, 2);
sparse_remove_columns(s1, 0);
sparse_remove_columns(s1, 13);
sparse_remove_columns(s1, 2, {-2});
sparse_remove_columns(s1, {2, 3});
sparse_remove_columns(s1, 2 .. -1);
sparse_remove_columns(s1, 1 .. 3, 3 .. 1);
sparse_remove_columns(s1, 3 .. 4);

;end;
