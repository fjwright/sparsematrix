m := mat((a,b),(c,d));
s := sparsify m;

sparse_matrix_p m;
sparse_matrix_p s;
sparse_square_matrix_p s;

sparse_symmetric_matrix_p s;
ss := s + sparse_tp s;
sparse_symmetric_matrix_p ss;
sas := s - sparse_tp s;
sparse_skew_symmetric_matrix_p sas;

m := mat((1+i,2),(3i,4+i));
s := sparsify m;
sh := s + sparse_tp << conj s >>;
sparse_hermitean_matrix_p sh;

sah := s - sparse_tp << conj s >>;
sparse_skew_hermitean_matrix_p sah;

;end;
