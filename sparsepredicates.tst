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
sh := s + << conj sparse_tp s >>;
sparse_hermitian_matrix_p sh;

sah := s - << conj sparse_tp s >>;
sparse_skew_hermitian_matrix_p sah;

% Lipshutz
sparse_matrix so(3,3);
so(1,1) := so(3,2) := 1/9$
so(1,2) := so(3,1) := 8/9$
so(1,3) := so(2,2) := -4/9$
so(2,1) := so(3,3) := 4/9$
so(2,3) := -7/9$
so;
so * sparse_tp so;
sparse_identity_matrix_p ws;
sparse_orthogonal_matrix_p so;

sparse_matrix so(2,2);
so(1,1) := so(2,2) := cos theta$
so(2,1) := -(so(1,2) := sin theta)$
so;
for all th let (cos th)^2 + (sin th)^2 => 1;
so * sparse_tp so;
sparse_identity_matrix_p ws;
sparse_orthogonal_matrix_p so;

so(2,1) := -so(2,1)$
so(2,2) := -so(2,2)$
so;
sparse_orthogonal_matrix_p so;

sparse_matrix su(3,3);
su(1,1) := su(2,2) := 1$
su(1,2) := -(su(2,1) := i)$
su(1,3) := su(3,2) := -1 + i$
su(2,3) := su(3,1) := 1 + i$
su := su/2;
sparse_unitary_matrix_p su;

;end;

% Note: sparse_tp conj s fails; error in how conj is mapped over a
% sparse matrix.
