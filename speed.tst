in in$

load_package sparse;

% Make an invertible random sparse matrix RSM:
size := 42$                             % 40 OK, 45 too big
sparse_random_matrix rsm(size,size);
for i := 1:size do rsm(i,i) := 1;
% Make a dense matrix copy RM:
rm := densify rsm$
% Make a sparse copy RS:
rs := rm$ transmat rs;

on time;

det rm;
sparse_det rsm;                % OK with size = 100
det rs;                        % FAILS: Heap exhausted with size = 100

rm^10$
rsm^10$
rs^10$

rm^(-1)$
rsm^(-1)$
rs^(-1)$

;end;
