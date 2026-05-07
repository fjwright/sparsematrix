in "sparsematrix.red"$

load_package sparse;

% Make a sparse matrix RSM:
sparse_random_matrix(rsm,1000,1000);
% Make a dense matrix copy RM:
rm := densify rsm$
% Make a sparse copy RS (transmat is destructive!):
rs := densify rsm$ transmat rs;

on time;

trace rm;
trace rs;
sparse_trace rsm;

tp rm$
tp rs$
sparse_tp rsm$

det rm;
det rs;
sparse_det rsm;

;end;
