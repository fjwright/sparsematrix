in "sparsematrix.red"$


sparse_matrix sd(5,5);

length sd;

sd;

for i := 1:5 do sd(i,i) := i;

sd;

for i := 1:5 do write sd(i,i) := sd(i,i);

densify sd;

m := mat((a,b),(c,d));

sparsify m;

s := ws;

random_sparse_matrix(sr, 10, 10);

sr;

%Trace:
sparse_trace sr;

% Transpose:
sparse_tp sr;                           % FAILS!


;end;
