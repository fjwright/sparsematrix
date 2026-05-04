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

sparse_random_matrix(sr, 5, 5);

%Trace:
sparse_trace sr;

densify sr;

% Transpose:
densify sparse_tp sr;

% Determinant:
m := mat((a,b,c),(d,e,f),(g,h,i));
det m;
s := sparsify m;
sparse_det s;

% Addition:
m1 := mat((a,b),(c,d));
m2 := mat((e,f),(g,h));
m1 + m2;
s1 := sparsify m1;
s2 := sparsify m2;
% symbolic;
% global '(ss1 ss2 ss);
% ss1 := sparse!-matsm 's1;
% ss2 := sparse!-matsm 's2;
% ss := sparse!-addm(ss1,ss2);
% sparse!-matpri sparse!-matsm!*1 ss;
s1 + s2;


;end;
