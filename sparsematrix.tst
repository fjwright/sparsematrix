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

% Addition and scalar multiplication:
m1 := mat((a,b),(c,d));
m2 := mat((e,f),(g,h));
s1 := sparsify m1;
s2 := sparsify m2;
% symbolic;
% global '(ss1 ss2 ss);
% ss1 := sparse!-matsm 's1;
% ss2 := sparse!-matsm 's2;
% ss := sparse!-addm(ss1,ss2);
% sparse!-matpri sparse!-matsm!*1 ss;
m1 + m2;
s1 + s2;
2m1 + 3m2;
2s1 + 3s2;

% Matrix multiplication:
m1 := mat((a,b,c),(d,e,f));
m2 := mat((g,h),(i,j),(k,l));
s1 := sparsify m1;
s2 := sparsify m2;
symbolic;
global '(ss1 ss2 ss);
ss1 := sparse!-matsm 's1;
ss2 := sparse!-matsm 's2;
ss := sparse!-multm(ss1,ss2);
sparse!-matpri sparse!-matsm!*1 ss;
algebraic;
m1*m2;

;end;
