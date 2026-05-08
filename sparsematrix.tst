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

% Aggregate property:
abs s;
depend {a,b,c,d}, x;
df(s,x);

% Substitution:
sub(a=aa,s);

sparse_random_matrix(sr, 5, 5);

% Trace:
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
m3 := mat((a,b,c),(d,e,f));
m4 := mat((g,h),(i,j),(k,l));
s3 := sparsify m1;
s4 := sparsify m2;
m3*m4;
s3*s4;

m4*m3;
s4*s3;

2m1*m3*m4;
2s1*s3*s4;

m3*m4*m1*3;
s3*s4*s1*3;

% Positive integer powers:
m := mat((1,2),(3,4));
s := sparsify m;
m^10;
s^10;
% Zero and negative powers currently fail!

;end;
