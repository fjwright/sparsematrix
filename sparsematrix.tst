on errcont;

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

% Implicit mapping:
abs s;
depend {a,b,c,d}, x;
df(s,x);

% Substitution:
sub(a=aa,s);

sparse_random_matrix sr(5,5);

% Trace:
sparse_trace sr;
trace sr;

densify sr;

% Transpose:
densify sparse_tp sr;
tp sr;
densify tp sr;

% Determinant:
m := mat((a,b,c),(d,e,f),(g,h,i));
s := sparsify m;
det m;
sparse_det s;
det s;
if det m = det s then "OK" else "***** ERROR *****";

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
s3 := sparsify m3;
s4 := sparsify m4;
m3*m4;
s3*s4;

m4*m3;
s4*s3;

2m1*m3*m4;
2s1*s3*s4;

m3*m4*m1*3;
s3*s4*s1*3;

% Positive integer powers:
m5 := mat((1,2),(3,4));
s5 := sparsify m5;
m5^10;
s5^10;

% Rank:
% REDUCE manual example
rank m3;
sparse_rank sparsify m3;
rank sparsify m3;

% Wikipedia example
m6 := mat((1,0,1),(0,1,1),(0,1,1));
rank m6;
sparse_rank sparsify m6;
rank sparsify m6;

% Wikipedia example
m6 := mat((1,1,0,2),(-1,-1,0,-2));
rank m6;
sparse_rank sparsify m6;
rank sparsify m6;

% Submatrices:
% m := mat((a,b,c),(d,e,f),(g,h,i));
% s := sparsify m;
m;
densify sparse_submatrix(s,2,2);        % not advertised!

% Cofactors:
sparse_cofactor(s,1,1);
cofactor(s,1,1);
% Determinant via cofactors of first row:
for j := 1 : 3 sum s(1,j)*cofactor(s,1,j);
if ws = det s then "OK" else "***** ERROR *****";

% Inverses and non-positive integer powers:
m5^0;
s5^0;
m5^(-1);
s5^(-1);
m5^(-1)*m5;
s5^(-1)*s5;
if ws = s5^0 then "OK" else "***** ERROR *****";
m5*m5^(-1);
s5*s5^(-1);
if ws = s5^0 then "OK" else "***** ERROR *****";
m5^(-2);
s5^(-2);

% 1*1 zero matrix:
matrix m0(1,1);
m0^0;
sparsify m0;
sparse_matrix s0(1,1);
s0;
densify s0;
s0^0;
m0^(-1);
s0^(-1);
m5*m0;
s5*s0;
m5/m0;
s5/s0;

% 1*1 nonzero matrix:
m0(1,1) := 42;
s0(1,1) := 42;
m5*m0;
s5*s0;
m5/m0;
s5/s0;

% Explicit mapping:
m := mat((x^2,x^5),(x^4,x^5));
s := sparsify m;
map(int(~w,x), m);
map(int(~w,x), s);
map(1 + ~w, m);
map(1 + ~w, s);

% Nullspace:
m := mat((1,2,3,4),(5,6,7,8));
s := sparsify m;
nullspace mat2list m;
sparse_nullspace s;
nullspace s;

sparse_random_matrix s(5,10);
m := densify s;
nullspace mat2list m;
if second length m - length ws = rank m then "OK" else "***** ERROR *****";
nullspace s;
if second length s - length ws = rank s then "OK" else "***** ERROR *****";

% Mixed matrix types in algebraic expressions:

off sparse_matrix_dense_print;          % to show matrix type

m := mat((1,2),(3,4));
s := sparsify m;

sparse_matrix_auto_convert_mode();
(2s+m)*m;
(2m+s)^2;

sparse_matrix_auto_convert_mode sparse;
(2s+m)*m;
(2m+s)^2;

sparse_matrix_auto_convert_mode none;
(2s+m)*m;
(2m+s)^2;

;end;
