load_package sparse;

size := 500$
% Make an invertible random sparse matrix SM
% (with some symbolic entries):
sparse_random_matrix sm(size,size);
for i := 1:size do sm(i,i) := x^i;
% Make a dense matrix copy M:
m := densify sm$
% Make a sparse copy S:
s := m$ transmat s;

on time;

% Addition:
matrix msum(size,size);
for count := 1:100 do msum := msum + m; % MATRIX is slow (24969 ms)
sparse_matrix smsum(size,size);
for count := 1:100 do smsum := smsum + sm; % SPARSEMATRIX is OK (94 ms)
sparse ssum(<<size>>,<<size>>);      % doesn't accept variable size!!!
for count := 1:100 do ssum := ssum + s; % SPARSE is fastest (47 ms !!!)

% Multiplication:
m^10$                                  % MATRIX is slow (16687 ms)
sm^10$                                 % SPARSEMATRIX is fast (891 ms)
s^10$                                  % SPARSE is fastest (547 ms !!!)

% Determinant:
det m;         % MATRIX is OK (468 ms)
sparse_det sm; % SPARSEMATRIX is fastest (16 ms)
% det s;       % SPARSE crashes!
% FAILS even with size = 50: Heap exhausted (eventually)

% Inverse:
m^(-1)$                              % MATRIX is OK (1438 ms)
sm^(-1)$                             % SPARSEMATRIX is fast (62 ms)
% s^(-1)$                              % SPARSE crashes
% Only appears to work with numerical entries!

if m^(-1)*m = m*m^(-1) then "OK" else "***** ERROR *****"; % OK
if sm^(-1)*sm = sm*sm^(-1) then "OK" else "***** ERROR *****"; % OK
% if s^(-1)*s = s*s^(-1) then "OK" else "***** ERROR *****"; % FAILS
% Both s^(-1)*s and s*s^(-1) fail for different reasons!

;end;
