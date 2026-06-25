% A sparse analogue of the LINALG random_matrix operator but with more
% flexibility without using switches.

% sparse_random_matrix id(m,n,options) where options can be one or
% more of:

% Number type:

%   number meaning limit
%   {integer, rational, real} (where real assumes on rounded)
%   complex

% Matrix type:

%   diagonal, band(number), upper, lower,
%   symmetric, anti_symmetric/skew_symmetric,
%   hermitian, anti_hermitian/skew_hermitian
