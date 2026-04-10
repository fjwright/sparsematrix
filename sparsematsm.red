module sparsematsm;               % Simplification of sparse matrices.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-04-10 16:34:03 franc>
% Created: April 2026

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
%  * Redistributions of source code must retain the relevant copyright
%    notice, this list of conditions and the following disclaimer.
%
%  * Redistributions in binary form must reproduce the relevant
%    copyright notice, this list of conditions and the following
%    disclaimer in the documentation and/or other materials provided
%    with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% $Id$

% This file is a reworking of "matrix/matsm.red" to use hash tables to
% represent sparse matrices.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation and simplification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

put('sparse!-matrix, 'evfn, 'sparse!-matsm!*);

symbolic procedure sparse!-matsm!*(u,v);
   % Sparse matrix expression simplification function.
   sparse!-matsm!*1 sparse!-matsm u;

symbolic procedure sparse!-matsm!*1 u;
   % Assume u evaluates to a sparse matrix internal form
   %   (<hash> <m> <n>).
   % Simplify each element to an algebraic expression and return
   %   (sparse!-mat <hash> <m> <n>).
   % *** TEMPORARY HACK TO CHECK SIMPLER FACILITIES! ***
   'sparse!-mat . u;

symbolic procedure sparse!-matsm u;
   % If u evaluates to a variable assigned a sparse matrix
   % then return a sparse matrix internal form
   %   (<hash> <m> <n> . u).
   % If u evaluates to the form
   %   (sparse!-mat <hash> <m> <n>)
   % then return a sparse matrix internal form
   %   (<hash> <m> <n>).
   % The elements are SQs.  If u evaluates to an operator expression
   % (op s) then apply op to (<hash> <m> <n>) and return (<hash> <m>
   % <n>) where the elements are SQs.
   % *** TEMPORARY HACK TO CHECK SIMPLER FACILITIES! ***
   begin scalar x;
      if idp u and (x := get(u, 'avalue))
         and eqcar(x, 'sparse!-matrix)
            and eqcar(x := cadr x, 'sparse!-mat) then <<
               % Set name to u:
               rplacd(cdddr x, u);
               return cdr x
            >>
      else if eqcar(u, 'sparse!-mat) then
         return cdr u
      else return apply(car u, {sparse!-matsm(cadr u)});
   end;

% %%%%%%%%%
% Transpose
% %%%%%%%%%

symbolic procedure sparsetp u; sparse!-tp1 sparse!-matsm u;

flag('(sparsetp), 'sparse!-matflg);
put('sparsetp, 'rtypefn, 'getrtypecar); % declares algebraic operator

symbolic procedure sparse!-tp1 u;
   % Return the transpose of the sparse matrix internal form U =
   % (<hash> <m> <n>) as a new sparse matrix internal form.
   begin scalar
      alist := hashcontents(car u),
      % Each alist element has the form ((i j) value).
      newhash := mk!-sparse!-matrix!-hash();
      for each el in alist do           % write transposed element
         puthash({cadar el,caar el}, newhash, cadr el);
      return {newhash, caddr u, cadr u}
   end;

endmodule;

end;
