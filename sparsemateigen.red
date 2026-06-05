module sparsemateigen; % Compute eigen-values & vectors of sparse matrix.

% Author: Eberhard Schruefer.
% Modification: James Davenport and Fran Burstall.
% Revised for sparse matrices represented as hash tables by FJW.
% Time-stamp: <2026-06-04 18:28:31 franc>

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

% This file is a reworking of part of "matrix/glmat.red" to use hash
% tables to represent sparse matrices.

fluid '(!*factor !*sqfree kord!*);

global '(!!arbint);

if null !!arbint then !!arbint := 0;

% algebraic operator arbcomplex;

% Done this way since it's also defined in the solve1 module.

deflist('((arbcomplex simpiden)),'simpfn);


flag('(sparse_mateigen),'opfn);

% flag('(sparse_mateigen),'noval);

symbolic procedure sparse_mateigen(u, eival);
   % U is an algebraic sparse matrix form,
   % EIVAL an indeterminate naming the eigenvalues.
   % Return a list of lists:
   %   {{eival-eq1,multiplicity1,eigenvector1},....},
   % where eival-eq is a polynomial and eigenvector is a matrix.
   begin scalar arbvars,exu,sgn,q,r,s,x,y,z,eivec,!*factor,!*sqfree,
         !*exp,!*rounded;
      integer l;
      !*exp := t;
      if not(getrtype u eq 'sparse!-matrix) then typerr(u,"sparse matrix");
      eival := !*a2k eival;
      kord!* := eival . kord!*;
      exu := sparse_mateigen1(sparse!-matsm u,eival);
      q := car exu;
      y := cadr exu;
      z := caddr exu;
      exu := cdddr exu;
      !*sqfree := t;
      for each j in cdr fctrf numr subs2(lc z ./ 1)
         do if null domainp car j and mvar car j eq eival
         then s := (if null red car j
         then !*k2f mvar car j . (ldeg car j*cdr j)
         else j) . s;
      for each j in q
         do (if x then rplacd(x,cdr x + cdr j)
         else s := (y . cdr j) . s)
            where x := assoc(y,s) where y := absf reorder car j;
      l := length s;
      % Build the return list:
      % S is a list (of length l) of pairs of the form
      %   <square-free factor> . <multiplicity>
      r := 'list .
         for each j in s collect <<  % square-free factor of char poly
            if null((cdr j = 1) and (l = 1)) then <<
               y := 1;
               for each k in exu do
                  if x := reduce!-mod!-eig(car j,c!:extmult(k,y))
                  then y := x
            >>;
            arbvars := nil;
            for each k in lpow z do
               if (y=1) or null(k member lpow y) then
                  arbvars := (k . makearbcomplex()) . arbvars;
            sgn := (y=1) or evenp length lpow y;
            eivec := 'mat . for each k in lpow z collect % REVISE
               {if x := assoc(k,arbvars) then mvar cdr x
               else prepsq!* mkgleig(k, y, sgn := not sgn, arbvars)};
            % {square-free factor, multiplicity, eigenvector}:
            {'list, prepsq!*(car j ./ 1), cdr j, eivec}
         >>;
      kord!* := cdr kord!*;
      return r
   end;

symbolic procedure sparse_mateigen1(u, eival);
   % U is a simplified sparse matrix canonical form, EIVAL an
   % indeterminate naming the eigenvalues.
   begin scalar q,x,y,z; integer l,lm,m;
      lm := length car u;               % column dimension
      z := 1;
      u := for each v in u collect <<   % for each row, v
         y := 1;            % find lcm of denominators of all elements, y
         for each w in v do y := lcm(y, denr w); % for each element, w
         m := lm;                                % col index
         l := l + 1;                    % row index
         x := nil;
         for each j in reverse v do <<  % for each element backwards, j
            if numr j or l = m then
               x := {m} .* multf(if l = m
               then
                  addf(numr j, negf multf(!*k2f eival, denr j))
               else numr j,
                     quotf(y, denr j)) .+ x;
            m := m - 1
         >>;
         y := z;
         z := c!:extmult(if null red x then <<
            q := (if p then (car p  . (cdr p + 1)) . delete(p,q)
            else (lc x  . 1) . q) where p = assoc(lc x,q);
            !*p2f lpow x>> else x, z);
         x
      >>;
      return q . y . z . u
   end;

endmodule;

end;
