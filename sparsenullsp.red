module sparsenullsp;  % Compute the nullspace (basis vectors) of a sparse matrix.

% Author: Herbert Melenk <melenk@sc.zib-berlin.de>.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the relevant copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%      documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNERS OR
% CONTRIBUTORS
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

% $Id: nullsp.red 5874 2021-07-30 21:08:56Z arthurcnorman $

% Algorithm: Rational Gaussian elimination with standard quotients.

put('sparse_nullspace, 'psopfn, 'sparse!-nullspace!-eval);

symbolic procedure sparse!-nullspace!-eval u;
   % Interface for the nullspace calculation.
   begin scalar v, matinput;
      v := reval car u;
      if eqcar(v, 'sparse!-mat) then
         << matinput:=t; v := cdr v >>
      else typerr ("sparse!-matrix", u);
      v := sparse!-nullspace!-alg v;
      return 'list . for each vect in v collect
         if matinput then 'mat . for each x in vect collect list x
         else 'list . vect;
   end;

symbolic procedure sparse!-nullspace!-alg u;
   % U is a sparse matrix, encoded as (<hash> <m> <n>) where the
   % matrix elements are algebraic expressions.
   % Result is the basis of the kernel of U in the same encoding.
   begin scalar mp, vars, rvars, r, res, oldorder;
      scalar hash := car u;
      integer n := caddr u;
      vars := for i := 1 : n collect gensym();
      rvars := reverse vars;
      oldorder := setkorder rvars;
      % Build a list of SQs, one for each row, where each SQ is a sum
      % of terms of the form V*<matrix element> and V is a new
      % variable (gensym) that is distinct for each column.
      mp := for i := 1 : cadr u collect <<
         r := nil ./ 1;                 % SQ zero
         for j := 1 : n do
            begin scalar v := vars,
                  el := gethash(i.j, hash);
               if el then r := addsq(r, simp {'times, car v, el});
               v := cdr v;
            end;
         r
      >>;
      res := nullspace!-elim(mp, rvars);
      setkorder oldorder;
      return reverse for each q in res collect
         for each x in vars collect
            cdr atsoc(x, q);
   end;

% symbolic procedure nullspace!-elim(m,vars) is defined in
% "matrix/nullsp.red".

endmodule;

end;
