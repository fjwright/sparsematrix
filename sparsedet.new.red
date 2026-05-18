module newsparsedet;                 % Determinant of a sparse matrix.

% Author: Francis J. Wright <https://sourceforge.net/u/fjwright>
% Time-stamp: <2026-05-18 15:02:55 franc>
% Created: April 2026

put('sparse_det, 'simpfn, 'simpsparse!-det);
flag('(sparse_det), 'immediate);

% Using expansion by cofactors.
% No support for Bareiss algorithm at present!

symbolic procedure simpsparse!-det u;
   % Return the determinant of a sparse matrix, cf. det.
   sparse!-detq sparse!-matsm carx(u, 'sparse_det);

% Access nonzero sparse cofactor elements via a jump-list containing
% elements of the form (i j j_hash), where i is a row index within the
% current cofactor, j is a column index within the current cofactor,
% and j_hash is the corresponding column index in the hash table.  All
% cofactor element values are accessed via the hash table.

% Could extract only the hash keys instead of calling hashcontents!

symbolic procedure sparse!-detq u;
   % Top level determinant function.
   % U is a sparse matrix canonical form (<hash> <m> <n>).
   begin scalar alist, len := cadr u;   % number of rows <m>
      if caddr u neq len then rederr "Non square sparse matrix";
      if len = 1 then return gethash(1 . 1, car u) or 0;
      alist := hashcontents car u;
      % ...a list of pairs of the form ((i . j) . value)
      alist := for each el in alist collect {caar el, cdar el, cdar el};
      % ...initially, j = j_hash; access values via hash
      return sparse!-detq1(car u, len, 1, alist);
   end;

symbolic procedure sparse!-detq1(hash, len, i, jlist);
   % HASH contains elements of a sparse square matrix of order LEN.
   % The elements are assumed to be standard quotients.
   % Return the determinant of the matrix.
   % Algorithm is recursive expansion by cofactors of first row.
   % I is the current "first" row index, initially 1, finally LEN.
   % JLIST is a list of indices for nonzero elements of the current
   % sub-matrix of the form ((i j j_hash)...), initially ((i j j)...)
   begin scalar i1, result := nil ./ 1; % zero SQ
      % Base case: last row, single element.
      if i = len then return
         if jlist then gethash(i . caddr car jlist, hash)
         else  result;
      % General case: sum of each element in row i times its cofactor.
      i1 := i + 1;                      % first row of sub-matrix
      for each ind in jlist do          % jlist is not sorted!
         if car ind = i then
         begin scalar
            j := cadr ind,
            j_hash := caddr ind,
            el := gethash(i . j_hash, hash),
            jlist_cofac := sparse!-cofactor!-jumplist(jlist, i, j);
            if evenp j then el := negsq el;
            result := addsq(result, multsq(el,
               % determinant of matrix excluding row i and column j:
               sparse!-detq1(hash, len, i1, jlist_cofac)));
         end;
      return result
   end;

symbolic procedure sparse!-cofactor!-jumplist(jlist, i, j);
   % Construct the new jump-list for the cofactor of element (I,J),
   % i.e. omitting row I (and implicitly higher rows) and column J,
   % and renumbering to the right of column J.
   % JLIST is a list of indices for nonzero elements of the current
   % sub-matrix of the form ((ii jj j_hash)...), where ii and j_hash
   % cannot be changed, but an element can be omitted.
   for each x in jlist join if car x > i then
   begin scalar jj := cadr x;
      return
         if jj < j then {x}
         else if jj > j then {{car x, jj-1, caddr x}};
   end;

endmodule;

end;
