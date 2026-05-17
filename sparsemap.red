module sparsemap;

% Revised version of the function map!-eval1 in "alg/map.red" taken
% from "eds/edspatch.red".

% Extend MAP/SELECT to other structures than list/matrix

symbolic procedure map!-eval1(o,q,fcn1,fcn2);
 % o       structure to be mapped.
 % q       map expression (univariate function).
 % fcn1    function for evaluating members of o.
 % fcn2    function computing results (e.g. aeval).
 map!-apply(map!-function(q,fcn1,fcn2),o);

symbolic procedure map!-function(q,fcn1,fcn2);
   begin scalar v,w;
   v := '!&!&x;
   if idp q
      and (get(q,'simpfn) or get(q,'number!-of!-args)=1)
   then <<w:=v; q:={q,v}>>
   else if eqcar(q,'replaceby) then
      <<w:=cadr q; q:=caddr q>>
   else
   <<w:=map!-frvarsof(q,nil);
      if null w then rederr "map/select: no free variable found" else
         if cdr w then rederr "map/select: free variable ambiguous";
      w := car w;
   >>;
   if eqcar(w,'!~) then w:=cadr w;
   q := sublis({w.v,{'!~,w}.v},q);
   return {'lambda,{'w},
      {'map!-eval2,'w,mkquote v,mkquote q,mkquote fcn1,mkquote fcn2}};
   end;

symbolic procedure map!-apply(f,o);
   if atom o then apply1(f,o)
   else (if m then apply2(m,f,o)
               else car o . for each w in cdr o collect apply1(f,w))
                   where m = get(car o,'mapfn);

symbolic procedure mapmat(f,o);
   'mat . for each row in cdr o collect
      for each w in row collect
         apply1(f,w);

put('mat,'mapfn,'mapmat);

% Additions by FJW:

symbolic procedure map!-sparse!-mat(f,o);
   'sparse!-mat . map!-sparse!-matrix(cdr o,
      (lambda w; apply1(f,w)));

put('sparse!-mat, 'mapfn, 'map!-sparse!-mat);

endmodule;

end;
