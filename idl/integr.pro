function integr,f,rev=rev,x=x
;
;  Integral of f(t)
;   9-jun-92
;
n=n_elements(f)
a=make_array(size=size(f),/nozero)
;
;  check whether we want to multiply by dx (if so, then assume uniform mesh)
;
if n_elements(x) ne 0 then dx=x(1)-x(0) else dx=1.
;
if keyword_set(rev) then begin
  a(n-1)=0.
  for i=n-2,0L,-1L do a(i)=a(i+1)+.5*(f(i)+f(i+1))
endif else begin
  a(0)=0.
  for i=1L,n-1L do a(i)=a(i-1)+.5*(f(i)+f(i-1))
endelse
return,a*dx
end
