pc_power_xy,var1,var2,lint_shell=1,v1='uz_xy',obj=obj
spec1=reform(obj.spec1)
s=size(spec1)
nt=2L^fix(alog(s[2])/alog(2.))
;
;  reduce temporal data points to powers of 2
;
print,'original and reduced time dimension:',s[2],nt
spec1b=spec1(*,s[2]-nt:*)
fko=fft(spec1b,-1,dim=2)
nk=s[1]
k=findgen(nk)/4.
;
;  compute time array
;
tt=obj.tt(s[2]-nt:*)
dt=(max(tt)-min(tt))/(nt-1)
Ttot=dt*nt
dom=2.*!pi/Ttot
om=findgen(nt)*dom
;
;  coarsegrain
;
nevery=10
pc_coarsegrain,om,fko,nevery,oom,ffko;,/aver
;
;  clip high omega values
;
oom_max=max(k)
nom_max=max(where(oom le oom_max))
ffko=ffko(*,0:nom_max-1)
oom=oom(0:nom_max-1)
contour,alog(abs(ffko)),/fil,nlev=30,k,oom,yr=[0,1]*max(k)
END
