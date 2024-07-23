;pc_read_grid,obj=grid

nz=256 & cs2max=10. & depth=2*!pi ;(256strati_n1)
nz=256 & cs2max=10. & depth=!pi ;(n7)
nz=128 & cs2max=4. & depth=!pi ;(128nomag)
;gr=0.5
gr=1.
;
;w=.1 & z2=2.*!pi/6 & z1=z2-depth
w=.04 & z2=!pi/6 & z1=z2-depth
; z2 is the coordinate of upper boundary
w=.06 & z2=!pi/8 & z1=z2-depth
gam=5./3.
csbase=1.
z=grange(z1,z2,nz)
;
cs2=1.+.5*(1.+tanh(z/w))*(cs2max-1.)
;cs2=csbase+.5*(1.+tanh(z/w))*(cs2max-csbase)
;=========
;zz=integr(1./cs2,x=z)-depth
zz=integr(gr/cs2,x=z)-depth
;=========

;plot,z,z
;oplot,z,zz,ps=-1
;
lncs2=alog(cs2)
lnrho=-gam*zz-lncs2
ss=(gam-1)*zz+lncs2
;
!p.multi=[0,1,2]
plot_io,z,exp(lnrho),psym=-2
plot,z,ss,psym=-2
!p.multi=0
;
openw,1,'stratification.dat'
for n=0,nz-1 do begin
  printf,1,z(n),lnrho(n),ss(n)
endfor
close,1
END
