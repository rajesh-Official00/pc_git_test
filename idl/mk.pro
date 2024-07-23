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
zr=integr(1./cs2,x=z)
zz=integr(gr/cs2,x=z)-depth
;=========

plot,z,zr, $
  linestyle=2, $
  xrange=[z[0], z[nz-1]], $
  yrange=[-3,3]
oplot,z,zz,ps=1
;

END
