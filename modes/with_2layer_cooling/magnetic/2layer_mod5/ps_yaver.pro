if !d.name eq 'PS' then begin
  device,filename='ko.eps',xsize=18,ysize=14,yoffset=3,/color,/encapsul
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end

!p.charsize=1.7
!x.margin=[6.5,1.0]
!y.margin=[3.2,.5]

@parameters
irev=1
if irev eq 1 then begin
 loadct,5
 tvlct,r,g,b,/get
 rr=reverse(r)
 gg=reverse(g)
 bb=reverse(b)
 tvlct,rr,gg,bb
endif else begin
 loadct,5
endelse
;
default,iread,1
if iread eq 1 then begin
 pc_read_yaver,obj=yav
 pc_read_param,obj=param,/param2
 s=size(yav.uzmxz)
 nx=s[1] & nz=s[2] & nt=s[3]
 tt=yav.t
 Lxyz=param.Lxyz
 gr=-param.gravz
 cp=param.cp
 lx=Lxyz[0]
 ly=Lxyz[1]
 lz=Lxyz[2]
endif
iread=0
;

default,iz1,225
default,t1,40
default,t2,max(tt)

gd=where(tt gt t1 and tt lt t2)
tgd=tt[gd]
t_tot=max(tgd)-min(tgd)

dkx=2.*!pi/lx
dom=2.*!pi/t_tot

ngd=n_elements(tgd)
is_ngd_even=(ngd mod 2) eq 0
if (is_ngd_even) then $
  om=[indgen(ngd/2+1),-reverse(indgen(ngd/2-1)+1)]*dom $
else $
  om=[indgen(ngd/2+1),-reverse(indgen(ngd/2)+1)]*dom

kx=[indgen(nx/2+1),-reverse(indgen(nx/2-1)+1)]*dkx

; data, its Fourier transform and spectrum
d1=reform(yav.uzmxz[*,iz1,gd])
f1=fft(d1,-1)
s1=alog(abs(f1)^2)
print,'minmax(s1)=',minmax(s1)
default,l1,-25
default,l2,-12
lev=grange(l1,l2,21)
; cc=clip(s1,minmax(lev))

tilde='!9!s!aA!n!r!6'
!x.title=tilde+'!8k!dx!n!6'
!y.title=tilde+'!7x!6'
xr=[0,10]
yr=[0,10]
yr=[0,8]

; contour,cc,kx,om,/fil,lev=lev,xr=xr,yr=yr,background=0,col=255
contour,s1,kx,om,/fil,lev=lev,xr=xr,yr=yr,background=0,col=255
om_f=sqrt(gr*kx)
oplot,kx,om_f,li=4,col=255,thick=5

;cwd,run
;print,'mv ko.eps kyo_'+str+run+'.eps'
; stop
end
