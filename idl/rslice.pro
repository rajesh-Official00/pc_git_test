;
;; pc_build -t read_videofiles
;; ./src/read_videofiles.x
;; uu3


pc_read_param,obj=param
pc_read_video,field='uu3',obj=uz,njump=4000,nt=4000
; top slice; xy2 plane

s=size(uz.xy2)
nx=s[1] & ny=s[2] & nt=s[3]

;f1=fft(uz.xy2,-1)

s1=alog10(abs(f1))

gr=-param.gravz
Lxyz=param.Lxyz
lx=Lxyz[0]
ly=Lxyz[1]
lz=Lxyz[2]

t_tot=max(uz.t)-min(uz.t)
dkx=2.*!pi/lx
dom=2.*!pi/t_tot

ngd=nt

is_ngd_even=(ngd mod 2) eq 0
if (is_ngd_even) then $
  om=[indgen(ngd/2+1),-reverse(indgen(ngd/2-1)+1)]*dom $
else $
  om=[indgen(ngd/2+1),-reverse(indgen(ngd/2)+1)]*dom

kx=[indgen(nx/2+1),-reverse(indgen(nx/2-1)+1)]*dkx

ky=kx

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

print,'minmax(s1)=',minmax(s1)
l1=-12
l2=-3.5
lev=grange(l1,l2,21)
iky=0
cc=clip(reform(s1[*,iky,*]),minmax(lev))

tilde='!9!s!aA!n!r!6'
!x.title=tilde+'!8k!dx!n!6'
!y.title=tilde+'!7x!6'
xr=[0,10]
yr=[0,10]
yr=[0,7.5]

contour,cc,kx,om,/fil,lev=lev,xr=xr,yr=yr,background=0,col=255
om_f=sqrt(gr*kx)
oplot,kx,om_f,li=4,col=255,thick=5

end
