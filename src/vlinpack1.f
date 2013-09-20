C Output from Public domain Ratfor, version 1.01
      subroutine vqrdca(x,ldx,n,p,fasrkub3,jpvt,work,xwdf5ltg,eps)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      double precision dsign, dabs, dmax1, dsqrt
      integer min0
      integer ldx,n,p,xwdf5ltg
      integer jpvt(*)
      integer j,jj,jp,l,lup,curpvt
      double precision x(ldx,p),fasrkub3(p),work(*),eps
      double precision vdnrm2,tt
      double precision ddot8,nrmxl,t
      do23000 j=1,p 
      fasrkub3(j) = vdnrm2(n,x(1,j),ldx,1)
      work(j) = fasrkub3(j)
23000 continue
23001 continue
      l=1
      lup = min0(n,p)
      curpvt = p
23002 if(l.le.lup)then
      fasrkub3(l) = 0.0d0
      nrmxl = vdnrm2(n-l+1, x(l,l), ldx, 1)
      if(nrmxl .lt. eps)then
      call dshift8(x,ldx,n,l,curpvt)
      jp = jpvt(l)
      t=fasrkub3(l)
      tt=work(l)
      j=l+1
23006 if(.not.(j.le.curpvt))goto 23008
      jj=j-1
      jpvt(jj)=jpvt(j)
      fasrkub3(jj)=fasrkub3(j)
      work(jj)=work(j)
23007 j=j+1
      goto 23006
23008 continue
      jpvt(curpvt)=jp
      fasrkub3(curpvt)=t
      work(curpvt)=tt
      curpvt=curpvt-1
      if(lup.gt.curpvt)then
      lup=curpvt
      endif
      else
      if(l.eq.n)then
      goto 23003
      endif
      if(x(l,l).ne.0.0d0)then
      nrmxl = dsign(nrmxl,x(l,l))
      endif
      call dscal8(n-l+1,1.0d0/nrmxl,x(l,l),1)
      x(l,l) = 1.0d0+x(l,l)
      j=l+1
23015 if(.not.(j.le.curpvt))goto 23017
      t = -ddot8(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
      call daxpy8(n-l+1,t,x(l,l),1,x(l,j),1)
      if(fasrkub3(j).ne.0.0d0)then
      tt = 1.0d0-(dabs(x(l,j))/fasrkub3(j))**2
      tt = dmax1(tt,0.0d0)
      t = tt
      tt = 1.0d0+0.05d0*tt*(fasrkub3(j)/work(j))**2
      if(tt.ne.1.0d0)then
      fasrkub3(j) = fasrkub3(j)*dsqrt(t)
      else
      fasrkub3(j) = vdnrm2(n-l,x(l+1,j),ldx,1)
      work(j) = fasrkub3(j)
      endif
      endif
23016 j=j+1
      goto 23015
23017 continue
      fasrkub3(l) = x(l,l)
      x(l,l) = -nrmxl
      l=l+1
      endif
      goto 23002
      endif
23003 continue
      xwdf5ltg = lup
      return
      end
