      subroutine dhkt9w(x,ldx,n,p,i0qvzl,jpvt,bgu6fw,cqui1v,eps)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      double precision dsign, dabs, dmax1, dsqrt
      integer min0
      integer ldx,n,p,cqui1v
      integer jpvt(1)
      integer j,jj,jp,l,lup,qxy4wd
      double precision x(ldx,p),i0qvzl(p),bgu6fw(1),eps
      double precision vdnrm2,tt
      double precision ddot8,nrmxl,t
      do 23000 j=1,p 
      i0qvzl(j) = vdnrm2(n,x(1,j),ldx,1)
      bgu6fw(j) = i0qvzl(j)
23000 continue
      l=1
      lup = min0(n,p)
      qxy4wd = p
23002 if(.not.(l.le.lup))goto 23003
      i0qvzl(l) = 0.0d0
      nrmxl = vdnrm2(n-l+1, x(l,l), ldx, 1)
      if(.not.(nrmxl .lt. eps))goto 23004
      call dshift8(x,ldx,n,l,qxy4wd)
      jp = jpvt(l)
      t=i0qvzl(l)
      tt=bgu6fw(l)
      j=l+1
23006 if(.not.(j.le.qxy4wd))goto 23008
      jj=j-1
      jpvt(jj)=jpvt(j)
      i0qvzl(jj)=i0qvzl(j)
      bgu6fw(jj)=bgu6fw(j)
       j=j+1
      goto 23006
23008 continue
      jpvt(qxy4wd)=jp
      i0qvzl(qxy4wd)=t
      bgu6fw(qxy4wd)=tt
      qxy4wd=qxy4wd-1
      if(.not.(lup.gt.qxy4wd))goto 23009
      lup=qxy4wd
23009 continue
      goto 23005
23004 continue
      if(.not.(l.eq.n))goto 23011
      goto 23003
23011 continue
      if(.not.(x(l,l).ne.0.0d0))goto 23013
      nrmxl = dsign(nrmxl,x(l,l))
23013 continue
      call dscal8(n-l+1,1.0d0/nrmxl,x(l,l),1)
      x(l,l) = 1.0d0+x(l,l)
      j=l+1
23015 if(.not.(j.le.qxy4wd))goto 23017
      t = -ddot8(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
      call daxpy8(n-l+1,t,x(l,l),1,x(l,j),1)
      if(.not.(i0qvzl(j).ne.0.0d0))goto 23018
      tt = 1.0d0-(dabs(x(l,j))/i0qvzl(j))**2
      tt = dmax1(tt,0.0d0)
      t = tt
      tt = 1.0d0+0.05d0*tt*(i0qvzl(j)/bgu6fw(j))**2
      if(.not.(tt.ne.1.0d0))goto 23020
      i0qvzl(j) = i0qvzl(j)*dsqrt(t)
      goto 23021
23020 continue
      i0qvzl(j) = vdnrm2(n-l,x(l+1,j),ldx,1)
      bgu6fw(j) = i0qvzl(j)
23021 continue
23018 continue
       j=j+1
      goto 23015
23017 continue
      i0qvzl(l) = x(l,l)
      x(l,l) = -nrmxl
      l=l+1
23005 continue
      goto 23002
23003 continue
      cqui1v = lup
      return
      end
