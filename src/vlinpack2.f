c  This file contains modified code from Hastie and Tibshirani's
c  GAMFIT code, as well as a rational cholesky function or two.
c  All code here derives from linpack 
c  T.Yee 7/10/99




c This function was formerly real function dnrm2, but now converted
c to double precision
c Nb. changed "float(n)" to "dfloat(n)"

      double precision function vdnrm2 ( n, dx,ldx, incx)
c
c added by tyee 23/9/00:
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      integer          next
      double precision dx(ldx), cutlo, cuthi, hitest, sum
      double precision xmax,zero,one
      data             zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() double precision and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181e-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438e19
c     data cutlo, cuthi / 8.232e-11,  1.304e19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232e-11,  1.304e19 /
c

      if(n .gt. 0) go to 10
         vdnrm2  = zero
         go to 300
c
   10 next = 30
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
c  20    go to next,(30, 50, 70, 110)
   20 if(next .eq. 30) go to 30
      if(next .eq. 50) go to 50
      if(next .eq. 70) go to 70
      if(next .eq. 110) go to 110
c An error!!!
      vdnrm2 = 0.0d0
      return

   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      next = 50 
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      next = 70
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      next = 110
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
c 11/4/01: replacing "**2.0d0" by "**2" (three times in this file) 
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi / dfloat( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
        if(dabs(dx(j)) .ge. hitest) go to 100
        sum = sum + dx(j)**2
   95 continue
      vdnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      vdnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end



c ==============================================================
c This is modified linpack Fortran code 
c Changes marked with yyy 
c 23/9/99
c Works


      subroutine vdpbfa7(abd,lda,n,m,info,d)
      integer lda,n,m,info
      double precision abd(lda,*), d(n)
c
c
c
c 20130419: orig.:
c     double precision abd(lda,1), d(n)
c
c
c
c     vdpbfa7 is dpbfa8 but uses Rational Cholesky instead of ordinary 
c     Cholesky
c
c     abd = t(u) d u  where u is unit upper triangular and d is diagonal
c     the diagonal of d is stored where the 1's of the u would be stored
c
c     See dpbfa8 for more information 
c     d(1:n) is assigned the values of diag(d), and abd(m+1,) <- 1
c
c     Improvement yet to do: 
c         delete d and put its contents into abd(m+1,) (intrinsic 1's)
c
c     internal variables
c
c     double precision ddot8
      double precision s,t
      integer ik,j,jk,k,mu, i,row
c     begin block with ...exits to 40
c
c
c yyy
         d(1) = abd(m+1,1)
c
         do 30 j = 1, n
c           print *, "j = ", j
            info = j
            s = 0.0d0
            ik = m + 1
            jk = max0(j-m,1)
            mu = max0(m+2-j,1)
            if (m .lt. mu) go to 20
            do 10 k = mu, m
c              print *, "    k = ", k
c              t = abd(k,j) - ddot8(k-mu,abd(ik,jk),1,abd(mu,j),1)
c
               t = abd(k,j)
               do 1 i = 1,k-mu
                  row = mu-2+i+j-m
                  t = t - d(row)*abd(ik-1+i,jk)*abd(mu-1+i,j)
c                 print *, "    row = ", row
   1           continue
c
c yyy
c              t = t/abd(m+1,jk)
               row = mu-2+(k-mu+1)+j-m
c              print *, "    row = ", row
               t = t/d(row)
c
               abd(k,j) = t
c
c yyy
c              print *, "    index  = ", mu-1+i+j-m
               s = s + t*t*d(row)
c
               ik = ik - 1
               jk = jk + 1
   10       continue
   20       continue
            s = abd(m+1,j) - s
c
c     ......exit
            if (s .le. 0.0d0) go to 40
c
c yyy
c           abd(m+1,j) = dsqrt(s)
            abd(m+1,j) = 1d0
            d(j) = s
c
   30    continue
         info = 0
   40 continue
      return
      end



      subroutine vdpbsl7(abd,lda,n,m,b,d)
      integer lda,n,m
      double precision abd(lda,*),b(*),d(*)
c
c
c
c 20130419: orig:
c     double precision abd(lda,1),b(1),d(1)
c
c
c
c     vdpbsl7 is dpbsl8 but uses Rational Cholesky instead of ordinary 
c     Cholesky
c
c     See dpbsl8 for more information 
c
c     Improvement yet to do: 
c         delete d and put its contents into abd(m+1,) (intrinsic 1's)
c
c     internal variables
c
      double precision ddot8,t
      integer k,kb,la,lb,lm
c
c     solve trans(r)*y = b
c
      do 10 k = 1, n
         lm = min0(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         t = ddot8(lm,abd(la,k),1,b(lb),1)
c
c yyy
c        b(k) = (b(k) - t)/abd(m+1,k)
         b(k) =  b(k) - t
c
   10 continue
c
c
c yyy
      do 15 k = 1, n
         b(k) =  b(k)/d(k)
   15 continue
c
c
c     solve r*x = y
c
      do 20 kb = 1, n
         k = n + 1 - kb
         lm = min0(k-1,m)
         la = m + 1 - lm
         lb = k - lm
c
c yyy
c        b(k) = b(k)/abd(m+1,k)
c
         t = -b(k)
         call daxpy8(lm,t,abd(la,k),1,b(lb),1)
   20 continue
      return
      end




