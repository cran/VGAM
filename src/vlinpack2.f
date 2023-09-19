









      double precision function vdnrm2 ( n, dx,ldx, incx)
c
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
c
      integer          next
      double precision dx(ldx), cutlo, cuthi, hitest, sum
      double precision xmax,zero,one
      data             zero, one /0.0d0, 1.0d0/
c
c
c
c
c
c
      data cutlo, cuthi / 8.232e-11,  1.304e19 /
c

      if(n .gt. 0) go to 10
         vdnrm2  = zero
         go to 300
c
   10 next = 30
      sum = zero
      nn = n * incx
      i = 1
   20 if(next .eq. 30) go to 30
      if(next .eq. 50) go to 50
      if(next .eq. 70) go to 70
      if(next .eq. 110) go to 110
      vdnrm2 = 0.0d0
      return

   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      next = 50 
      xmax = zero
c
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
      next = 70
      go to 105
c
c
  100 i = j
      next = 110
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c
   75 sum = (sum * xmax) * xmax
c
c
c
   85 hitest = cuthi / DBLE( n )
c
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
c
c
      vdnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end





      subroutine vdpbfa7(abd,lda,n,m,info,d)
      integer lda,n,m,info
      double precision abd(lda,*), d(n)
c
c
c
c
c
c
c
c
c
c
c
      double precision s,t
      integer ik,j,jk,k,mu, i,row
c
c
         d(1) = abd(m+1,1)
c
         do 30 j = 1, n
            info = j
            s = 0.0d0
            ik = m + 1
            jk = max0(j-m,1)
            mu = max0(m+2-j,1)
            if (m .lt. mu) go to 20
            do 10 k = mu, m
c
               t = abd(k,j)
               do 1 i = 1,k-mu
                  row = mu-2+i+j-m
                  t = t - d(row)*abd(ik-1+i,jk)*abd(mu-1+i,j)
   1           continue
c
               row = mu-2+(k-mu+1)+j-m
               t = t/d(row)
c
               abd(k,j) = t
c
               s = s + t*t*d(row)
c
               ik = ik - 1
               jk = jk + 1
   10       continue
   20       continue
            s = abd(m+1,j) - s
c
            if (s .le. 0.0d0) go to 40
c
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
c
c
c
c
c
c
c
      double precision ddot8,t
      integer k,kb,la,lb,lm
c
c
      do 10 k = 1, n
         lm = min0(k-1,m)
         la = m + 1 - lm
         lb = k - lm
         t = ddot8(lm,abd(la,k),1,b(lb),1)
c
         b(k) =  b(k) - t
c
   10 continue
c
c
      do 15 k = 1, n
         b(k) =  b(k)/d(k)
   15 continue
c
c
c
      do 20 kb = 1, n
         k = n + 1 - kb
         lm = min0(k-1,m)
         la = m + 1 - lm
         lb = k - lm
c
c
         t = -b(k)
         call daxpy8(lm,t,abd(la,k),1,b(lb),1)
   20 continue
      return
      end




