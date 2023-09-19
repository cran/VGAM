









      subroutine daxpy8(n,da,dx,incx,dy,incy)
      implicit logical (a-z) 
c
c
c
c
c
c
c
c
      double precision dx(*),dy(*),da
      integer          i,incx,incy,m,mp1,n

      integer          ix, iy

c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c
c
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end



      subroutine  dcopy8(n,dx,incx,dy,incy)
      implicit logical (a-z) 
c
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c
c
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end



      double precision function ddot8(n,dx,incx,dy,incy)
c

      implicit logical (a-z) 

c
c
      double precision dx(*),dy(*),dtemp
      integer          i,incx,incy,ix,iy,m,mp1,n
c
      ddot8 = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot8 = dtemp
      return
c
c
c
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3)+dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot8 = dtemp
      return
      end



      double precision function dnrm28 ( n, dx,ldx, incx)
      implicit logical (a-z) 

      integer            n, ldx, incx, i, j, nn

      integer            next
      double precision   dx(ldx), cutlo, cuthi, hitest, sum,
     *                   xmax,zero,one

      data   zero, one /0.0d0, 1.0d0/
c
c
c
c
c
c
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm28  = zero
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
      dnrm28 = 0.0d0
      return
c
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
      dnrm28 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c
c
      dnrm28 = xmax * dsqrt(sum)
  300 continue
      return
      end


      subroutine  dscal8(n,da,dx,incx)
      implicit logical (a-z) 
c
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c
c
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end



      subroutine dshift8(x,ldx,n,j,k)
      implicit logical (a-z) 
      integer ldx,n,j,k
      double precision x(ldx,k), tt
      integer i,jj
      if(k.le.j)return
      do 100 i=1,n
      tt=x(i,j)
      do 50 jj=j+1,k
      x(i,jj-1)=x(i,jj)
 50   continue
      x(i,k)=tt
100   continue
      return
      end




      subroutine vdqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      implicit logical (a-z) 
      integer ldx,n,k,job,info
      double precision x(ldx,*),qraux(*),y(*),qy(*),qty(*),b(*),rsd(*),
     *                 xb(*)
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
      integer i,j,jj,ju,kp1
      double precision ddot8,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c
      info = 0
c
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0d0
      go to 250
   40 continue
c
c
         if (cqy) call dcopy8(n,y,1,qy,1)
         if (cqty) call dcopy8(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot8(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy8(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot8(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy8(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c
         if (cb) call dcopy8(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy8(k,qty,1,xb,1)
         if(cr .and. k .lt. n) call dcopy8(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy8(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot8(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy8(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -ddot8(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy8(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end


