      subroutine vdgam1(x, ghry8z, c4uxow)
      implicit logical (a-z)
      double precision x, ghry8z
      integer c4uxow
      double precision w, sqap4b, temp
      if(.not.(x .le. 0.0d0))goto 23000
      c4uxow = 0
      return
23000 continue
      c4uxow = 1
23001 continue
      if(.not.(x .lt. 6.0d0))goto 23002
      call vdgam2(x + 6.0d0, temp, c4uxow)
      ghry8z = temp - 1.0d0/x - 1.0d0/(x + 1.0d0) - 1.0d0/(x + 2.0d0) - 
     &1.0d0/(x + 3.0d0) - 1.0d0/(x + 4.0d0) - 1.0d0/(x + 5.0d0)
      return
23002 continue
      w = 1.0d0 / (x * x)
      sqap4b = ((w * (-1.0d0/12.0d0 + ((w * (1.0d0/120.0d0 + ((w * (-1.
     &0d0/252.0d0 + ((w * (1.0d0/240.0d0 + ((w * (-1.0d0/132.0d0 + ((w *
     & (691.0d0/32760.0d0 + ((w * (-1.0d0/12.0d0 + (3617.0d0 * w)/8160.
     &0d0)))))))))))))))))))))
      ghry8z = ( dlog(x) - 0.5d0/x + sqap4b )
      return
      end
      subroutine vdgam2(x, ghry8z, c4uxow)
      implicit logical (a-z)
      double precision x, ghry8z
      integer c4uxow
      double precision w, sqap4b, temp
      if(.not.(x .le. 0.0d0))goto 23004
      c4uxow = 0
      return
23004 continue
      c4uxow = 1
23005 continue
      if(.not.(x .lt. 6.0d0))goto 23006
      call vdgam1(x + 6.0d0, temp, c4uxow)
      ghry8z = temp - 1.0d0/x - 1.0d0/(x + 1.0d0) - 1.0d0/(x + 2.0d0) - 
     &1.0d0/(x + 3.0d0) - 1.0d0/(x + 4.0d0) - 1.0d0/(x + 5.0d0)
      return
23006 continue
      w = 1.0d0 / (x * x)
      sqap4b = ((w * (-1.0d0/12.0d0 + ((w * (1.0d0/120.0d0 + ((w * (-1.
     &0d0/252.0d0 + ((w * (1.0d0/240.0d0 + ((w * (-1.0d0/132.0d0 + ((w *
     & (691.0d0/32760.0d0 + ((w * (-1.0d0/12.0d0 + (3617.0d0 * w)/8160.
     &0d0)))))))))))))))))))))
      ghry8z = ( dlog(x) - 0.5d0/x + sqap4b )
      return
      end
      subroutine vtgam1(x, ghry8z, c4uxow)
      implicit logical (a-z)
      double precision x, ghry8z
      integer c4uxow
      double precision w, sqap4b, temp
      if(.not.(x .le. 0.0d0))goto 23008
      c4uxow = 0
      return
23008 continue
      c4uxow = 1
23009 continue
      if(.not.(x .lt. 6.0d0))goto 23010
      call vtgam2(x + 6.0d0, temp, c4uxow)
      ghry8z = temp + 1.0d0/x**2 + 1.0d0/(x + 1.0d0)**2 + 1.0d0/(x + 2.
     &0d0)**2 + 1.0d0/(x + 3.0d0)**2 + 1.0d0/(x + 4.0d0)**2 + 1.0d0/(x +
     & 5.0d0)**2
      return
23010 continue
      w = 1.0d0 / (x * x)
      sqap4b = 1.0d0 + (w * (1.0d0/6.0d0 + (w * (-1.0d0/30.0d0 + (w * (
     &1.0d0/42.0d0 + (w * (-1.0d0/30.0d0 + (w * (5.0d0/66.0d0 + (w * (-
     &691.0d0/2370.0d0 + (w * (7.0d0/6.0d0 - (3617.0d0 * w)/510.0d0)))))
     &)))))))))
      ghry8z = 0.5d0 * w + sqap4b / x
      return
      end
      subroutine vtgam2(x, ghry8z, c4uxow)
      implicit logical (a-z)
      double precision x, ghry8z
      integer c4uxow
      double precision w, sqap4b, temp
      if(.not.(x .le. 0.0d0))goto 23012
      c4uxow = 0
      return
23012 continue
      c4uxow = 1
23013 continue
      if(.not.(x .lt. 6.0d0))goto 23014
      call vtgam1(x + 6.0d0, temp, c4uxow)
      ghry8z = temp + 1.0d0/x**2 + 1.0d0/(x + 1.0d0)**2 + 1.0d0/(x + 2.
     &0d0)**2 + 1.0d0/(x + 3.0d0)**2 + 1.0d0/(x + 4.0d0)**2 + 1.0d0/(x +
     & 5.0d0)**2
      return
23014 continue
      w = 1.0d0 / (x * x)
      sqap4b = 1.0d0 + (w * (1.0d0/6.0d0 + (w * (-1.0d0/30.0d0 + (w * (
     &1.0d0/42.0d0 + (w * (-1.0d0/30.0d0 + (w * (5.0d0/66.0d0 + (w * (-
     &691.0d0/2370.0d0 + (w * (7.0d0/6.0d0 - (3617.0d0 * w)/510.0d0)))))
     &)))))))))
      ghry8z = 0.5d0 * w + sqap4b / x
      return
      end
      subroutine dgam1w(x, ghry8z, n, c4uxow)
      implicit logical (a-z)
      integer n, c4uxow
      double precision x(n), ghry8z(n)
      integer i, lqhm2g
      c4uxow = 1
      do 23016 i=1,n 
      call vdgam1(x(i), ghry8z(i), lqhm2g)
      if(.not.(lqhm2g .ne. 1))goto 23018
      c4uxow = lqhm2g
23018 continue
23016 continue
      return
      end
      subroutine tgam1w(x, ghry8z, n, c4uxow)
      implicit logical (a-z)
      integer n, c4uxow
      double precision x(n), ghry8z(n)
      integer i, lqhm2g
      c4uxow = 1
      do 23020 i=1,n 
      call vtgam1(x(i), ghry8z(i), lqhm2g)
      if(.not.(lqhm2g .ne. 1))goto 23022
      c4uxow = lqhm2g
23022 continue
23020 continue
      return
      end
