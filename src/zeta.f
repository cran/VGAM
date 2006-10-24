      subroutine vzetawr(x, ans, deriv, nn)
      double precision x(nn), ans(nn), zeta8, dzeta8, ddzeta8
      integer deriv, nn, i
      double precision b2(12)
      call vbecoef(b2)
      if(.not.(deriv.eq.0))goto 23000
      do 23002 i=1,nn 
      ans(i) = zeta8(x(i), b2)
23002 continue
23000 continue
      if(.not.(deriv.eq.1))goto 23004
      do 23006 i=1,nn 
      ans(i) = dzeta8(x(i), b2)
23006 continue
23004 continue
      if(.not.(deriv.eq.2))goto 23008
      do 23010 i=1,nn 
      ans(i) = ddzeta8(x(i), b2)
23010 continue
23008 continue
      return
      end
      double precision function zeta8(s, b2)
      double precision s
      double precision b2(12)
      integer a, k
      integer m, n, m2
      double precision sum, p, a2, fred
      a = 12
      k = 8
      a2 = a * a
      p = s / 2.0d0 / a2
      sum = 1.0d0 / (s - 1.0d0) + 0.5d0 / a + b2(1) * p
      do 23012 m = 2,k 
      m2 = m + m
      p = p * (s + m2 - 3.0d0) * (s + m2 - 2.0d0) / (m2 - 1.0d0) / m2 / 
     &a2
      sum = sum + p * b2(m)
23012 continue
      fred = dexp((s - 1.0d0) * dlog( dfloat(a) ))
      sum = 1.0d0 + sum / fred
      do 23014 n = 2,(a - 1) 
      sum = sum + dexp(-s * dlog( dfloat(n) ))
23014 continue
      zeta8 = sum
      return
      end
      double precision function dzeta8(s, b2)
      double precision s
      double precision b2(12)
      integer a, k
      integer m, n, m2
      double precision sum, p, q, a2, loga, logn
      double precision fred
      a = 12
      k = 8
      loga = dlog( dfloat(a) )
      a2 = a * a
      p = s / 2.0d0 / a2
      q = 1.0d0 / s - loga
      sum = b2(1) * p * q
      do 23016 m = 2,k 
      m2 = m + m
      p = p * (s + m2 - 3.0d0) * (s + m2 - 2.0d0) / (m2 - 1.0d0) / m2 / 
     &a2
      q = q + 1.0d0 / (s + m2 - 3.0d0) + 1.0d0 / (s + m2 - 2.0d0)
      sum = sum + b2(m) * p * q
23016 continue
      fred = dexp((1.0d0 - s) * loga)
      sum = (sum - 1.0d0/ (s - 1.0d0)**2 - loga * (1.0d0/(s - 1.0d0) + 
     &0.50d0/a)) * fred
      do 23018 n = 2, (a - 1) 
      logn = dlog( dfloat(n) )
      sum = sum - logn / dexp(logn * s)
23018 continue
      dzeta8 = sum
      return
      end
      double precision function upsilon8(s, b2)
      double precision s, dzeta8, zeta8, b2(12)
      upsilon8 = -dzeta8(s, b2) / zeta8(s, b2)
      return
      end
      double precision function ddzeta8(s, b2)
      double precision s, b2(12)
      integer a, k
      integer m, n, m2
      double precision sum, p, q, r, a2, loga, logn
      double precision fred, fred2
      a = 12
      k = 8
      loga = dlog( dfloat(a) )
      a2 = a * a
      p = s / 2.0d0 / a2
      q = 1.0d0 / s - loga
      r = 1.0d0 / s / s
      sum = b2(1) * p * (q * q - r)
      do 23020 m = 2,k 
      m2 = m + m
      p = p * (s + m2 - 3.0d0) * (s + m2 - 2.0d0) / (m2 - 1.0d0) / m2 / 
     &a2
      q = q + 1.0d0 / (s + m2 - 3.0d0) + 1.0d0 / (s + m2 - 2.0d0)
      r = r + 1.0d0 / (s + m2 - 3.0d0)**2 + 1.0d0 / (s + m2 - 2.0d0)**2
      sum = sum + b2(m) * p * (q * q - r)
23020 continue
      fred = dexp((1.0d0 - s) * loga)
      fred2 = loga**2 * (1.0d0/(s - 1.0d0) + 0.5d0/a)
      sum = (sum + 2.0d0 / (s - 1.0d0)**3 + 2.0d0 * loga / (s - 1.0d0)**
     &2 + fred2) * fred
      do 23022 n = 2,(a - 1) 
      logn = dlog( dfloat(n) )
      sum = sum + (logn)**2 / dexp(logn * s)
23022 continue
      ddzeta8 = sum
      return
      end
      double precision function duds(s, b2)
      double precision s, b2(12), zs, zeta8, dzeta8, ddzeta8
      zs = zeta8(s, b2)
      duds = (dzeta8(s, b2) / zs)**2 - ddzeta8(s, b2) / zs
      return
      end
      subroutine vbecoef(b2)
      double precision b2(12)
      b2(1) = 1.0d0 / 6.0d0
      b2(2) = -1.0d0 / 30.0d0
      b2(3) = 1.0d0 / 42.0d0
      b2(4) = -1.0d0 / 30.0d0
      b2(5) = 5.0d0 / 66.0d0
      b2(6) = -691.0d0 / 2730.0d0
      b2(7) = 7.0d0 / 6.0d0
      b2(8) = -3617.0d0 / 510.0d0
      b2(9) = 4386.7d0 / 79.8d0
      b2(10) = -1746.11d0 / 3.30d0
      b2(11) = 8545.13d0 / 1.38d0
      b2(12) = -2363.64091d0 / 0.02730d0
      return
      end
