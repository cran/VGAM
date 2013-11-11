C Output from Public domain Ratfor, version 1.01
      subroutine vdgam1(x, lfu2qhid, dvhw1ulq)
      implicit logical (a-z)
      double precision x, lfu2qhid
      integer dvhw1ulq
      double precision w, series, obr6tcex
      dvhw1ulq = 1
      if(x .le. 0.0d0)then
      dvhw1ulq = 0
      return
      endif
      if(x .lt. 6.0d0)then
      call vdgam2(x + 6.0d0, obr6tcex, dvhw1ulq)
      lfu2qhid = obr6tcex - 1.0d0/x - 1.0d0/(x + 1.0d0) - 1.0d0/(x + 2.0
     *d0) - 1.0d0/(x + 3.0d0) - 1.0d0/(x + 4.0d0) - 1.0d0/(x + 5.0d0)
      return
      endif
      w = 1.0d0 / (x * x)
      series = ((w * (-1.0d0/12.0d0 + ((w * (1.0d0/120.0d0 + ((w * (-1.0
     *d0/252.0d0 + ((w * (1.0d0/240.0d0 + ((w * (-1.0d0/132.0d0 + ((w * 
     *(691.0d0/32760.0d0 + ((w * (-1.0d0/12.0d0 + (3617.0d0 * w)/8160.0d
     *0)))))))))))))))))))))
      lfu2qhid = ( dlog(x) - 0.5d0/x + series )
      return
      end
      subroutine vdgam2(x, lfu2qhid, dvhw1ulq)
      implicit logical (a-z)
      double precision x, lfu2qhid
      integer dvhw1ulq
      double precision w, series, obr6tcex
      dvhw1ulq = 1
      if(x .le. 0.0d0)then
      dvhw1ulq = 0
      return
      endif
      if(x .lt. 6.0d0)then
      call vdgam1(x + 6.0d0, obr6tcex, dvhw1ulq)
      lfu2qhid = obr6tcex - 1.0d0/x - 1.0d0/(x + 1.0d0) - 1.0d0/(x + 2.0
     *d0) - 1.0d0/(x + 3.0d0) - 1.0d0/(x + 4.0d0) - 1.0d0/(x + 5.0d0)
      return
      endif
      w = 1.0d0 / (x * x)
      series = ((w * (-1.0d0/12.0d0 + ((w * (1.0d0/120.0d0 + ((w * (-1.0
     *d0/252.0d0 + ((w * (1.0d0/240.0d0 + ((w * (-1.0d0/132.0d0 + ((w * 
     *(691.0d0/32760.0d0 + ((w * (-1.0d0/12.0d0 + (3617.0d0 * w)/8160.0d
     *0)))))))))))))))))))))
      lfu2qhid = ( dlog(x) - 0.5d0/x + series )
      return
      end
      subroutine vtgam1(x, lfu2qhid, dvhw1ulq)
      implicit logical (a-z)
      double precision x, lfu2qhid
      integer dvhw1ulq
      double precision w, series, obr6tcex
      dvhw1ulq = 1
      if(x .le. 0.0d0)then
      dvhw1ulq = 0
      return
      endif
      if(x .lt. 6.0d0)then
      call vtgam2(x + 6.0d0, obr6tcex, dvhw1ulq)
      lfu2qhid = obr6tcex + 1.0d0/x**2 + 1.0d0/(x + 1.0d0)**2 + 1.0d0/(x
     * + 2.0d0)**2 + 1.0d0/(x + 3.0d0)**2 + 1.0d0/(x + 4.0d0)**2 + 1.0d0
     */(x + 5.0d0)**2
      return
      endif
      w = 1.0d0 / (x * x)
      series = 1.0d0 + (w * (1.0d0/6.0d0 + (w * (-1.0d0/30.0d0 + (w * (1
     *.0d0/42.0d0 + (w * (-1.0d0/30.0d0 + (w * (5.0d0/66.0d0 + (w * (-69
     *1.0d0/2370.0d0 + (w * (7.0d0/6.0d0 - (3617.0d0 * w)/510.0d0)))))))
     *)))))))
      lfu2qhid = 0.5d0 * w + series / x
      return
      end
      subroutine vtgam2(x, lfu2qhid, dvhw1ulq)
      implicit logical (a-z)
      double precision x, lfu2qhid
      integer dvhw1ulq
      double precision w, series, obr6tcex
      dvhw1ulq = 1
      if(x .le. 0.0d0)then
      dvhw1ulq = 0
      return
      endif
      if(x .lt. 6.0d0)then
      call vtgam1(x + 6.0d0, obr6tcex, dvhw1ulq)
      lfu2qhid = obr6tcex + 1.0d0/x**2 + 1.0d0/(x + 1.0d0)**2 + 1.0d0/(x
     * + 2.0d0)**2 + 1.0d0/(x + 3.0d0)**2 + 1.0d0/(x + 4.0d0)**2 + 1.0d0
     */(x + 5.0d0)**2
      return
      endif
      w = 1.0d0 / (x * x)
      series = 1.0d0 + (w * (1.0d0/6.0d0 + (w * (-1.0d0/30.0d0 + (w * (1
     *.0d0/42.0d0 + (w * (-1.0d0/30.0d0 + (w * (5.0d0/66.0d0 + (w * (-69
     *1.0d0/2370.0d0 + (w * (7.0d0/6.0d0 - (3617.0d0 * w)/510.0d0)))))))
     *)))))))
      lfu2qhid = 0.5d0 * w + series / x
      return
      end
      subroutine dgam1w(x, lfu2qhid, n, dvhw1ulq)
      implicit logical (a-z)
      integer n, dvhw1ulq
      double precision x(n), lfu2qhid(n)
      integer i, okobr6tcex
      dvhw1ulq = 1
      do23016 i=1,n 
      call vdgam1(x(i), lfu2qhid(i), okobr6tcex)
      if(okobr6tcex .ne. 1)then
      dvhw1ulq = okobr6tcex
      endif
23016 continue
23017 continue
      return
      end
      subroutine tgam1w(x, lfu2qhid, n, dvhw1ulq)
      implicit logical (a-z)
      integer n, dvhw1ulq
      double precision x(n), lfu2qhid(n)
      integer i, okobr6tcex
      dvhw1ulq = 1
      do23020 i=1,n 
      call vtgam1(x(i), lfu2qhid(i), okobr6tcex)
      if(okobr6tcex .ne. 1)then
      dvhw1ulq = okobr6tcex
      endif
23020 continue
23021 continue
      return
      end
      subroutine cum8sum(ci1oyxas, lfu2qhid, nlfu2qhid, valong, ntot, no
     *tdvhw1ulq)
      implicit logical (a-z)
      integer nlfu2qhid, ntot, notdvhw1ulq
      double precision ci1oyxas(ntot), lfu2qhid(nlfu2qhid), valong(ntot)
      integer ayfnwr1v, iii
      iii = 1
      lfu2qhid(iii) = ci1oyxas(iii)
      do23024 ayfnwr1v=2,ntot 
      if(valong(ayfnwr1v) .gt. valong(ayfnwr1v-1))then
      lfu2qhid(iii) = lfu2qhid(iii) + ci1oyxas(ayfnwr1v)
      else
      iii = iii + 1
      lfu2qhid(iii) = ci1oyxas(ayfnwr1v)
      endif
23024 continue
23025 continue
      if(iii .eq. nlfu2qhid)then
      notdvhw1ulq = 0
      else
      notdvhw1ulq = 1
      endif
      return
      end
