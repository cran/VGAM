C Output from Public domain Ratfor, version 1.01
      subroutine dpdlyjn(psi, i9mwnvqt, mymu, sigma, kpzavbj3ative, lfu2
     *qhid)
      implicit logical (a-z)
      integer kpzavbj3ative
      double precision psi, i9mwnvqt, mymu, sigma, lfu2qhid(3)
      integer hbsl0gto, izero0
      double precision aa, bb, uqnkc6zg, n3iasxug
      logical cc, pos
      hbsl0gto = 1
      izero0 = 0
      n3iasxug = 1.0d-04
      mymu = 0.0d0
      sigma = 1.0d0
      cc = (psi .ge. 0.0d0)
      if(cc)then
      bb = i9mwnvqt
      pos = (dabs(i9mwnvqt) .le. n3iasxug)
      else
      bb = -2.0d0 + i9mwnvqt
      pos = (dabs(i9mwnvqt-2.0d0) .le. n3iasxug)
      endif
      aa = 1.0d0 + psi * bb
      if(kpzavbj3ative .ge. 0)then
      if(pos)then
      lfu2qhid(1) = psi
      else
      lfu2qhid(1) = aa / bb
      endif
      endif
      if(kpzavbj3ative .ge. 1)then
      if(pos)then
      lfu2qhid(2) = (lfu2qhid(1)**2) / 2
      else
      uqnkc6zg = lfu2qhid(1)
      lfu2qhid(2) = (aa * (dlog(aa)/bb) - uqnkc6zg) / bb
      endif
      endif
      if(kpzavbj3ative .ge. 2)then
      if(pos)then
      lfu2qhid(3) = (lfu2qhid(1)**3) / 3
      else
      uqnkc6zg = lfu2qhid(2) * 2.0d0
      lfu2qhid(3) = (aa * (dlog(aa)/bb) ** 2 - uqnkc6zg) / bb
      endif
      endif
      return
      end
      subroutine gleg11(ghz9vuba, i9mwnvqt, mymu, sigma, kpzavbj3mat, le
     *nkpzavbj3mat, lfu2qhid)
      implicit logical (a-z)
      integer lenkpzavbj3mat
      double precision ghz9vuba, i9mwnvqt, mymu, sigma, kpzavbj3mat(4), 
     *lfu2qhid
      integer hbsl0gto, itwo2, three3
      double precision psi, pim12, o3jyipdf(3), two12
      three3 = 3
      itwo2 = 2
      hbsl0gto = 1
      two12 = 1.41421356237309515d0
      if(lenkpzavbj3mat .gt. 0)then
      lfu2qhid = kpzavbj3mat(4) * (kpzavbj3mat(2)**2 + two12 * sigma * g
     *hz9vuba * kpzavbj3mat(3))
      else
      pim12 = 0.564189583547756279d0
      psi = mymu + two12 * sigma * ghz9vuba
      call dpdlyjn(psi, i9mwnvqt, mymu, sigma, itwo2, o3jyipdf)
      lfu2qhid = (dexp(-ghz9vuba*ghz9vuba) * pim12) * (o3jyipdf(2)**2 + 
     *(psi - mymu) * o3jyipdf(3)) / sigma**2
      endif
      return
      end
      subroutine gleg12(ghz9vuba, i9mwnvqt, mymu, sigma, kpzavbj3mat, le
     *nkpzavbj3mat, lfu2qhid)
      implicit logical (a-z)
      integer lenkpzavbj3mat
      double precision ghz9vuba, i9mwnvqt, mymu, sigma, kpzavbj3mat(4), 
     *lfu2qhid
      integer hbsl0gto, itwo2
      double precision psi, pim12, two12
      double precision tad5vhsu(3)
      itwo2 = 2
      hbsl0gto = 1
      if(lenkpzavbj3mat .gt. 0)then
      lfu2qhid = kpzavbj3mat(4) * (-kpzavbj3mat(2))
      else
      pim12 = 0.564189583547756279d0
      two12 = 1.41421356237309515d0
      psi = mymu + two12 * sigma * ghz9vuba
      call dpdlyjn(psi, i9mwnvqt, mymu, sigma, hbsl0gto, tad5vhsu)
      lfu2qhid = (dexp(-ghz9vuba*ghz9vuba) * pim12) * (-tad5vhsu(2)) / s
     *igma**2
      endif
      return
      end
      subroutine gleg13(ghz9vuba, i9mwnvqt, mymu, sigma, kpzavbj3mat, le
     *nkpzavbj3mat, lfu2qhid)
      implicit logical (a-z)
      integer lenkpzavbj3mat
      double precision ghz9vuba, i9mwnvqt, mymu, sigma, kpzavbj3mat(4), 
     *lfu2qhid
      integer hbsl0gto, itwo2
      double precision psi, mtpim12, two12
      double precision tad5vhsu(3)
      itwo2 = 2
      hbsl0gto = 1
      if(lenkpzavbj3mat .gt. 0)then
      lfu2qhid = kpzavbj3mat(4) * (-kpzavbj3mat(2)) * dsqrt(8.0d0) * ghz
     *9vuba
      else
      mtpim12 = -1.12837916709551256d0
      two12 = 1.41421356237309515d0
      psi = mymu + two12 * sigma * ghz9vuba
      call dpdlyjn(psi, i9mwnvqt, mymu, sigma, hbsl0gto, tad5vhsu)
      lfu2qhid = dexp(-ghz9vuba*ghz9vuba) * mtpim12 * tad5vhsu(2) * (psi
     * - mymu) / sigma**3
      endif
      return
      end
      subroutine gint3(minx, maxx, wts, ahl0onwx, i9mwnvqt, mymu, sigma,
     * kk, lfu2qhid, elemnt)
      implicit logical (a-z)
      integer kk, elemnt
      double precision minx, maxx, wts(kk), ahl0onwx(kk), lfu2qhid, i9mw
     *nvqt, mymu, sigma
      integer gp1jxzuh, lenkpzavbj3mat
      double precision atx, dint, tint, kpzavbj3mat(4), midpt, range12
      lenkpzavbj3mat = 0
      midpt = 0.50d0 * (minx + maxx)
      range12 = 0.50d0 * (maxx - minx)
      dint = 0.0d0
      if(elemnt .eq. 1)then
      do23022 gp1jxzuh=1,kk 
      atx = midpt + range12 * ahl0onwx(gp1jxzuh)
      call gleg11(atx, i9mwnvqt, mymu, sigma, kpzavbj3mat, lenkpzavbj3ma
     *t, tint)
      dint = dint + tint * wts(gp1jxzuh)
23022 continue
23023 continue
      else
      if(elemnt .eq. 2)then
      do23026 gp1jxzuh=1,kk 
      atx = midpt + range12 * ahl0onwx(gp1jxzuh)
      call gleg12(atx, i9mwnvqt, mymu, sigma, kpzavbj3mat, lenkpzavbj3ma
     *t, tint)
      dint = dint + tint * wts(gp1jxzuh)
23026 continue
23027 continue
      else
      if(elemnt .eq. 3)then
      do23030 gp1jxzuh=1,kk 
      atx = midpt + range12 * ahl0onwx(gp1jxzuh)
      call gleg13(atx, i9mwnvqt, mymu, sigma, kpzavbj3mat, lenkpzavbj3ma
     *t, tint)
      dint = dint + tint * wts(gp1jxzuh)
23030 continue
23031 continue
      endif
      endif
      endif
      lfu2qhid = lfu2qhid + range12 * dint
      return
      end
      subroutine yjngintf(minx, maxx, ahl0onwx, wts, kuzxj1lo, kk, i9mwn
     *vqt, mymu, sigma, lfu2qhid, qaltf0nz)
      implicit logical (a-z)
      integer kuzxj1lo, kk
      double precision minx(kuzxj1lo), maxx(kuzxj1lo), wts(kk), ahl0onwx
     *(kk), i9mwnvqt(kuzxj1lo), mymu(kuzxj1lo), sigma(kuzxj1lo), lfu2qhi
     *d(3,kuzxj1lo), qaltf0nz
      integer ayfnwr1v, iii, gp1jxzuh, lencomp, ipzbcvw3, hmayv1xt, elem
     *nt, hbsl0gto, itwo2
      double precision xd4mybgj, j4qgxvlk, wiptsjx8
      hbsl0gto = 1
      itwo2 = 2
      lencomp = 12
      do23032 ayfnwr1v = 1,kuzxj1lo 
      do23034 elemnt=1,3 
      j4qgxvlk = -10.0d0
      do23036 iii=2,lencomp 
      ipzbcvw3 = 2 ** iii
      xd4mybgj = (maxx(ayfnwr1v) - minx(ayfnwr1v)) / ipzbcvw3
      lfu2qhid(elemnt,ayfnwr1v) = 0.0d0
      do23038 gp1jxzuh=1,ipzbcvw3 
      call gint3(minx(ayfnwr1v)+(gp1jxzuh-1)*xd4mybgj, minx(ayfnwr1v)+gp
     *1jxzuh*xd4mybgj, wts, ahl0onwx, i9mwnvqt(ayfnwr1v), mymu(ayfnwr1v)
     *, sigma(ayfnwr1v), kk, lfu2qhid(elemnt,ayfnwr1v), elemnt)
23038 continue
23039 continue
      wiptsjx8 = dabs(lfu2qhid(elemnt,ayfnwr1v) - j4qgxvlk) / (1.0d0 + d
     *abs(lfu2qhid(elemnt,ayfnwr1v)))
      if(wiptsjx8 .lt. qaltf0nz)then
      goto 234
      else
      j4qgxvlk = lfu2qhid(elemnt,ayfnwr1v)
      endif
23036 continue
23037 continue
234   hmayv1xt = 0
23034 continue
23035 continue
23032 continue
23033 continue
      return
      end
