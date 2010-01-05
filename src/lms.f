      subroutine dpdlyjn(psi, dwgkz6, sfnva0, fpqrt7, g8jieq, ghry8z)
      implicit logical (a-z)
      integer g8jieq
      double precision psi, dwgkz6, sfnva0, fpqrt7, ghry8z(3)
      integer uxzze7, oht3ga
      double precision aa, bb, ig5cma, fiumb4
      logical cc, pos
      uxzze7 = 1
      oht3ga = 0
      fiumb4 = 1.0d-04
      sfnva0 = 0.0d0
      fpqrt7 = 1.0d0
      cc = (psi .ge. 0.0d0)
      if(.not.(cc))goto 23000
      bb = dwgkz6
      pos = (dabs(dwgkz6) .le. fiumb4)
      goto 23001
23000 continue
      bb = -2.0d0 + dwgkz6
      pos = (dabs(dwgkz6-2.0d0) .le. fiumb4)
23001 continue
      aa = 1.0d0 + psi * bb
      if(.not.(g8jieq .ge. 0))goto 23002
      if(.not.(pos))goto 23004
      ghry8z(1) = psi
      goto 23005
23004 continue
      ghry8z(1) = aa / bb
23005 continue
23002 continue
      if(.not.(g8jieq .ge. 1))goto 23006
      if(.not.(pos))goto 23008
      ghry8z(2) = (ghry8z(1)**2) / 2
      goto 23009
23008 continue
      ig5cma = ghry8z(1)
      ghry8z(2) = (aa * (dlog(aa)/bb) - ig5cma) / bb
23009 continue
23006 continue
      if(.not.(g8jieq .ge. 2))goto 23010
      if(.not.(pos))goto 23012
      ghry8z(3) = (ghry8z(1)**3) / 3
      goto 23013
23012 continue
      ig5cma = ghry8z(2) * 2.0d0
      ghry8z(3) = (aa * (dlog(aa)/bb) ** 2 - ig5cma) / bb
23013 continue
23010 continue
      return
      end
      subroutine gleg11(hr83e, dwgkz6, sfnva0, fpqrt7, uvf4mp, ws5jcg, 
     &ghry8z)
      implicit logical (a-z)
      integer ws5jcg
      double precision hr83e, dwgkz6, sfnva0, fpqrt7, uvf4mp(4), ghry8z
      integer uxzze7, itwo2, ynmzp6
      double precision psi, hc0tub, xkwp2m(3), dq3rxy
      ynmzp6 = 3
      itwo2 = 2
      uxzze7 = 1
      dq3rxy = 1.41421356237309515d0
      if(.not.(ws5jcg .gt. 0))goto 23014
      ghry8z = uvf4mp(4) * (uvf4mp(2)**2 + dq3rxy * fpqrt7 * hr83e * 
     &uvf4mp(3))
      goto 23015
23014 continue
      hc0tub = 0.564189583547756279d0
      psi = sfnva0 + dq3rxy * fpqrt7 * hr83e
      call dpdlyjn(psi, dwgkz6, sfnva0, fpqrt7, itwo2, xkwp2m)
      ghry8z = (dexp(-hr83e*hr83e) * hc0tub) * (xkwp2m(2)**2 + (psi - 
     &sfnva0) * xkwp2m(3)) / fpqrt7**2
23015 continue
      return
      end
      subroutine zuqx1p(hr83e, dwgkz6, sfnva0, fpqrt7, uvf4mp, ws5jcg, 
     &ghry8z)
      implicit logical (a-z)
      integer ws5jcg
      double precision hr83e, dwgkz6, sfnva0, fpqrt7, uvf4mp(4), ghry8z
      integer uxzze7, itwo2
      double precision psi, hc0tub, mw6reg(2), dq3rxy
      itwo2 = 2
      uxzze7 = 1
      if(.not.(ws5jcg .gt. 0))goto 23016
      ghry8z = uvf4mp(4) * (-uvf4mp(2))
      goto 23017
23016 continue
      hc0tub = 0.564189583547756279d0
      dq3rxy = 1.41421356237309515d0
      psi = sfnva0 + dq3rxy * fpqrt7 * hr83e
      call dpdlyjn(psi, dwgkz6, sfnva0, fpqrt7, uxzze7, mw6reg)
      ghry8z = (dexp(-hr83e*hr83e) * hc0tub) * (-mw6reg(2)) / fpqrt7**2
23017 continue
      return
      end
      subroutine gleg13(hr83e, dwgkz6, sfnva0, fpqrt7, uvf4mp, ws5jcg, 
     &ghry8z)
      implicit logical (a-z)
      integer ws5jcg
      double precision hr83e, dwgkz6, sfnva0, fpqrt7, uvf4mp(4), ghry8z
      integer uxzze7, itwo2
      double precision psi, oaqng6, mw6reg(2), dq3rxy
      itwo2 = 2
      uxzze7 = 1
      if(.not.(ws5jcg .gt. 0))goto 23018
      ghry8z = uvf4mp(4) * (-uvf4mp(2)) * dsqrt(8.0d0) * hr83e
      goto 23019
23018 continue
      oaqng6 = -1.12837916709551256d0
      dq3rxy = 1.41421356237309515d0
      psi = sfnva0 + dq3rxy * fpqrt7 * hr83e
      call dpdlyjn(psi, dwgkz6, sfnva0, fpqrt7, uxzze7, mw6reg)
      ghry8z = dexp(-hr83e*hr83e) * oaqng6 * mw6reg(2) * (psi - sfnva0) 
     &/ fpqrt7**3
23019 continue
      return
      end
      subroutine rnvz5t(r7zvis, bd8olv, wts, oqie8v, dwgkz6, sfnva0, 
     &fpqrt7, kk, ghry8z, nepms8)
      implicit logical (a-z)
      integer kk, nepms8
      double precision r7zvis, bd8olv, wts(kk), oqie8v(kk), ghry8z, 
     &dwgkz6, sfnva0, fpqrt7
      integer nd6mep, ws5jcg
      double precision atx, tns0gf, dy3ljx, uvf4mp(4), byn1gh, k8ousd
      ws5jcg = 0
      byn1gh = 0.50d0 * (r7zvis + bd8olv)
      k8ousd = 0.50d0 * (bd8olv - r7zvis)
      tns0gf = 0.0d0
      if(.not.(nepms8 .eq. 1))goto 23020
      do 23022 nd6mep=1,kk 
      atx = byn1gh + k8ousd * oqie8v(nd6mep)
      call gleg11(atx, dwgkz6, sfnva0, fpqrt7, uvf4mp, ws5jcg, dy3ljx)
      tns0gf = tns0gf + dy3ljx * wts(nd6mep)
23022 continue
      goto 23021
23020 continue
      if(.not.(nepms8 .eq. 2))goto 23024
      do 23026 nd6mep=1,kk 
      atx = byn1gh + k8ousd * oqie8v(nd6mep)
      call zuqx1p(atx, dwgkz6, sfnva0, fpqrt7, uvf4mp, ws5jcg, dy3ljx)
      tns0gf = tns0gf + dy3ljx * wts(nd6mep)
23026 continue
      goto 23025
23024 continue
      if(.not.(nepms8 .eq. 3))goto 23028
      do 23030 nd6mep=1,kk 
      atx = byn1gh + k8ousd * oqie8v(nd6mep)
      call gleg13(atx, dwgkz6, sfnva0, fpqrt7, uvf4mp, ws5jcg, dy3ljx)
      tns0gf = tns0gf + dy3ljx * wts(nd6mep)
23030 continue
23028 continue
23025 continue
23021 continue
      ghry8z = ghry8z + k8ousd * tns0gf
      return
      end
      subroutine yjngintf(r7zvis, bd8olv, oqie8v, wts, nfiumb4, kk, 
     &dwgkz6, sfnva0, fpqrt7, ghry8z, kqoy6w)
      implicit logical (a-z)
      integer nfiumb4, kk
      double precision r7zvis(nfiumb4), bd8olv(nfiumb4), wts(kk), 
     &oqie8v(kk), dwgkz6(nfiumb4), sfnva0(nfiumb4), fpqrt7(nfiumb4), 
     &ghry8z(3,nfiumb4), kqoy6w
      integer w3gohz, p1rifj, nd6mep, o2yadh, btip7u, epx9jf, nepms8, 
     &uxzze7, itwo2
      double precision mu4ygk, azgts7, xmr7cj
      uxzze7 = 1
      itwo2 = 2
      o2yadh = 12
      do 23032 w3gohz = 1,nfiumb4 
      do 23034 nepms8=1,3 
      azgts7 = -10.0d0
      do 23036 p1rifj=2,o2yadh 
      btip7u = 2 ** p1rifj
      mu4ygk = (bd8olv(w3gohz) - r7zvis(w3gohz)) / btip7u
      ghry8z(nepms8,w3gohz) = 0.0d0
      do 23038 nd6mep=1,btip7u 
      call rnvz5t(r7zvis(w3gohz)+(nd6mep-1)*mu4ygk, r7zvis(w3gohz)+
     &nd6mep*mu4ygk, wts, oqie8v, dwgkz6(w3gohz), sfnva0(w3gohz), 
     &fpqrt7(w3gohz), kk, ghry8z(nepms8,w3gohz), nepms8)
23038 continue
      xmr7cj = dabs(ghry8z(nepms8,w3gohz) - azgts7) / (1.0d0 + dabs(
     &ghry8z(nepms8,w3gohz)))
      if(.not.(xmr7cj .lt. kqoy6w))goto 23040
      goto 234
      goto 23041
23040 continue
      azgts7 = ghry8z(nepms8,w3gohz)
23041 continue
23036 continue
234   epx9jf = 0
23034 continue
23032 continue
      return
      end
