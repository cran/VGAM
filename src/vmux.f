C Output from Public domain Ratfor, version 1.01
      subroutine qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      implicit logical (a-z)
      integer wy1vqfzu, tgiyxdw1(*), dufozmt7(*)
      integer urohxe6t, bpvaqm5z, ayfnwr1v
      ayfnwr1v = 1
      urohxe6t = wy1vqfzu
23000 if(.not.(urohxe6t .ge. 1))goto 23002
      do23003 bpvaqm5z=1,urohxe6t 
      tgiyxdw1(ayfnwr1v) = bpvaqm5z
      ayfnwr1v = ayfnwr1v+1
23003 continue
23004 continue
23001 urohxe6t=urohxe6t-1
      goto 23000
23002 continue
      ayfnwr1v = 1
      do23005 urohxe6t=1,wy1vqfzu 
      do23007 bpvaqm5z=urohxe6t,wy1vqfzu 
      dufozmt7(ayfnwr1v) = bpvaqm5z
      ayfnwr1v = ayfnwr1v+1
23007 continue
23008 continue
23005 continue
23006 continue
      return
      end
      integer function viamf(cz8qdfyj, rvy1fpli, wy1vqfzu, tgiyxdw1, duf
     *ozmt7)
      integer cz8qdfyj, rvy1fpli, wy1vqfzu, tgiyxdw1(*), dufozmt7(*)
      integer urohxe6t, imk5wjxg
      imk5wjxg = wy1vqfzu*(wy1vqfzu+1)/2
      do23009 urohxe6t=1,imk5wjxg 
      if((tgiyxdw1(urohxe6t).eq.cz8qdfyj .and. dufozmt7(urohxe6t).eq.rvy
     *1fpli) .or. (tgiyxdw1(urohxe6t).eq.rvy1fpli .and. dufozmt7(urohxe6
     *t).eq.cz8qdfyj))then
      viamf = urohxe6t
      return
      endif
23009 continue
23010 continue
      viamf = 0
      return
      end
      subroutine vm2af(mat, a, dimm, tgiyxdw1, dufozmt7, kuzxj1lo, wy1vq
     *fzu, rb1onzwu)
      implicit logical (a-z)
      integer dimm, tgiyxdw1(dimm), dufozmt7(dimm), kuzxj1lo, wy1vqfzu, 
     *rb1onzwu
      double precision mat(dimm,kuzxj1lo), a(wy1vqfzu,wy1vqfzu,kuzxj1lo)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, imk5wjxg
      imk5wjxg = wy1vqfzu * (wy1vqfzu + 1) / 2
      if(rb1onzwu .eq. 1 .or. dimm .ne. imk5wjxg)then
      ayfnwr1v = 1
23015 if(.not.(ayfnwr1v .le. kuzxj1lo))goto 23017
      yq6lorbx = 1
23018 if(.not.(yq6lorbx .le. wy1vqfzu))goto 23020
      gp1jxzuh = 1
23021 if(.not.(gp1jxzuh .le. wy1vqfzu))goto 23023
      a(gp1jxzuh,yq6lorbx,ayfnwr1v) = 0.0d0
23022 gp1jxzuh=gp1jxzuh+1
      goto 23021
23023 continue
23019 yq6lorbx=yq6lorbx+1
      goto 23018
23020 continue
23016 ayfnwr1v=ayfnwr1v+1
      goto 23015
23017 continue
      endif
      do23024 ayfnwr1v=1,kuzxj1lo 
      do23026 yq6lorbx=1,dimm 
      a(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx),ayfnwr1v) = mat(yq6lorbx,a
     *yfnwr1v)
      if(rb1onzwu .eq. 0)then
      a(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx),ayfnwr1v) = mat(yq6lorbx,a
     *yfnwr1v)
      endif
23026 continue
23027 continue
23024 continue
23025 continue
      return
      end
      subroutine mux22f(wpuarq2m, tlgduey8, lfu2qhid, dimu, tgiyxdw1, du
     *fozmt7, kuzxj1lo, wy1vqfzu, wk1200)
      implicit logical (a-z)
      integer dimu, tgiyxdw1(*), dufozmt7(*), kuzxj1lo, wy1vqfzu
      double precision wpuarq2m(dimu,kuzxj1lo), tlgduey8(kuzxj1lo,wy1vqf
     *zu), lfu2qhid(wy1vqfzu,kuzxj1lo), wk1200(wy1vqfzu,wy1vqfzu)
      double precision q6zdcwxk
      integer ayfnwr1v, yq6lorbx, bpvaqm5z, one, rb1onzwu
      one = 1
      rb1onzwu = 1
      ayfnwr1v = 1
23030 if(.not.(ayfnwr1v .le. kuzxj1lo))goto 23032
      call vm2af(wpuarq2m(1,ayfnwr1v), wk1200, dimu, tgiyxdw1, dufozmt7,
     * one, wy1vqfzu, rb1onzwu)
      yq6lorbx = 1
23033 if(.not.(yq6lorbx .le. wy1vqfzu))goto 23035
      q6zdcwxk = 0.0d0
      bpvaqm5z = yq6lorbx
23036 if(.not.(bpvaqm5z .le. wy1vqfzu))goto 23038
      q6zdcwxk = q6zdcwxk + wk1200(yq6lorbx,bpvaqm5z) * tlgduey8(ayfnwr1
     *v,bpvaqm5z)
23037 bpvaqm5z=bpvaqm5z+1
      goto 23036
23038 continue
      lfu2qhid(yq6lorbx,ayfnwr1v) = q6zdcwxk
23034 yq6lorbx=yq6lorbx+1
      goto 23033
23035 continue
23031 ayfnwr1v=ayfnwr1v+1
      goto 23030
23032 continue
      return
      end
      subroutine vbksf(wpuarq2m, bvecto, wy1vqfzu, kuzxj1lo, wk1200, tgi
     *yxdw1, dufozmt7, dimu)
      implicit logical (a-z)
      integer wy1vqfzu, kuzxj1lo, tgiyxdw1(*), dufozmt7(*), dimu
      double precision wpuarq2m(dimu,kuzxj1lo), bvecto(wy1vqfzu,kuzxj1lo
     *), wk1200(wy1vqfzu,wy1vqfzu)
      double precision q6zdcwxk
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, rb1onzwu, one
      rb1onzwu = 1
      one = 1
      ayfnwr1v = 1
23039 if(.not.(ayfnwr1v .le. kuzxj1lo))goto 23041
      call vm2af(wpuarq2m(1,ayfnwr1v), wk1200, dimu, tgiyxdw1, dufozmt7,
     * one, wy1vqfzu, rb1onzwu)
      yq6lorbx = wy1vqfzu
23042 if(.not.(yq6lorbx .ge. 1))goto 23044
      q6zdcwxk = bvecto(yq6lorbx,ayfnwr1v)
      gp1jxzuh = yq6lorbx+1
23045 if(.not.(gp1jxzuh .le. wy1vqfzu))goto 23047
      q6zdcwxk = q6zdcwxk - wk1200(yq6lorbx,gp1jxzuh) * bvecto(gp1jxzuh,
     *ayfnwr1v)
23046 gp1jxzuh=gp1jxzuh+1
      goto 23045
23047 continue
      bvecto(yq6lorbx,ayfnwr1v) = q6zdcwxk / wk1200(yq6lorbx,yq6lorbx)
23043 yq6lorbx=yq6lorbx-1
      goto 23042
23044 continue
23040 ayfnwr1v=ayfnwr1v+1
      goto 23039
23041 continue
      return
      end
      subroutine vcholf(wmat, bvecto, wy1vqfzu, dvhw1ulq, isolve)
      implicit logical (a-z)
      integer isolve
      integer wy1vqfzu, dvhw1ulq
      double precision wmat(wy1vqfzu,wy1vqfzu), bvecto(wy1vqfzu)
      double precision q6zdcwxk, dsqrt
      integer ayfnwr1v, yq6lorbx, gp1jxzuh
      dvhw1ulq=1
      do23048 ayfnwr1v=1,wy1vqfzu
      q6zdcwxk = 0d0
      do23050 gp1jxzuh=1,ayfnwr1v-1 
      q6zdcwxk = q6zdcwxk + wmat(gp1jxzuh,ayfnwr1v) * wmat(gp1jxzuh,ayfn
     *wr1v)
23050 continue
23051 continue
      wmat(ayfnwr1v,ayfnwr1v) = wmat(ayfnwr1v,ayfnwr1v) - q6zdcwxk
      if(wmat(ayfnwr1v,ayfnwr1v) .le. 0d0)then
      dvhw1ulq = 0
      return
      endif
      wmat(ayfnwr1v,ayfnwr1v) = dsqrt(wmat(ayfnwr1v,ayfnwr1v))
      do23054 yq6lorbx=ayfnwr1v+1,wy1vqfzu
      q6zdcwxk = 0d0
      do23056 gp1jxzuh=1,ayfnwr1v-1 
      q6zdcwxk = q6zdcwxk + wmat(gp1jxzuh,ayfnwr1v) * wmat(gp1jxzuh,yq6l
     *orbx)
23056 continue
23057 continue
      wmat(ayfnwr1v,yq6lorbx) = (wmat(ayfnwr1v,yq6lorbx) - q6zdcwxk) / w
     *mat(ayfnwr1v,ayfnwr1v)
23054 continue
23055 continue
23048 continue
23049 continue
      if(isolve .eq. 0)then
      do23060 ayfnwr1v=2,wy1vqfzu 
      do23062 yq6lorbx=1,ayfnwr1v-1 
      wmat(ayfnwr1v,yq6lorbx) = 0.0d0
23062 continue
23063 continue
      return
23060 continue
23061 continue
      endif
      do23064 yq6lorbx=1,wy1vqfzu 
      q6zdcwxk = bvecto(yq6lorbx)
      do23066 gp1jxzuh=1,yq6lorbx-1 
      q6zdcwxk = q6zdcwxk - wmat(gp1jxzuh,yq6lorbx) * bvecto(gp1jxzuh)
23066 continue
23067 continue
      bvecto(yq6lorbx) = q6zdcwxk / wmat(yq6lorbx,yq6lorbx)
23064 continue
23065 continue
      yq6lorbx = wy1vqfzu
23068 if(.not.(yq6lorbx .ge. 1))goto 23070
      q6zdcwxk = bvecto(yq6lorbx)
      gp1jxzuh = yq6lorbx+1
23071 if(.not.(gp1jxzuh .le. wy1vqfzu))goto 23073
      q6zdcwxk = q6zdcwxk - wmat(yq6lorbx,gp1jxzuh) * bvecto(gp1jxzuh)
23072 gp1jxzuh=gp1jxzuh+1
      goto 23071
23073 continue
      bvecto(yq6lorbx) = q6zdcwxk / wmat(yq6lorbx,yq6lorbx)
23069 yq6lorbx=yq6lorbx-1
      goto 23068
23070 continue
      return
      end
      subroutine mux17f(wpuarq2m, he7mqnvy, wy1vqfzu, xjc4ywlh, kuzxj1lo
     *, wk1200, wk3400, tgiyxdw1, dufozmt7, dimu, rutyk8mg)
      implicit logical (a-z)
      integer dimu, wy1vqfzu, xjc4ywlh, kuzxj1lo, tgiyxdw1(*), dufozmt7(
     **), rutyk8mg
      double precision wpuarq2m(dimu,kuzxj1lo), he7mqnvy(rutyk8mg,xjc4yw
     *lh), wk1200(wy1vqfzu,wy1vqfzu), wk3400(wy1vqfzu,xjc4ywlh)
      double precision q6zdcwxk
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, bpvaqm5z
      do23074 yq6lorbx=1,wy1vqfzu 
      do23076 ayfnwr1v=1,wy1vqfzu 
      wk1200(ayfnwr1v,yq6lorbx) = 0.0d0
23076 continue
23077 continue
23074 continue
23075 continue
      do23078 ayfnwr1v=1,kuzxj1lo 
      do23080 bpvaqm5z=1,dimu 
      wk1200(tgiyxdw1(bpvaqm5z), dufozmt7(bpvaqm5z)) = wpuarq2m(bpvaqm5z
     *,ayfnwr1v)
23080 continue
23081 continue
      do23082 gp1jxzuh=1,xjc4ywlh 
      do23084 yq6lorbx=1,wy1vqfzu 
      wk3400(yq6lorbx,gp1jxzuh) = he7mqnvy((ayfnwr1v-1)*wy1vqfzu+yq6lorb
     *x,gp1jxzuh)
23084 continue
23085 continue
23082 continue
23083 continue
      do23086 gp1jxzuh=1,xjc4ywlh 
      do23088 yq6lorbx=1,wy1vqfzu 
      q6zdcwxk = 0d0
      do23090 bpvaqm5z=yq6lorbx,wy1vqfzu 
      q6zdcwxk = q6zdcwxk + wk1200(yq6lorbx,bpvaqm5z) * wk3400(bpvaqm5z,
     *gp1jxzuh)
23090 continue
23091 continue
      he7mqnvy((ayfnwr1v-1)*wy1vqfzu+yq6lorbx,gp1jxzuh) = q6zdcwxk
23088 continue
23089 continue
23086 continue
23087 continue
23078 continue
23079 continue
      return
      end
      subroutine vrinvf9(wpuarq2m, ldr, wy1vqfzu, dvhw1ulq, ks3wejcv, wo
     *rk)
      implicit logical (a-z)
      integer ldr, wy1vqfzu, dvhw1ulq
      double precision wpuarq2m(ldr,wy1vqfzu), ks3wejcv(wy1vqfzu,wy1vqfz
     *u), work(wy1vqfzu,wy1vqfzu)
      double precision q6zdcwxk
      integer yq6lorbx, gp1jxzuh, col, uaoynef0
      dvhw1ulq = 1
      yq6lorbx = 1
23092 if(.not.(yq6lorbx .le. wy1vqfzu))goto 23094
      col = 1
23095 if(.not.(col .le. wy1vqfzu))goto 23097
      work(yq6lorbx,col) = 0.0d0
23096 col=col+1
      goto 23095
23097 continue
23093 yq6lorbx=yq6lorbx+1
      goto 23092
23094 continue
      col = 1
23098 if(.not.(col .le. wy1vqfzu))goto 23100
      yq6lorbx = col
23101 if(.not.(yq6lorbx .ge. 1))goto 23103
      if(yq6lorbx .eq. col)then
      q6zdcwxk = 1.0d0
      else
      q6zdcwxk = 0.0d0
      endif
      gp1jxzuh = yq6lorbx+1
23106 if(.not.(gp1jxzuh .le. col))goto 23108
      q6zdcwxk = q6zdcwxk - wpuarq2m(yq6lorbx,gp1jxzuh) * work(gp1jxzuh,
     *col)
23107 gp1jxzuh=gp1jxzuh+1
      goto 23106
23108 continue
      if(wpuarq2m(yq6lorbx,yq6lorbx) .eq. 0.0d0)then
      dvhw1ulq = 0
      else
      work(yq6lorbx,col) = q6zdcwxk / wpuarq2m(yq6lorbx,yq6lorbx)
      endif
23102 yq6lorbx=yq6lorbx-1
      goto 23101
23103 continue
23099 col=col+1
      goto 23098
23100 continue
      yq6lorbx = 1
23111 if(.not.(yq6lorbx .le. wy1vqfzu))goto 23113
      col = yq6lorbx
23114 if(.not.(col .le. wy1vqfzu))goto 23116
      if(yq6lorbx .lt. col)then
      uaoynef0 = col
      else
      uaoynef0 = yq6lorbx
      endif
      q6zdcwxk = 0.0d0
      gp1jxzuh = uaoynef0
23119 if(.not.(gp1jxzuh .le. wy1vqfzu))goto 23121
      q6zdcwxk = q6zdcwxk + work(yq6lorbx,gp1jxzuh) * work(col,gp1jxzuh)
23120 gp1jxzuh=gp1jxzuh+1
      goto 23119
23121 continue
      ks3wejcv(yq6lorbx,col) = q6zdcwxk
      ks3wejcv(col,yq6lorbx) = q6zdcwxk
23115 col=col+1
      goto 23114
23116 continue
23112 yq6lorbx=yq6lorbx+1
      goto 23111
23113 continue
      return
      end
      subroutine tldz5ion(xx, lfu2qhid)
      implicit logical (a-z)
      double precision xx, lfu2qhid
      double precision x, y, hofjnx2e, q6zdcwxk, xd4mybgj(6)
      integer yq6lorbx
      xd4mybgj(1)= 76.18009172947146d0
      xd4mybgj(2)= -86.50532032941677d0
      xd4mybgj(3)= 24.01409824083091d0
      xd4mybgj(4)= -1.231739572450155d0
      xd4mybgj(5)= 0.1208650973866179d-2
      xd4mybgj(6)= -0.5395239384953d-5
      x = xx
      y = xx
      hofjnx2e = x+5.50d0
      hofjnx2e = hofjnx2e - (x+0.50d0) * dlog(hofjnx2e)
      q6zdcwxk=1.000000000190015d0
      yq6lorbx=1
23122 if(.not.(yq6lorbx .le. 6))goto 23124
      y = y + 1.0d0
      q6zdcwxk = q6zdcwxk + xd4mybgj(yq6lorbx)/y
23123 yq6lorbx=yq6lorbx+1
      goto 23122
23124 continue
      lfu2qhid = -hofjnx2e + dlog(2.5066282746310005d0 * q6zdcwxk / x)
      return
      end
      subroutine enbin9(bzmd6ftv, hdqsx7bk, nm0eljqk, n2kersmx, n, dvhw1
     *ulq, zy1mchbf, ux3nadiw, rsynp1go, sguwj9ty)
      implicit logical (a-z)
      integer n, dvhw1ulq, zy1mchbf, sguwj9ty
      double precision bzmd6ftv(n, zy1mchbf), hdqsx7bk(n, zy1mchbf), nm0
     *eljqk(n, zy1mchbf), n2kersmx, ux3nadiw, rsynp1go
      integer ayfnwr1v, kij0gwer
      double precision oxjgzv0e, btiehdm2, ydb, vjz5sxty, esql7umk, pvcj
     *l2na, mwuvskg1, ft3ijqmy, hmayv1xt, q6zdcwxk, plo6hkdr
      real csi9ydge
      if(n2kersmx .le. 0.80d0 .or. n2kersmx .ge. 1.0d0)then
      dvhw1ulq = 0
      return
      endif
      btiehdm2 = 100.0d0 * rsynp1go
      oxjgzv0e = 0.001d0
      dvhw1ulq = 1
      kij0gwer=1
23127 if(.not.(kij0gwer.le.zy1mchbf))goto 23129
      ayfnwr1v=1
23130 if(.not.(ayfnwr1v.le.n))goto 23132
      vjz5sxty = nm0eljqk(ayfnwr1v,kij0gwer) / hdqsx7bk(ayfnwr1v,kij0gwe
     *r)
      if((vjz5sxty .lt. oxjgzv0e) .or. (nm0eljqk(ayfnwr1v,kij0gwer) .gt.
     * 1.0d5))then
      bzmd6ftv(ayfnwr1v,kij0gwer) = -nm0eljqk(ayfnwr1v,kij0gwer) * (1.0d
     *0 + hdqsx7bk(ayfnwr1v,kij0gwer)/(hdqsx7bk(ayfnwr1v,kij0gwer) + nm0
     *eljqk(ayfnwr1v,kij0gwer))) / hdqsx7bk(ayfnwr1v,kij0gwer)**2
      if(bzmd6ftv(ayfnwr1v,kij0gwer) .gt. -btiehdm2)then
      bzmd6ftv(ayfnwr1v,kij0gwer) = -btiehdm2
      endif
      goto 20
      endif
      q6zdcwxk = 0.0d0
      pvcjl2na = hdqsx7bk(ayfnwr1v,kij0gwer) / (hdqsx7bk(ayfnwr1v,kij0gw
     *er) + nm0eljqk(ayfnwr1v,kij0gwer))
      mwuvskg1 = 1.0d0 - pvcjl2na
      csi9ydge = hdqsx7bk(ayfnwr1v,kij0gwer)
      if(pvcjl2na .lt. btiehdm2)then
      pvcjl2na = btiehdm2
      endif
      if(mwuvskg1 .lt. btiehdm2)then
      mwuvskg1 = btiehdm2
      endif
      esql7umk = 100.0d0 + 15.0d0 * nm0eljqk(ayfnwr1v,kij0gwer)
      if(esql7umk .lt. sguwj9ty)then
      esql7umk = sguwj9ty
      endif
      ft3ijqmy = pvcjl2na ** csi9ydge
      ux3nadiw = ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / hdqsx7bk(ayfnwr1v,kij0gwer)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = 1.0d0
      ft3ijqmy = hdqsx7bk(ayfnwr1v,kij0gwer) * mwuvskg1 * ft3ijqmy
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + ydb
     *)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = 2.0d0
23143 if(((ux3nadiw .le. n2kersmx) .or. (plo6hkdr .gt. 1.0d-4)) .and. (y
     *db .lt. esql7umk))then
      ft3ijqmy = (hdqsx7bk(ayfnwr1v,kij0gwer) - 1.0d0 + ydb) * mwuvskg1 
     ** ft3ijqmy / ydb
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + ydb
     *)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = ydb + 1.0d0
      goto 23143
      endif
23144 continue
      bzmd6ftv(ayfnwr1v,kij0gwer) = -q6zdcwxk
20    hmayv1xt = 0.0d0
23131 ayfnwr1v=ayfnwr1v+1
      goto 23130
23132 continue
23128 kij0gwer=kij0gwer+1
      goto 23127
23129 continue
      return
      end
      subroutine enbin8(bzmd6ftv, hdqsx7bk, hsj9bzaq, n2kersmx, kuzxj1lo
     *, dvhw1ulq, zy1mchbf, ux3nadiw, rsynp1go)
      implicit logical (a-z)
      integer kuzxj1lo, dvhw1ulq, zy1mchbf
      double precision bzmd6ftv(kuzxj1lo, zy1mchbf), hdqsx7bk(kuzxj1lo, 
     *zy1mchbf), hsj9bzaq(kuzxj1lo, zy1mchbf), n2kersmx, ux3nadiw, rsynp
     *1go
      integer ayfnwr1v, kij0gwer, esql7umk
      double precision ft3ijqmy, tad5vhsu, o3jyipdf, pq0hfucn, q6zdcwxk,
     * d1, d2, plo6hkdr, hnu1vjyw
      logical pok1, pok2, pok12
      double precision oxjgzv0e, onemse, nm0eljqk, btiehdm2, ydb, kbig
      d1 = 0.0d0
      d2 = 0.0d0
      btiehdm2 = -100.0d0 * rsynp1go
      esql7umk = 3000
      if(n2kersmx .le. 0.80d0 .or. n2kersmx .ge. 1.0d0)then
      dvhw1ulq = 0
      return
      endif
      kbig = 1.0d4
      oxjgzv0e = 0.001d0
      hnu1vjyw = 1.0d0 - rsynp1go
      onemse = 1.0d0 / (1.0d0 + oxjgzv0e)
      dvhw1ulq = 1
      kij0gwer=1
23147 if(.not.(kij0gwer.le.zy1mchbf))goto 23149
      ayfnwr1v=1
23150 if(.not.(ayfnwr1v.le.kuzxj1lo))goto 23152
      if(hdqsx7bk(ayfnwr1v,kij0gwer) .gt. kbig)then
      hdqsx7bk(ayfnwr1v,kij0gwer) = kbig
      endif
      if(hsj9bzaq(ayfnwr1v,kij0gwer) .lt. oxjgzv0e)then
      hsj9bzaq(ayfnwr1v,kij0gwer) = oxjgzv0e
      endif
      if((hsj9bzaq(ayfnwr1v,kij0gwer) .gt. onemse))then
      nm0eljqk = hdqsx7bk(ayfnwr1v,kij0gwer) * (1.0d0/hsj9bzaq(ayfnwr1v,
     *kij0gwer) - 1.0d0)
      bzmd6ftv(ayfnwr1v,kij0gwer) = -nm0eljqk * (1.0d0 + hdqsx7bk(ayfnwr
     *1v,kij0gwer)/(hdqsx7bk(ayfnwr1v,kij0gwer) + nm0eljqk)) / hdqsx7bk(
     *ayfnwr1v,kij0gwer)**2
      if(bzmd6ftv(ayfnwr1v,kij0gwer) .gt. btiehdm2)then
      bzmd6ftv(ayfnwr1v,kij0gwer) = btiehdm2
      endif
      goto 20
      endif
      q6zdcwxk = 0.0d0
      pok1 = .true.
      pok2 = hsj9bzaq(ayfnwr1v,kij0gwer) .lt. (1.0d0-rsynp1go)
      pok12 = pok1 .and. pok2
      if(pok12)then
      d2 = hdqsx7bk(ayfnwr1v,kij0gwer) * dlog(hsj9bzaq(ayfnwr1v,kij0gwer
     *))
      ux3nadiw = dexp(d2)
      else
      ux3nadiw = 0.0d0
      endif
      plo6hkdr = (1.0d0 - ux3nadiw) / hdqsx7bk(ayfnwr1v,kij0gwer)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      call tldz5ion(hdqsx7bk(ayfnwr1v,kij0gwer), o3jyipdf)
      ydb = 1.0d0
      call tldz5ion(ydb + hdqsx7bk(ayfnwr1v,kij0gwer), tad5vhsu)
      pq0hfucn = 0.0d0
      if(pok12)then
      d1 = dlog(1.0d0 - hsj9bzaq(ayfnwr1v,kij0gwer))
      ft3ijqmy = dexp(ydb * d1 + d2 + tad5vhsu - o3jyipdf - pq0hfucn)
      else
      ft3ijqmy = 0.0d0
      endif
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + ydb
     *)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = 2.0d0
23165 if((ux3nadiw .le. n2kersmx) .or. (plo6hkdr .gt. 1.0d-4))then
      tad5vhsu = tad5vhsu + dlog(ydb + hdqsx7bk(ayfnwr1v,kij0gwer) - 1.0
     *d0)
      pq0hfucn = pq0hfucn + dlog(ydb)
      if(pok12)then
      ft3ijqmy = dexp(ydb * d1 + d2 + tad5vhsu - o3jyipdf - pq0hfucn)
      else
      ft3ijqmy = 0.0d0
      endif
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + ydb
     *)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = ydb + 1.0d0
      if(ydb .gt. 1.0d3)then
      goto 21
      endif
      goto 23165
      endif
23166 continue
21    bzmd6ftv(ayfnwr1v,kij0gwer) = -q6zdcwxk
20    tad5vhsu = 0.0d0
23151 ayfnwr1v=ayfnwr1v+1
      goto 23150
23152 continue
23148 kij0gwer=kij0gwer+1
      goto 23147
23149 continue
      return
      end
      subroutine mbessi0(bvecto, kuzxj1lo, kpzavbj3, d0, d1, d2, zjkrtol
     *8, qaltf0nz)
      implicit logical (a-z)
      integer kuzxj1lo, kpzavbj3, zjkrtol8, c5aesxkus
      double precision bvecto(kuzxj1lo), d0(kuzxj1lo), d1(kuzxj1lo), d2(
     *kuzxj1lo), qaltf0nz
      integer ayfnwr1v, gp1jxzuh
      double precision f0, t0, m0, f1, t1, m1, f2, t2, m2
      double precision toobig
      toobig = 20.0d0
      zjkrtol8 = 0
      if(.not.(kpzavbj3 .eq. 0 .or. kpzavbj3 .eq. 1 .or. kpzavbj3 .eq. 2
     *))then
      zjkrtol8 = 1
      return
      endif
      do23173 gp1jxzuh=1,kuzxj1lo 
      if(dabs(bvecto(gp1jxzuh)) .gt. toobig)then
      zjkrtol8 = 1
      return
      endif
      t1 = bvecto(gp1jxzuh) / 2.0d0
      f1 = t1
      t0 = t1 * t1
      f0 = 1.0d0 + t0
      t2 = 0.50d0
      f2 = t2
      c5aesxkus = 15
      if(dabs(bvecto(gp1jxzuh)) .gt. 10)then
      c5aesxkus = 25
      endif
      if(dabs(bvecto(gp1jxzuh)) .gt. 15)then
      c5aesxkus = 35
      endif
      if(dabs(bvecto(gp1jxzuh)) .gt. 20)then
      c5aesxkus = 40
      endif
      if(dabs(bvecto(gp1jxzuh)) .gt. 30)then
      c5aesxkus = 55
      endif
      do23185 ayfnwr1v=1,c5aesxkus 
      m0 = (bvecto(gp1jxzuh) / (2.0d0*(ayfnwr1v+1.0d0))) ** 2.0
      m1 = m0 * (1.0d0 + 1.0d0/ayfnwr1v)
      m2 = m1 * (2.0d0*ayfnwr1v + 1.0d0) / (2.0d0*ayfnwr1v - 1.0d0)
      t0 = t0 * m0
      t1 = t1 * m1
      t2 = t2 * m2
      f0 = f0 + t0
      f1 = f1 + t1
      f2 = f2 + t2
      if((dabs(t0) .lt. qaltf0nz) .and. (dabs(t1) .lt. qaltf0nz) .and. (
     *dabs(t2) .lt. qaltf0nz))then
      goto 23186
      endif
23185 continue
23186 continue
      if(0 .le. kpzavbj3)then
      d0(gp1jxzuh) = f0
      endif
      if(1 .le. kpzavbj3)then
      d1(gp1jxzuh) = f1
      endif
      if(2 .le. kpzavbj3)then
      d2(gp1jxzuh) = f2
      endif
23173 continue
23174 continue
      return
      end
