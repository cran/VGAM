      subroutine qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      implicit logical (a-z)
      integer wy1vqfzu, tgiyxdw1(1), dufozmt7(1)
      integer urohxe6t, bpvaqm5z, ayfnwr1v
      ayfnwr1v = 1
      urohxe6t = wy1vqfzu
23000 if(.not.(urohxe6t.ge.1))goto 23002
      do 23003 bpvaqm5z=1,urohxe6t 
      tgiyxdw1(ayfnwr1v) = bpvaqm5z
      ayfnwr1v = ayfnwr1v+1
23003 continue
       urohxe6t=urohxe6t-1
      goto 23000
23002 continue
      ayfnwr1v = 1
      do 23005 urohxe6t=1,wy1vqfzu 
      do 23007 bpvaqm5z=urohxe6t,wy1vqfzu 
      dufozmt7(ayfnwr1v) = bpvaqm5z
      ayfnwr1v = ayfnwr1v+1
23007 continue
23005 continue
      return
      end
      integer function viamf(cz8qdfyj, rvy1fpli, wy1vqfzu, tgiyxdw1, 
     &dufozmt7)
      integer cz8qdfyj, rvy1fpli, wy1vqfzu, tgiyxdw1(1), dufozmt7(1)
      integer urohxe6t, imk5wjxg
      imk5wjxg = wy1vqfzu*(wy1vqfzu+1)/2
      do 23009 urohxe6t=1,imk5wjxg 
      if(.not.((tgiyxdw1(urohxe6t).eq.cz8qdfyj .and. dufozmt7(urohxe6t)
     &.eq.rvy1fpli) .or.(tgiyxdw1(urohxe6t).eq.rvy1fpli .and. dufozmt7(
     &urohxe6t).eq.cz8qdfyj)))goto 23011
      viamf = urohxe6t
      return
23011 continue
23009 continue
      viamf = 0
      return
      end
      subroutine vm2af(mat, a, dimm, tgiyxdw1, dufozmt7, kuzxj1lo, 
     &wy1vqfzu, upper)
      implicit logical (a-z)
      integer dimm, tgiyxdw1(dimm), dufozmt7(dimm), kuzxj1lo, wy1vqfzu, 
     &upper
      double precision mat(dimm,kuzxj1lo), a(wy1vqfzu,wy1vqfzu,kuzxj1lo)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, imk5wjxg
      imk5wjxg = wy1vqfzu * (wy1vqfzu + 1) / 2
      if(.not.(upper .eq. 1 .or. dimm .ne. imk5wjxg))goto 23013
      ayfnwr1v = 1
23015 if(.not.(ayfnwr1v.le.kuzxj1lo))goto 23017
      yq6lorbx = 1
23018 if(.not.(yq6lorbx.le.wy1vqfzu))goto 23020
      gp1jxzuh = 1
23021 if(.not.(gp1jxzuh.le.wy1vqfzu))goto 23023
      a(gp1jxzuh,yq6lorbx,ayfnwr1v) = 0.0d0
       gp1jxzuh=gp1jxzuh+1
      goto 23021
23023 continue
       yq6lorbx=yq6lorbx+1
      goto 23018
23020 continue
       ayfnwr1v=ayfnwr1v+1
      goto 23015
23017 continue
23013 continue
      do 23024 ayfnwr1v=1,kuzxj1lo 
      do 23026 yq6lorbx=1,dimm 
      a(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx),ayfnwr1v) = mat(yq6lorbx,
     &ayfnwr1v)
      if(.not.(upper .eq. 0))goto 23028
      a(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx),ayfnwr1v) = mat(yq6lorbx,
     &ayfnwr1v)
23028 continue
23026 continue
23024 continue
      return
      end
      subroutine nudh6szqf(wpuarq2m, tlgduey8, lfu2qhid, dimu, tgiyxdw1,
     & dufozmt7, kuzxj1lo, wy1vqfzu, wk1200)
      implicit logical (a-z)
      integer dimu, tgiyxdw1(1), dufozmt7(1), kuzxj1lo, wy1vqfzu
      double precision wpuarq2m(dimu,kuzxj1lo), tlgduey8(kuzxj1lo,
     &wy1vqfzu), lfu2qhid(wy1vqfzu,kuzxj1lo), wk1200(wy1vqfzu,wy1vqfzu)
      double precision q6zdcwxk
      integer ayfnwr1v, yq6lorbx, bpvaqm5z, one, upper
      one = 1
      upper = 1
      ayfnwr1v = 1
23030 if(.not.(ayfnwr1v.le.kuzxj1lo))goto 23032
      call vm2af(wpuarq2m(1,ayfnwr1v), wk1200, dimu, tgiyxdw1, dufozmt7,
     & one, wy1vqfzu, upper)
      yq6lorbx = 1
23033 if(.not.(yq6lorbx.le.wy1vqfzu))goto 23035
      q6zdcwxk = 0.0d0
      bpvaqm5z = yq6lorbx
23036 if(.not.(bpvaqm5z.le.wy1vqfzu))goto 23038
      q6zdcwxk = q6zdcwxk + wk1200(yq6lorbx,bpvaqm5z) * tlgduey8(
     &ayfnwr1v,bpvaqm5z)
       bpvaqm5z=bpvaqm5z+1
      goto 23036
23038 continue
      lfu2qhid(yq6lorbx,ayfnwr1v) = q6zdcwxk
       yq6lorbx=yq6lorbx+1
      goto 23033
23035 continue
       ayfnwr1v=ayfnwr1v+1
      goto 23030
23032 continue
      return
      end
      subroutine vbksf(wpuarq2m, bvecto, wy1vqfzu, kuzxj1lo, wk1200, 
     &tgiyxdw1, dufozmt7, dimu)
      implicit logical (a-z)
      integer wy1vqfzu, kuzxj1lo, tgiyxdw1(1), dufozmt7(1), dimu
      double precision wpuarq2m(dimu,kuzxj1lo), bvecto(wy1vqfzu,
     &kuzxj1lo), wk1200(wy1vqfzu,wy1vqfzu)
      double precision q6zdcwxk
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, upper, one
      upper = 1
      one = 1
      ayfnwr1v = 1
23039 if(.not.(ayfnwr1v.le.kuzxj1lo))goto 23041
      call vm2af(wpuarq2m(1,ayfnwr1v), wk1200, dimu, tgiyxdw1, dufozmt7,
     & one, wy1vqfzu, upper)
      yq6lorbx = wy1vqfzu
23042 if(.not.(yq6lorbx.ge.1))goto 23044
      q6zdcwxk = bvecto(yq6lorbx,ayfnwr1v)
      gp1jxzuh = yq6lorbx+1
23045 if(.not.(gp1jxzuh.le.wy1vqfzu))goto 23047
      q6zdcwxk = q6zdcwxk - wk1200(yq6lorbx,gp1jxzuh) * bvecto(gp1jxzuh,
     &ayfnwr1v)
       gp1jxzuh=gp1jxzuh+1
      goto 23045
23047 continue
      bvecto(yq6lorbx,ayfnwr1v) = q6zdcwxk / wk1200(yq6lorbx,yq6lorbx)
       yq6lorbx=yq6lorbx-1
      goto 23042
23044 continue
       ayfnwr1v=ayfnwr1v+1
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
      do 23048 ayfnwr1v=1,wy1vqfzu
      q6zdcwxk = 0d0
      do 23050 gp1jxzuh=1,ayfnwr1v-1 
      q6zdcwxk = q6zdcwxk + wmat(gp1jxzuh,ayfnwr1v) * wmat(gp1jxzuh,
     &ayfnwr1v)
23050 continue
      wmat(ayfnwr1v,ayfnwr1v) = wmat(ayfnwr1v,ayfnwr1v) - q6zdcwxk
      if(.not.(wmat(ayfnwr1v,ayfnwr1v) .le. 0d0))goto 23052
      dvhw1ulq = 0
      return
23052 continue
      wmat(ayfnwr1v,ayfnwr1v) = dsqrt(wmat(ayfnwr1v,ayfnwr1v))
      do 23054 yq6lorbx=ayfnwr1v+1,wy1vqfzu
      q6zdcwxk = 0d0
      do 23056 gp1jxzuh=1,ayfnwr1v-1 
      q6zdcwxk = q6zdcwxk + wmat(gp1jxzuh,ayfnwr1v) * wmat(gp1jxzuh,
     &yq6lorbx)
23056 continue
      wmat(ayfnwr1v,yq6lorbx) = (wmat(ayfnwr1v,yq6lorbx) - q6zdcwxk) / 
     &wmat(ayfnwr1v,ayfnwr1v)
23054 continue
23048 continue
      if(.not.(isolve .eq. 0))goto 23058
      do 23060 ayfnwr1v=2,wy1vqfzu 
      do 23062 yq6lorbx=1,ayfnwr1v-1 
      wmat(ayfnwr1v,yq6lorbx) = 0.0d0
23062 continue
      return
23060 continue
23058 continue
      do 23064 yq6lorbx=1,wy1vqfzu 
      q6zdcwxk = bvecto(yq6lorbx)
      do 23066 gp1jxzuh=1,yq6lorbx-1 
      q6zdcwxk = q6zdcwxk - wmat(gp1jxzuh,yq6lorbx) * bvecto(gp1jxzuh)
23066 continue
      bvecto(yq6lorbx) = q6zdcwxk / wmat(yq6lorbx,yq6lorbx)
23064 continue
      yq6lorbx = wy1vqfzu
23068 if(.not.(yq6lorbx.ge.1))goto 23070
      q6zdcwxk = bvecto(yq6lorbx)
      gp1jxzuh = yq6lorbx+1
23071 if(.not.(gp1jxzuh.le.wy1vqfzu))goto 23073
      q6zdcwxk = q6zdcwxk - wmat(yq6lorbx,gp1jxzuh) * bvecto(gp1jxzuh)
       gp1jxzuh=gp1jxzuh+1
      goto 23071
23073 continue
      bvecto(yq6lorbx) = q6zdcwxk / wmat(yq6lorbx,yq6lorbx)
       yq6lorbx=yq6lorbx-1
      goto 23068
23070 continue
      return
      end
      subroutine mxrbkut0f(wpuarq2m, he7mqnvy, wy1vqfzu, xjc4ywlh, 
     &kuzxj1lo, wk1200, wk3400, tgiyxdw1, dufozmt7, dimu, rutyk8mg)
      implicit logical (a-z)
      integer dimu, wy1vqfzu, xjc4ywlh, kuzxj1lo, tgiyxdw1(1), dufozmt7(
     &1), rutyk8mg
      double precision wpuarq2m(dimu,kuzxj1lo), he7mqnvy(rutyk8mg,
     &xjc4ywlh), wk1200(wy1vqfzu,wy1vqfzu), wk3400(wy1vqfzu,xjc4ywlh)
      double precision q6zdcwxk
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, bpvaqm5z
      do 23074 yq6lorbx=1,wy1vqfzu 
      do 23076 ayfnwr1v=1,wy1vqfzu 
      wk1200(ayfnwr1v,yq6lorbx) = 0.0d0
23076 continue
23074 continue
      do 23078 ayfnwr1v=1,kuzxj1lo 
      do 23080 bpvaqm5z=1,dimu 
      wk1200(tgiyxdw1(bpvaqm5z), dufozmt7(bpvaqm5z)) = wpuarq2m(
     &bpvaqm5z,ayfnwr1v)
23080 continue
      do 23082 gp1jxzuh=1,xjc4ywlh 
      do 23084 yq6lorbx=1,wy1vqfzu 
      wk3400(yq6lorbx,gp1jxzuh) = he7mqnvy((ayfnwr1v-1)*wy1vqfzu+
     &yq6lorbx,gp1jxzuh)
23084 continue
23082 continue
      do 23086 gp1jxzuh=1,xjc4ywlh 
      do 23088 yq6lorbx=1,wy1vqfzu 
      q6zdcwxk = 0d0
      do 23090 bpvaqm5z=yq6lorbx,wy1vqfzu 
      q6zdcwxk = q6zdcwxk + wk1200(yq6lorbx,bpvaqm5z) * wk3400(bpvaqm5z,
     &gp1jxzuh)
23090 continue
      he7mqnvy((ayfnwr1v-1)*wy1vqfzu+yq6lorbx,gp1jxzuh) = q6zdcwxk
23088 continue
23086 continue
23078 continue
      return
      end
      subroutine vrinvf9(wpuarq2m, ldr, wy1vqfzu, dvhw1ulq, ks3wejcv, 
     &work)
      implicit logical (a-z)
      integer ldr, wy1vqfzu, dvhw1ulq
      double precision wpuarq2m(ldr,wy1vqfzu), ks3wejcv(wy1vqfzu,
     &wy1vqfzu), work(wy1vqfzu,wy1vqfzu)
      double precision q6zdcwxk
      integer yq6lorbx, gp1jxzuh, col, uaoynef0
      dvhw1ulq = 1
      yq6lorbx = 1
23092 if(.not.(yq6lorbx.le.wy1vqfzu))goto 23094
      col = 1
23095 if(.not.(col.le.wy1vqfzu))goto 23097
      work(yq6lorbx,col) = 0.0d0
       col=col+1
      goto 23095
23097 continue
       yq6lorbx=yq6lorbx+1
      goto 23092
23094 continue
      col = 1
23098 if(.not.(col.le.wy1vqfzu))goto 23100
      yq6lorbx = col
23101 if(.not.(yq6lorbx.ge.1))goto 23103
      if(.not.(yq6lorbx .eq. col))goto 23104
      q6zdcwxk = 1.0d0
      goto 23105
23104 continue
      q6zdcwxk = 0.0d0
23105 continue
      gp1jxzuh = yq6lorbx+1
23106 if(.not.(gp1jxzuh.le.col))goto 23108
      q6zdcwxk = q6zdcwxk - wpuarq2m(yq6lorbx,gp1jxzuh) * work(gp1jxzuh,
     &col)
       gp1jxzuh=gp1jxzuh+1
      goto 23106
23108 continue
      if(.not.(wpuarq2m(yq6lorbx,yq6lorbx) .eq. 0.0d0))goto 23109
      dvhw1ulq = 0
      goto 23110
23109 continue
      work(yq6lorbx,col) = q6zdcwxk / wpuarq2m(yq6lorbx,yq6lorbx)
23110 continue
       yq6lorbx=yq6lorbx-1
      goto 23101
23103 continue
       col=col+1
      goto 23098
23100 continue
      yq6lorbx = 1
23111 if(.not.(yq6lorbx.le.wy1vqfzu))goto 23113
      col = yq6lorbx
23114 if(.not.(col.le.wy1vqfzu))goto 23116
      if(.not.(yq6lorbx .lt. col))goto 23117
      uaoynef0 = col
      goto 23118
23117 continue
      uaoynef0 = yq6lorbx
23118 continue
      q6zdcwxk = 0.0d0
      gp1jxzuh = uaoynef0
23119 if(.not.(gp1jxzuh.le.wy1vqfzu))goto 23121
      q6zdcwxk = q6zdcwxk + work(yq6lorbx,gp1jxzuh) * work(col,gp1jxzuh)
       gp1jxzuh=gp1jxzuh+1
      goto 23119
23121 continue
      ks3wejcv(yq6lorbx,col) = q6zdcwxk
      ks3wejcv(col,yq6lorbx) = q6zdcwxk
       col=col+1
      goto 23114
23116 continue
       yq6lorbx=yq6lorbx+1
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
23122 if(.not.(yq6lorbx.le.6))goto 23124
      y = y + 1.0d0
      q6zdcwxk = q6zdcwxk + xd4mybgj(yq6lorbx)/y
       yq6lorbx=yq6lorbx+1
      goto 23122
23124 continue
      lfu2qhid = -hofjnx2e + dlog(2.5066282746310005d0 * q6zdcwxk / x)
      return
      end
      subroutine enbin9(bzmd6ftv, hdqsx7bk, nm0eljqk, n2kersmx, n, 
     &dvhw1ulq, zy1mchbf, ux3nadiw, rsynp1go, sguwj9ty)
      implicit logical (a-z)
      integer n, dvhw1ulq, zy1mchbf, sguwj9ty
      double precision bzmd6ftv(n, zy1mchbf), hdqsx7bk(n, zy1mchbf), 
     &nm0eljqk(n, zy1mchbf), n2kersmx, ux3nadiw, rsynp1go
      integer ayfnwr1v, kij0gwer
      double precision oxjgzv0e, btiehdm2, ydb, vjz5sxty, esql7umk, 
     &pvcjl2na, mwuvskg1, ft3ijqmy, hmayv1xt, q6zdcwxk, plo6hkdr
      real csi9ydge
      if(.not.(n2kersmx .le. 0.80d0 .or. n2kersmx .ge. 1.0d0))goto 23125
      dvhw1ulq = 0
      return
23125 continue
      btiehdm2 = 100.0d0 * rsynp1go
      oxjgzv0e = 0.001d0
      dvhw1ulq = 1
      kij0gwer=1
23127 if(.not.(kij0gwer.le.zy1mchbf))goto 23129
      ayfnwr1v=1
23130 if(.not.(ayfnwr1v.le.n))goto 23132
      vjz5sxty = nm0eljqk(ayfnwr1v,kij0gwer) / hdqsx7bk(ayfnwr1v,
     &kij0gwer)
      if(.not.((vjz5sxty .lt. oxjgzv0e) .or. (nm0eljqk(ayfnwr1v,
     &kij0gwer) .gt. 1.0d5)))goto 23133
      bzmd6ftv(ayfnwr1v,kij0gwer) = -nm0eljqk(ayfnwr1v,kij0gwer) * (1.
     &0d0 + hdqsx7bk(ayfnwr1v,kij0gwer)/(hdqsx7bk(ayfnwr1v,kij0gwer) + 
     &nm0eljqk(ayfnwr1v,kij0gwer))) / hdqsx7bk(ayfnwr1v,kij0gwer)**2
      if(.not.(bzmd6ftv(ayfnwr1v,kij0gwer) .gt. -btiehdm2))goto 23135
      bzmd6ftv(ayfnwr1v,kij0gwer) = -btiehdm2
23135 continue
      goto 20
23133 continue
      q6zdcwxk = 0.0d0
      pvcjl2na = hdqsx7bk(ayfnwr1v,kij0gwer) / (hdqsx7bk(ayfnwr1v,
     &kij0gwer) + nm0eljqk(ayfnwr1v,kij0gwer))
      mwuvskg1 = 1.0d0 - pvcjl2na
      csi9ydge = hdqsx7bk(ayfnwr1v,kij0gwer)
      if(.not.(pvcjl2na .lt. btiehdm2))goto 23137
      pvcjl2na = btiehdm2
23137 continue
      if(.not.(mwuvskg1 .lt. btiehdm2))goto 23139
      mwuvskg1 = btiehdm2
23139 continue
      esql7umk = 100.0d0 + 15.0d0 * nm0eljqk(ayfnwr1v,kij0gwer)
      if(.not.(esql7umk .lt. sguwj9ty))goto 23141
      esql7umk = sguwj9ty
23141 continue
      ft3ijqmy = pvcjl2na ** csi9ydge
      ux3nadiw = ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / hdqsx7bk(ayfnwr1v,kij0gwer)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = 1.0d0
      ft3ijqmy = hdqsx7bk(ayfnwr1v,kij0gwer) * mwuvskg1 * ft3ijqmy
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + 
     &ydb)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = 2.0d0
23143 if(.not.(((ux3nadiw .le. n2kersmx) .or. (plo6hkdr .gt. 1.0d-4)) .
     &and.(ydb .lt. esql7umk)))goto 23144
      ft3ijqmy = (hdqsx7bk(ayfnwr1v,kij0gwer) - 1.0d0 + ydb) * mwuvskg1 
     &* ft3ijqmy / ydb
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + 
     &ydb)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = ydb + 1.0d0
      goto 23143
23144 continue
      bzmd6ftv(ayfnwr1v,kij0gwer) = -q6zdcwxk
20    hmayv1xt = 0.0d0
       ayfnwr1v=ayfnwr1v+1
      goto 23130
23132 continue
       kij0gwer=kij0gwer+1
      goto 23127
23129 continue
      return
      end
      subroutine enbin8(bzmd6ftv, hdqsx7bk, hsj9bzaq, n2kersmx, 
     &kuzxj1lo, dvhw1ulq, zy1mchbf, ux3nadiw, rsynp1go)
      implicit logical (a-z)
      integer kuzxj1lo, dvhw1ulq, zy1mchbf
      double precision bzmd6ftv(kuzxj1lo, zy1mchbf), hdqsx7bk(kuzxj1lo, 
     &zy1mchbf), hsj9bzaq(kuzxj1lo, zy1mchbf), n2kersmx, ux3nadiw, 
     &rsynp1go
      integer ayfnwr1v, kij0gwer, esql7umk
      double precision ft3ijqmy, tad5vhsu, o3jyipdf, pq0hfucn, q6zdcwxk,
     & d1, d2, plo6hkdr, hnu1vjyw
      logical pok1, pok2, pok12
      double precision oxjgzv0e, onemse, nm0eljqk, btiehdm2, ydb, kbig
      d1 = 0.0d0
      d2 = 0.0d0
      btiehdm2 = -100.0d0 * rsynp1go
      esql7umk = 3000
      if(.not.(n2kersmx .le. 0.80d0 .or. n2kersmx .ge. 1.0d0))goto 23145
      dvhw1ulq = 0
      return
23145 continue
      kbig = 1.0d4
      oxjgzv0e = 0.001d0
      hnu1vjyw = 1.0d0 - rsynp1go
      onemse = 1.0d0 / (1.0d0 + oxjgzv0e)
      dvhw1ulq = 1
      kij0gwer=1
23147 if(.not.(kij0gwer.le.zy1mchbf))goto 23149
      ayfnwr1v=1
23150 if(.not.(ayfnwr1v.le.kuzxj1lo))goto 23152
      if(.not.(hdqsx7bk(ayfnwr1v,kij0gwer) .gt. kbig))goto 23153
      hdqsx7bk(ayfnwr1v,kij0gwer) = kbig
23153 continue
      if(.not.(hsj9bzaq(ayfnwr1v,kij0gwer) .lt. oxjgzv0e))goto 23155
      hsj9bzaq(ayfnwr1v,kij0gwer) = oxjgzv0e
23155 continue
      if(.not.((hsj9bzaq(ayfnwr1v,kij0gwer) .gt. onemse)))goto 23157
      nm0eljqk = hdqsx7bk(ayfnwr1v,kij0gwer) * (1.0d0/hsj9bzaq(ayfnwr1v,
     &kij0gwer) - 1.0d0)
      bzmd6ftv(ayfnwr1v,kij0gwer) = -nm0eljqk * (1.0d0 + hdqsx7bk(
     &ayfnwr1v,kij0gwer)/(hdqsx7bk(ayfnwr1v,kij0gwer) + nm0eljqk)) / 
     &hdqsx7bk(ayfnwr1v,kij0gwer)**2
      if(.not.(bzmd6ftv(ayfnwr1v,kij0gwer) .gt. btiehdm2))goto 23159
      bzmd6ftv(ayfnwr1v,kij0gwer) = btiehdm2
23159 continue
      goto 20
23157 continue
      q6zdcwxk = 0.0d0
      pok1 = .true.
      pok2 = hsj9bzaq(ayfnwr1v,kij0gwer) .lt. (1.0d0-rsynp1go)
      pok12 = pok1 .and. pok2
      if(.not.(pok12))goto 23161
      d2 = hdqsx7bk(ayfnwr1v,kij0gwer) * dlog(hsj9bzaq(ayfnwr1v,
     &kij0gwer))
      ux3nadiw = dexp(d2)
      goto 23162
23161 continue
      ux3nadiw = 0.0d0
23162 continue
      plo6hkdr = (1.0d0 - ux3nadiw) / hdqsx7bk(ayfnwr1v,kij0gwer)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      call tldz5ion(hdqsx7bk(ayfnwr1v,kij0gwer), o3jyipdf)
      ydb = 1.0d0
      call tldz5ion(ydb + hdqsx7bk(ayfnwr1v,kij0gwer), tad5vhsu)
      pq0hfucn = 0.0d0
      if(.not.(pok12))goto 23163
      d1 = dlog(1.0d0 - hsj9bzaq(ayfnwr1v,kij0gwer))
      ft3ijqmy = dexp(ydb * d1 + d2 + tad5vhsu - o3jyipdf - pq0hfucn)
      goto 23164
23163 continue
      ft3ijqmy = 0.0d0
23164 continue
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + 
     &ydb)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = 2.0d0
23165 if(.not.((ux3nadiw .le. n2kersmx) .or. (plo6hkdr .gt. 1.0d-4)))
     &goto 23166
      tad5vhsu = tad5vhsu + dlog(ydb + hdqsx7bk(ayfnwr1v,kij0gwer) - 1.
     &0d0)
      pq0hfucn = pq0hfucn + dlog(ydb)
      if(.not.(pok12))goto 23167
      ft3ijqmy = dexp(ydb * d1 + d2 + tad5vhsu - o3jyipdf - pq0hfucn)
      goto 23168
23167 continue
      ft3ijqmy = 0.0d0
23168 continue
      ux3nadiw = ux3nadiw + ft3ijqmy
      plo6hkdr = (1.0d0 - ux3nadiw) / (hdqsx7bk(ayfnwr1v,kij0gwer) + 
     &ydb)**2
      q6zdcwxk = q6zdcwxk + plo6hkdr
      ydb = ydb + 1.0d0
      if(.not.(ydb .gt. 1.0d3))goto 23169
      goto 21
23169 continue
      goto 23165
23166 continue
21    bzmd6ftv(ayfnwr1v,kij0gwer) = -q6zdcwxk
20    tad5vhsu = 0.0d0
       ayfnwr1v=ayfnwr1v+1
      goto 23150
23152 continue
       kij0gwer=kij0gwer+1
      goto 23147
23149 continue
      return
      end
      subroutine mbessi0(bvecto, kuzxj1lo, kpzavbj3, d0, d1, d2, 
     &zjkrtol8, qaltf0nz)
      implicit logical (a-z)
      integer kuzxj1lo, kpzavbj3, zjkrtol8, c5aesxkus
      double precision bvecto(kuzxj1lo), d0(kuzxj1lo), d1(kuzxj1lo), d2(
     &kuzxj1lo), qaltf0nz
      integer ayfnwr1v, gp1jxzuh
      double precision f0, t0, m0, f1, t1, m1, f2, t2, m2
      double precision toobig
      toobig = 20.0d0
      zjkrtol8 = 0
      if(.not.(.not.(kpzavbj3 .eq. 0 .or. kpzavbj3 .eq. 1 .or. kpzavbj3 
     &.eq. 2)))goto 23171
      zjkrtol8 = 1
      return
23171 continue
      do 23173 gp1jxzuh=1,kuzxj1lo 
      if(.not.(dabs(bvecto(gp1jxzuh)) .gt. toobig))goto 23175
      zjkrtol8 = 1
      return
23175 continue
      t1 = bvecto(gp1jxzuh) / 2.0d0
      f1 = t1
      t0 = t1 * t1
      f0 = 1.0d0 + t0
      t2 = 0.50d0
      f2 = t2
      c5aesxkus = 15
      if(.not.(dabs(bvecto(gp1jxzuh)) .gt. 10))goto 23177
      c5aesxkus = 25
23177 continue
      if(.not.(dabs(bvecto(gp1jxzuh)) .gt. 15))goto 23179
      c5aesxkus = 35
23179 continue
      if(.not.(dabs(bvecto(gp1jxzuh)) .gt. 20))goto 23181
      c5aesxkus = 40
23181 continue
      if(.not.(dabs(bvecto(gp1jxzuh)) .gt. 30))goto 23183
      c5aesxkus = 55
23183 continue
      do 23185 ayfnwr1v=1,c5aesxkus 
      m0 = (bvecto(gp1jxzuh) / (2.0d0*(ayfnwr1v+1.0d0))) ** 2.0
      m1 = m0 * (1.0d0 + 1.0d0/ayfnwr1v)
      m2 = m1 * (2.0d0*ayfnwr1v + 1.0d0) / (2.0d0*ayfnwr1v - 1.0d0)
      t0 = t0 * m0
      t1 = t1 * m1
      t2 = t2 * m2
      f0 = f0 + t0
      f1 = f1 + t1
      f2 = f2 + t2
      if(.not.((dabs(t0) .lt. qaltf0nz) .and. (dabs(t1) .lt. qaltf0nz) 
     &.and. (dabs(t2) .lt. qaltf0nz)))goto 23187
      goto 23186
23187 continue
23185 continue
23186 continue
      if(.not.(0 .le. kpzavbj3))goto 23189
      d0(gp1jxzuh) = f0
23189 continue
      if(.not.(1 .le. kpzavbj3))goto 23191
      d1(gp1jxzuh) = f1
23191 continue
      if(.not.(2 .le. kpzavbj3))goto 23193
      d2(gp1jxzuh) = f2
23193 continue
23173 continue
      return
      end
