      subroutine qh4ulb(zqve1l, vvl1li, lku8xq)
      implicit logical (a-z)
      integer lku8xq, zqve1l(1), vvl1li(1)
      integer xi1mqb, i1nkrb, w3gohz
      w3gohz = 1
      xi1mqb = lku8xq
23000 if(.not.(xi1mqb.ge.1))goto 23002
      do 23003 i1nkrb=1,xi1mqb 
      zqve1l(w3gohz) = i1nkrb
      w3gohz = w3gohz+1
23003 continue
       xi1mqb=xi1mqb-1
      goto 23000
23002 continue
      w3gohz = 1
      do 23005 xi1mqb=1,lku8xq 
      do 23007 i1nkrb=xi1mqb,lku8xq 
      vvl1li(w3gohz) = i1nkrb
      w3gohz = w3gohz+1
23007 continue
23005 continue
      return
      end
      integer function viamf(s17te9, xl6qgm, lku8xq, zqve1l, vvl1li)
      integer s17te9, xl6qgm, lku8xq, zqve1l(1), vvl1li(1)
      integer xi1mqb, j0qwtz
      j0qwtz = lku8xq*(lku8xq+1)/2
      do 23009 xi1mqb=1,j0qwtz 
      if(.not.((zqve1l(xi1mqb).eq.s17te9 .and. vvl1li(xi1mqb).eq.xl6qgm)
     & .or.(zqve1l(xi1mqb).eq.xl6qgm .and. vvl1li(xi1mqb).eq.s17te9)))
     &goto 23011
      viamf = xi1mqb
      return
23011 continue
23009 continue
      viamf = 0
      return
      end
      subroutine vm2af(mat, a, p1yjqz, zqve1l, vvl1li, nfiumb4, lku8xq, 
     &teola6)
      implicit logical (a-z)
      integer p1yjqz, zqve1l(p1yjqz), vvl1li(p1yjqz), nfiumb4, lku8xq, 
     &teola6
      double precision mat(p1yjqz,nfiumb4), a(lku8xq,lku8xq,nfiumb4)
      integer w3gohz, d9rjek, nd6mep, j0qwtz
      j0qwtz = lku8xq * (lku8xq + 1) / 2
      if(.not.(teola6 .eq. 1 .or. p1yjqz .ne. j0qwtz))goto 23013
      w3gohz = 1
23015 if(.not.(w3gohz.le.nfiumb4))goto 23017
      d9rjek = 1
23018 if(.not.(d9rjek.le.lku8xq))goto 23020
      nd6mep = 1
23021 if(.not.(nd6mep.le.lku8xq))goto 23023
      a(nd6mep,d9rjek,w3gohz) = 0.0d0
       nd6mep=nd6mep+1
      goto 23021
23023 continue
       d9rjek=d9rjek+1
      goto 23018
23020 continue
       w3gohz=w3gohz+1
      goto 23015
23017 continue
23013 continue
      do 23024 w3gohz=1,nfiumb4 
      do 23026 d9rjek=1,p1yjqz 
      a(zqve1l(d9rjek),vvl1li(d9rjek),w3gohz) = mat(d9rjek,w3gohz)
      if(.not.(teola6 .eq. 0))goto 23028
      a(vvl1li(d9rjek),zqve1l(d9rjek),w3gohz) = mat(d9rjek,w3gohz)
23028 continue
23026 continue
23024 continue
      return
      end
      subroutine mux22f(jrxg6l, jmwo0z, ghry8z, zkjqhi, zqve1l, vvl1li, 
     &nfiumb4, lku8xq, mbd8lk)
      implicit logical (a-z)
      integer zkjqhi, zqve1l(1), vvl1li(1), nfiumb4, lku8xq
      double precision jrxg6l(zkjqhi,nfiumb4), jmwo0z(nfiumb4,lku8xq), 
     &ghry8z(lku8xq,nfiumb4), mbd8lk(lku8xq,lku8xq)
      double precision qnk4zf
      integer w3gohz, d9rjek, i1nkrb, one, teola6
      one = 1
      teola6 = 1
      w3gohz = 1
23030 if(.not.(w3gohz.le.nfiumb4))goto 23032
      call vm2af(jrxg6l(1,w3gohz), mbd8lk, zkjqhi, zqve1l, vvl1li, one, 
     &lku8xq, teola6)
      d9rjek = 1
23033 if(.not.(d9rjek.le.lku8xq))goto 23035
      qnk4zf = 0.0d0
      i1nkrb = d9rjek
23036 if(.not.(i1nkrb.le.lku8xq))goto 23038
      qnk4zf = qnk4zf + mbd8lk(d9rjek,i1nkrb) * jmwo0z(w3gohz,i1nkrb)
       i1nkrb=i1nkrb+1
      goto 23036
23038 continue
      ghry8z(d9rjek,w3gohz) = qnk4zf
       d9rjek=d9rjek+1
      goto 23033
23035 continue
       w3gohz=w3gohz+1
      goto 23030
23032 continue
      return
      end
      subroutine vbksf(jrxg6l, yg1jzv, lku8xq, nfiumb4, mbd8lk, zqve1l, 
     &vvl1li, zkjqhi)
      implicit logical (a-z)
      integer lku8xq, nfiumb4, zqve1l, vvl1li, zkjqhi
      double precision jrxg6l(zkjqhi,nfiumb4), yg1jzv(lku8xq,nfiumb4), 
     &mbd8lk(lku8xq,lku8xq)
      double precision qnk4zf
      integer w3gohz, d9rjek, nd6mep, teola6, one
      teola6 = 1
      one = 1
      w3gohz = 1
23039 if(.not.(w3gohz.le.nfiumb4))goto 23041
      call vm2af(jrxg6l(1,w3gohz), mbd8lk, zkjqhi, zqve1l, vvl1li, one, 
     &lku8xq, teola6)
      d9rjek = lku8xq
23042 if(.not.(d9rjek.ge.1))goto 23044
      qnk4zf = yg1jzv(d9rjek,w3gohz)
      nd6mep = d9rjek+1
23045 if(.not.(nd6mep.le.lku8xq))goto 23047
      qnk4zf = qnk4zf - mbd8lk(d9rjek,nd6mep) * yg1jzv(nd6mep,w3gohz)
       nd6mep=nd6mep+1
      goto 23045
23047 continue
      yg1jzv(d9rjek,w3gohz) = qnk4zf / mbd8lk(d9rjek,d9rjek)
       d9rjek=d9rjek-1
      goto 23042
23044 continue
       w3gohz=w3gohz+1
      goto 23039
23041 continue
      return
      end
      subroutine vcholf(w8xfic, yg1jzv, lku8xq, c4uxow, ex7hfo)
      implicit logical (a-z)
      integer ex7hfo
      integer lku8xq, c4uxow
      double precision w8xfic(lku8xq,lku8xq), yg1jzv(lku8xq)
      double precision qnk4zf, dsqrt
      integer w3gohz, d9rjek, nd6mep
      c4uxow=1
      do 23048 w3gohz=1,lku8xq
      qnk4zf = 0d0
      do 23050 nd6mep=1,w3gohz-1 
      qnk4zf = qnk4zf + w8xfic(nd6mep,w3gohz) * w8xfic(nd6mep,w3gohz)
23050 continue
      w8xfic(w3gohz,w3gohz) = w8xfic(w3gohz,w3gohz) - qnk4zf
      if(.not.(w8xfic(w3gohz,w3gohz) .le. 0d0))goto 23052
      c4uxow = 0
      return
23052 continue
      w8xfic(w3gohz,w3gohz) = dsqrt(w8xfic(w3gohz,w3gohz))
      do 23054 d9rjek=w3gohz+1,lku8xq
      qnk4zf = 0d0
      do 23056 nd6mep=1,w3gohz-1 
      qnk4zf = qnk4zf + w8xfic(nd6mep,w3gohz) * w8xfic(nd6mep,d9rjek)
23056 continue
      w8xfic(w3gohz,d9rjek) = (w8xfic(w3gohz,d9rjek) - qnk4zf) / w8xfic(
     &w3gohz,w3gohz)
23054 continue
23048 continue
      if(.not.(ex7hfo .eq. 0))goto 23058
      do 23060 w3gohz=2,lku8xq 
      do 23062 d9rjek=1,w3gohz-1 
      w8xfic(w3gohz,d9rjek) = 0.0d0
23062 continue
      return
23060 continue
23058 continue
      do 23064 d9rjek=1,lku8xq 
      qnk4zf = yg1jzv(d9rjek)
      do 23066 nd6mep=1,d9rjek-1 
      qnk4zf = qnk4zf - w8xfic(nd6mep,d9rjek) * yg1jzv(nd6mep)
23066 continue
      yg1jzv(d9rjek) = qnk4zf / w8xfic(d9rjek,d9rjek)
23064 continue
      d9rjek = lku8xq
23068 if(.not.(d9rjek.ge.1))goto 23070
      qnk4zf = yg1jzv(d9rjek)
      nd6mep = d9rjek+1
23071 if(.not.(nd6mep.le.lku8xq))goto 23073
      qnk4zf = qnk4zf - w8xfic(d9rjek,nd6mep) * yg1jzv(nd6mep)
       nd6mep=nd6mep+1
      goto 23071
23073 continue
      yg1jzv(d9rjek) = qnk4zf / w8xfic(d9rjek,d9rjek)
       d9rjek=d9rjek-1
      goto 23068
23070 continue
      return
      end
      subroutine mux17f(jrxg6l, p3vlea, lku8xq, o9ljyn, nfiumb4, mbd8lk,
     & cm6nof, zqve1l, vvl1li, zkjqhi, c4bdmu)
      implicit logical (a-z)
      integer zkjqhi, lku8xq, o9ljyn, nfiumb4, zqve1l(1), vvl1li(1), 
     &c4bdmu
      double precision jrxg6l(zkjqhi,nfiumb4), p3vlea(c4bdmu,o9ljyn), 
     &mbd8lk(lku8xq,lku8xq), cm6nof(lku8xq,o9ljyn)
      double precision qnk4zf
      integer w3gohz, d9rjek, nd6mep, i1nkrb
      do 23074 d9rjek=1,lku8xq 
      do 23076 w3gohz=1,lku8xq 
      mbd8lk(w3gohz,d9rjek) = 0.0d0
23076 continue
23074 continue
      do 23078 w3gohz=1,nfiumb4 
      do 23080 i1nkrb=1,zkjqhi 
      mbd8lk(zqve1l(i1nkrb), vvl1li(i1nkrb)) = jrxg6l(i1nkrb,w3gohz)
23080 continue
      do 23082 nd6mep=1,o9ljyn 
      do 23084 d9rjek=1,lku8xq 
      cm6nof(d9rjek,nd6mep) = p3vlea((w3gohz-1)*lku8xq+d9rjek,nd6mep)
23084 continue
23082 continue
      do 23086 nd6mep=1,o9ljyn 
      do 23088 d9rjek=1,lku8xq 
      qnk4zf = 0d0
      do 23090 i1nkrb=d9rjek,lku8xq 
      qnk4zf = qnk4zf + mbd8lk(d9rjek,i1nkrb) * cm6nof(i1nkrb,nd6mep)
23090 continue
      p3vlea((w3gohz-1)*lku8xq+d9rjek,nd6mep) = qnk4zf
23088 continue
23086 continue
23078 continue
      return
      end
      subroutine vrinvf9(jrxg6l, ldr, lku8xq, c4uxow, ku0goz, bgu6fw)
      implicit logical (a-z)
      integer ldr, lku8xq, c4uxow
      double precision jrxg6l(ldr,lku8xq), ku0goz(lku8xq,lku8xq), 
     &bgu6fw(lku8xq,lku8xq)
      double precision qnk4zf
      integer d9rjek, nd6mep, col, mavy5hmod
      c4uxow = 1
      d9rjek = 1
23092 if(.not.(d9rjek.le.lku8xq))goto 23094
      col = 1
23095 if(.not.(col.le.lku8xq))goto 23097
      bgu6fw(d9rjek,col) = 0.0d0
       col=col+1
      goto 23095
23097 continue
       d9rjek=d9rjek+1
      goto 23092
23094 continue
      col = 1
23098 if(.not.(col.le.lku8xq))goto 23100
      d9rjek = col
23101 if(.not.(d9rjek.ge.1))goto 23103
      if(.not.(d9rjek .eq. col))goto 23104
      qnk4zf = 1.0d0
      goto 23105
23104 continue
      qnk4zf = 0.0d0
23105 continue
      nd6mep = d9rjek+1
23106 if(.not.(nd6mep.le.col))goto 23108
      qnk4zf = qnk4zf - jrxg6l(d9rjek,nd6mep) * bgu6fw(nd6mep,col)
       nd6mep=nd6mep+1
      goto 23106
23108 continue
      if(.not.(jrxg6l(d9rjek,d9rjek) .eq. 0.0d0))goto 23109
      c4uxow = 0
      goto 23110
23109 continue
      bgu6fw(d9rjek,col) = qnk4zf / jrxg6l(d9rjek,d9rjek)
23110 continue
       d9rjek=d9rjek-1
      goto 23101
23103 continue
       col=col+1
      goto 23098
23100 continue
      d9rjek = 1
23111 if(.not.(d9rjek.le.lku8xq))goto 23113
      col = d9rjek
23114 if(.not.(col.le.lku8xq))goto 23116
      if(.not.(d9rjek .lt. col))goto 23117
      mavy5hmod = col
      goto 23118
23117 continue
      mavy5hmod = d9rjek
23118 continue
      qnk4zf = 0.0d0
      nd6mep = mavy5hmod
23119 if(.not.(nd6mep.le.lku8xq))goto 23121
      qnk4zf = qnk4zf + bgu6fw(d9rjek,nd6mep) * bgu6fw(col,nd6mep)
       nd6mep=nd6mep+1
      goto 23119
23121 continue
      ku0goz(d9rjek,col) = qnk4zf
      ku0goz(col,d9rjek) = qnk4zf
       col=col+1
      goto 23114
23116 continue
       d9rjek=d9rjek+1
      goto 23111
23113 continue
      return
      end
      subroutine atez9d(xx, ghry8z)
      implicit logical (a-z)
      double precision xx, ghry8z
      double precision x, y, j0izmn, qnk4zf, mu4ygk(6)
      integer d9rjek
      mu4ygk(1)= 76.18009172947146d0
      mu4ygk(2)= -86.50532032941677d0
      mu4ygk(3)= 24.01409824083091d0
      mu4ygk(4)= -1.231739572450155d0
      mu4ygk(5)= 0.1208650973866179d-2
      mu4ygk(6)= -0.5395239384953d-5
      x = xx
      y = xx
      j0izmn = x+5.50d0
      j0izmn = j0izmn - (x+0.50d0) * dlog(j0izmn)
      qnk4zf=1.000000000190015d0
      d9rjek=1
23122 if(.not.(d9rjek.le.6))goto 23124
      y = y + 1.0d0
      qnk4zf = qnk4zf + mu4ygk(d9rjek)/y
       d9rjek=d9rjek+1
      goto 23122
23124 continue
      ghry8z = -j0izmn + dlog(2.5066282746310005d0 * qnk4zf / x)
      return
      end
      subroutine enbin9(qe3jcd, gq815b, xkcm3b, ogq67o, n, c4uxow, 
     &lzgs0f, hmr3dx, jftq1, nh2qxl)
      implicit logical (a-z)
      integer n, c4uxow, lzgs0f, nh2qxl
      double precision qe3jcd(n, lzgs0f), gq815b(n, lzgs0f), xkcm3b(n, 
     &lzgs0f), ogq67o, hmr3dx, jftq1
      integer w3gohz, myx3od
      double precision nh5zwa, pohw8d, ydb, dyb3po1, qbca1x, pvl5mc, 
     &pdjzm4, rxe0so, epx9jf, qnk4zf, scnrp6
      real yogfz6
      if(.not.(ogq67o .le. 0.80d0 .or. ogq67o .ge. 1.0d0))goto 23125
      c4uxow = 0
      return
23125 continue
      pohw8d = 100.0d0 * jftq1
      nh5zwa = 0.001d0
      c4uxow = 1
      myx3od=1
23127 if(.not.(myx3od.le.lzgs0f))goto 23129
      w3gohz=1
23130 if(.not.(w3gohz.le.n))goto 23132
      dyb3po1 = xkcm3b(w3gohz,myx3od) / gq815b(w3gohz,myx3od)
      if(.not.((dyb3po1 .lt. nh5zwa) .or. (xkcm3b(w3gohz,myx3od) .gt. 1.
     &0d5)))goto 23133
      qe3jcd(w3gohz,myx3od) = -xkcm3b(w3gohz,myx3od) * (1.0d0 + gq815b(
     &w3gohz,myx3od)/(gq815b(w3gohz,myx3od) + xkcm3b(w3gohz,myx3od))) / 
     &gq815b(w3gohz,myx3od)**2
      if(.not.(qe3jcd(w3gohz,myx3od) .gt. -pohw8d))goto 23135
      qe3jcd(w3gohz,myx3od) = -pohw8d
23135 continue
      goto 20
23133 continue
      qnk4zf = 0.0d0
      pvl5mc = gq815b(w3gohz,myx3od) / (gq815b(w3gohz,myx3od) + xkcm3b(
     &w3gohz,myx3od))
      pdjzm4 = 1.0d0 - pvl5mc
      yogfz6 = gq815b(w3gohz,myx3od)
      if(.not.(pvl5mc .lt. pohw8d))goto 23137
      pvl5mc = pohw8d
23137 continue
      if(.not.(pdjzm4 .lt. pohw8d))goto 23139
      pdjzm4 = pohw8d
23139 continue
      qbca1x = 100.0d0 + 15.0d0 * xkcm3b(w3gohz,myx3od)
      if(.not.(qbca1x .lt. nh2qxl))goto 23141
      qbca1x = nh2qxl
23141 continue
      rxe0so = pvl5mc ** yogfz6
      hmr3dx = rxe0so
      scnrp6 = (1.0d0 - hmr3dx) / gq815b(w3gohz,myx3od)**2
      qnk4zf = qnk4zf + scnrp6
      ydb = 1.0d0
      rxe0so = gq815b(w3gohz,myx3od) * pdjzm4 * rxe0so
      hmr3dx = hmr3dx + rxe0so
      scnrp6 = (1.0d0 - hmr3dx) / (gq815b(w3gohz,myx3od) + ydb)**2
      qnk4zf = qnk4zf + scnrp6
      ydb = 2.0d0
23143 if(.not.(((hmr3dx .le. ogq67o) .or. (scnrp6 .gt. 1.0d-4)) .and.(
     &ydb .lt. qbca1x)))goto 23144
      rxe0so = (gq815b(w3gohz,myx3od) - 1.0d0 + ydb) * pdjzm4 * rxe0so /
     & ydb
      hmr3dx = hmr3dx + rxe0so
      scnrp6 = (1.0d0 - hmr3dx) / (gq815b(w3gohz,myx3od) + ydb)**2
      qnk4zf = qnk4zf + scnrp6
      ydb = ydb + 1.0d0
      goto 23143
23144 continue
      qe3jcd(w3gohz,myx3od) = -qnk4zf
20    epx9jf = 0.0d0
       w3gohz=w3gohz+1
      goto 23130
23132 continue
       myx3od=myx3od+1
      goto 23127
23129 continue
      return
      end
      subroutine enbin8(qe3jcd, gq815b, ncrb2f, ogq67o, nfiumb4, c4uxow,
     & lzgs0f, hmr3dx, jftq1)
      implicit logical (a-z)
      integer nfiumb4, c4uxow, lzgs0f
      double precision qe3jcd(nfiumb4, lzgs0f), gq815b(nfiumb4, lzgs0f),
     & ncrb2f(nfiumb4, lzgs0f), ogq67o, hmr3dx, jftq1
      integer w3gohz, myx3od, qbca1x
      double precision rxe0so, mw6reg, xkwp2m, xndw5e, qnk4zf, d1, d2, 
     &scnrp6, onemeps
      logical pok1, pok2, pok12
      double precision nh5zwa, hntu8v, xkcm3b, pohw8d, ydb, kbig
      d1 = 0.0d0
      d2 = 0.0d0
      pohw8d = -100.0d0 * jftq1
      qbca1x = 3000
      if(.not.(ogq67o .le. 0.80d0 .or. ogq67o .ge. 1.0d0))goto 23145
      c4uxow = 0
      return
23145 continue
      kbig = 1.0d4
      nh5zwa = 0.001d0
      onemeps = 1.0d0 - jftq1
      hntu8v = 1.0d0 / (1.0d0 + nh5zwa)
      c4uxow = 1
      myx3od=1
23147 if(.not.(myx3od.le.lzgs0f))goto 23149
      w3gohz=1
23150 if(.not.(w3gohz.le.nfiumb4))goto 23152
      if(.not.(gq815b(w3gohz,myx3od) .gt. kbig))goto 23153
      gq815b(w3gohz,myx3od) = kbig
23153 continue
      if(.not.(ncrb2f(w3gohz,myx3od) .lt. nh5zwa))goto 23155
      ncrb2f(w3gohz,myx3od) = nh5zwa
23155 continue
      if(.not.((ncrb2f(w3gohz,myx3od) .gt. hntu8v)))goto 23157
      xkcm3b = gq815b(w3gohz,myx3od) * (1.0d0/ncrb2f(w3gohz,myx3od) - 1.
     &0d0)
      qe3jcd(w3gohz,myx3od) = -xkcm3b * (1.0d0 + gq815b(w3gohz,myx3od)/(
     &gq815b(w3gohz,myx3od) + xkcm3b)) / gq815b(w3gohz,myx3od)**2
      if(.not.(qe3jcd(w3gohz,myx3od) .gt. pohw8d))goto 23159
      qe3jcd(w3gohz,myx3od) = pohw8d
23159 continue
      goto 20
23157 continue
      qnk4zf = 0.0d0
      pok1 = .true.
      pok2 = ncrb2f(w3gohz,myx3od) .lt. (1.0d0-jftq1)
      pok12 = pok1 .and. pok2
      if(.not.(pok12))goto 23161
      d2 = gq815b(w3gohz,myx3od) * dlog(ncrb2f(w3gohz,myx3od))
      hmr3dx = dexp(d2)
      goto 23162
23161 continue
      hmr3dx = 0.0d0
23162 continue
      scnrp6 = (1.0d0 - hmr3dx) / gq815b(w3gohz,myx3od)**2
      qnk4zf = qnk4zf + scnrp6
      call atez9d(gq815b(w3gohz,myx3od), xkwp2m)
      ydb = 1.0d0
      call atez9d(ydb + gq815b(w3gohz,myx3od), mw6reg)
      xndw5e = 0.0d0
      if(.not.(pok12))goto 23163
      d1 = dlog(1.0d0 - ncrb2f(w3gohz,myx3od))
      rxe0so = dexp(ydb * d1 + d2 + mw6reg - xkwp2m - xndw5e)
      goto 23164
23163 continue
      rxe0so = 0.0d0
23164 continue
      hmr3dx = hmr3dx + rxe0so
      scnrp6 = (1.0d0 - hmr3dx) / (gq815b(w3gohz,myx3od) + ydb)**2
      qnk4zf = qnk4zf + scnrp6
      ydb = 2.0d0
23165 if(.not.((hmr3dx .le. ogq67o) .or. (scnrp6 .gt. 1.0d-4)))goto 2316
     &6
      mw6reg = mw6reg + dlog(ydb + gq815b(w3gohz,myx3od) - 1.0d0)
      xndw5e = xndw5e + dlog(ydb)
      if(.not.(pok12))goto 23167
      rxe0so = dexp(ydb * d1 + d2 + mw6reg - xkwp2m - xndw5e)
      goto 23168
23167 continue
      rxe0so = 0.0d0
23168 continue
      hmr3dx = hmr3dx + rxe0so
      scnrp6 = (1.0d0 - hmr3dx) / (gq815b(w3gohz,myx3od) + ydb)**2
      qnk4zf = qnk4zf + scnrp6
      ydb = ydb + 1.0d0
      if(.not.(ydb .gt. 1.0d3))goto 23169
      goto 21
23169 continue
      goto 23165
23166 continue
21    qe3jcd(w3gohz,myx3od) = -qnk4zf
20    mw6reg = 0.0d0
       w3gohz=w3gohz+1
      goto 23150
23152 continue
       myx3od=myx3od+1
      goto 23147
23149 continue
      return
      end
      subroutine mbessi0(yg1jzv, nfiumb4, xt3fko, d0, d1, d2, gqxvz8, 
     &kqoy6w)
      implicit logical (a-z)
      integer nfiumb4, xt3fko, gqxvz8, pga6nus
      double precision yg1jzv(nfiumb4), d0(nfiumb4), d1(nfiumb4), d2(
     &nfiumb4), kqoy6w
      integer w3gohz, nd6mep
      double precision f0, t0, m0, f1, t1, m1, f2, t2, m2
      double precision toobig
      toobig = 20.0d0
      if(.not.(.not.(xt3fko .eq. 0 .or. xt3fko .eq. 1 .or. xt3fko .eq. 
     &2)))goto 23171
      gqxvz8 = 1
      return
23171 continue
      gqxvz8 = 0
23172 continue
      do 23173 nd6mep=1,nfiumb4 
      if(.not.(dabs(yg1jzv(nd6mep)) .gt. toobig))goto 23175
      gqxvz8 = 1
      return
23175 continue
      t1 = yg1jzv(nd6mep) / 2.0d0
      f1 = t1
      t0 = t1 * t1
      f0 = 1.0d0 + t0
      t2 = 0.50d0
      f2 = t2
      pga6nus = 15
      if(.not.(dabs(yg1jzv(nd6mep)) .gt. 10))goto 23177
      pga6nus = 25
23177 continue
      if(.not.(dabs(yg1jzv(nd6mep)) .gt. 15))goto 23179
      pga6nus = 35
23179 continue
      if(.not.(dabs(yg1jzv(nd6mep)) .gt. 20))goto 23181
      pga6nus = 40
23181 continue
      if(.not.(dabs(yg1jzv(nd6mep)) .gt. 30))goto 23183
      pga6nus = 55
23183 continue
      do 23185 w3gohz=1,pga6nus 
      m0 = (yg1jzv(nd6mep) / (2.0d0*(w3gohz+1.0d0))) ** 2.0
      m1 = m0 * (1.0d0 + 1.0d0/w3gohz)
      m2 = m1 * (2.0d0*w3gohz + 1.0d0) / (2.0d0*w3gohz - 1.0d0)
      t0 = t0 * m0
      t1 = t1 * m1
      t2 = t2 * m2
      f0 = f0 + t0
      f1 = f1 + t1
      f2 = f2 + t2
      if(.not.((dabs(t0) .lt. kqoy6w) .and. (dabs(t1) .lt. kqoy6w) 
     &.and. (dabs(t2) .lt. kqoy6w)))goto 23187
      goto 23186
23187 continue
23185 continue
23186 continue
      if(.not.(0 .le. xt3fko))goto 23189
      d0(nd6mep) = f0
23189 continue
      if(.not.(1 .le. xt3fko))goto 23191
      d1(nd6mep) = f1
23191 continue
      if(.not.(2 .le. xt3fko))goto 23193
      d2(nd6mep) = f2
23193 continue
23173 continue
      return
      end
