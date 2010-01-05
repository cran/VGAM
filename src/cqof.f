      subroutine nw22ca(te4qac, ghry8z)
      implicit logical (a-z)
      double precision te4qac, ghry8z
      integer sn
      double precision r1, r2, y, y2, y3, y4, y5, y6, y7
      double precision erf, yhg7o7, z, z2, z3, z4
      double precision gln1k1, hqt8l8, a66epi, p10,p11,p12,p13, q10,q11,
     &q12,q13
      double precision p20,p21,p22,p23,p24,p25,p26,p27
      double precision q20,q21,q22,q23,q24,q25,q26,q27
      double precision p30,p31,p32,p33,p34
      double precision q30,q31,q32,q33,q34
      gln1k1 = 1.414213562373095049d0
      hqt8l8 = 1.772453850905516027d0
      a66epi = 20.0d0
      p10 = 242.66795523053175d0
      p11 = 21.979261618294152d0
      p12 = 6.9963834886191355d0
      p13 = -.035609843701815385d0
      q10 = 215.05887586986120d0
      q11 = 91.164905404514901d0
      q12 = 15.082797630407787d0
      q13 = 1.0d0
      p20 = 300.4592610201616005d0
      p21 = 451.9189537118729422d0
      p22 = 339.3208167343436870d0
      p23 = 152.9892850469404039d0
      p24 = 43.16222722205673530d0
      p25 = 7.211758250883093659d0
      p26 = .5641955174789739711d0
      p27 = -.0000001368648573827167067d0
      q20 = 300.4592609569832933d0
      q21 = 790.9509253278980272d0
      q22 = 931.3540948506096211d0
      q23 = 638.9802644656311665d0
      q24 = 277.5854447439876434d0
      q25 = 77.00015293522947295d0
      q26 = 12.78272731962942351d0
      q27 = 1.0d0
      p30 = -.00299610707703542174d0
      p31 = -.0494730910623250734d0
      p32 = -.226956593539686930d0
      p33 = -.278661308609647788d0
      p34 = -.0223192459734184686d0
      q30 = .0106209230528467918d0
      q31 = .191308926107829841d0
      q32 = 1.05167510706793207d0
      q33 = 1.98733201817135256d0
      q34 = 1.0d0
      if(.not.(te4qac .lt. -a66epi))goto 23000
      ghry8z = 2.753624d-89
      return
23000 continue
      if(.not.(te4qac .gt. a66epi))goto 23002
      ghry8z = 1.0d0
      return
23002 continue
      y = te4qac / gln1k1
      if(.not.(y .lt. 0.0d0))goto 23004
      y = -y
      sn = -1
      goto 23005
23004 continue
      sn = 1
23005 continue
      y2 = y * y
      y4 = y2 * y2
      y6 = y4 * y2
      if(.not.(y .lt. 0.46875d0))goto 23006
      r1 = p10 + p11 * y2 + p12 * y4 + p13 * y6
      r2 = q10 + q11 * y2 + q12 * y4 + q13 * y6
      erf = y * r1 / r2
      if(.not.(sn .eq. 1))goto 23008
      ghry8z = 0.5d0 + 0.5*erf
      goto 23009
23008 continue
      ghry8z = 0.5d0 - 0.5*erf
23009 continue
      goto 23007
23006 continue
      if(.not.(y .lt. 4.0d0))goto 23010
      y3 = y2 * y
      y5 = y4 * y
      y7 = y6 * y
      r1 = p20 + p21 * y + p22 * y2 + p23 * y3 + p24 * y4 + p25 * y5 + 
     &p26 * y6 + p27 * y7
      r2 = q20 + q21 * y + q22 * y2 + q23 * y3 + q24 * y4 + q25 * y5 + 
     &q26 * y6 + q27 * y7
      yhg7o7 = dexp(-y2) * r1 / r2
      if(.not.(sn .eq. 1))goto 23012
      ghry8z = 1.0 - 0.5*yhg7o7
      goto 23013
23012 continue
      ghry8z = 0.5*yhg7o7
23013 continue
      goto 23011
23010 continue
      z = y4
      z2 = z * z
      z3 = z2 * z
      z4 = z2 * z2
      r1 = p30 + p31 * z + p32 * z2 + p33 * z3 + p34 * z4
      r2 = q30 + q31 * z + q32 * z2 + q33 * z3 + q34 * z4
      yhg7o7 = (dexp(-y2)/y) * (1.0 / hqt8l8 + r1 / (r2 * y2))
      if(.not.(sn .eq. 1))goto 23014
      ghry8z = 1.0d0 - 0.5*yhg7o7
      goto 23015
23014 continue
      ghry8z = 0.5*yhg7o7
23015 continue
23011 continue
23007 continue
      return
      end
      subroutine pnm1ow(te4qac, ghry8z, nfiumb4)
      implicit logical (a-z)
      integer nfiumb4, w3gohz
      double precision te4qac(nfiumb4), ghry8z(nfiumb4)
      do 23016 w3gohz=1,nfiumb4 
      call nw22ca(te4qac(w3gohz), ghry8z(w3gohz))
23016 continue
      return
      end
      subroutine q4cgho(te4qac, dwgkz6, ghry8z)
      implicit logical (a-z)
      double precision te4qac, dwgkz6, ghry8z
      double precision mu4ygka
      if(.not.(1.0d0 - te4qac .ge. 1.0d0))goto 23018
      ghry8z = -8.12589d0 / (3.0*dsqrt(dwgkz6))
      goto 23019
23018 continue
      if(.not.(1.0d0 - te4qac .le. 0.0d0))goto 23020
      ghry8z = 8.12589d0 / (3.0*dsqrt(dwgkz6))
      goto 23021
23020 continue
      call nw22ca(1.0d0-te4qac, mu4ygka)
      mu4ygka = mu4ygka / (3.0*dsqrt(dwgkz6))
      ghry8z = -3.0d0 * dlog(1.0d0 + mu4ygka)
23021 continue
23019 continue
      return
      end
      subroutine wgf0al(te4qac, ghry8z)
      implicit logical (a-z)
      double precision te4qac, ghry8z
      if(.not.(1.0d0 - te4qac .ge. 1.0d0))goto 23022
      ghry8z = -35.0d0
      goto 23023
23022 continue
      if(.not.(1.0d0 - te4qac .le. 0.0d0))goto 23024
      ghry8z = 3.542106d0
      goto 23025
23024 continue
      ghry8z = dlog(-dlog(1.0d0 - te4qac))
23025 continue
23023 continue
      return
      end
      subroutine u10e3o(te4qac, ghry8z)
      implicit logical (a-z)
      double precision te4qac, ghry8z
      if(.not.(1.0d0 - te4qac .ge. 1.0d0))goto 23026
      ghry8z = -34.53958d0
      goto 23027
23026 continue
      if(.not.(1.0d0 - te4qac .le. 0.0d0))goto 23028
      ghry8z = 34.53958d0
      goto 23029
23028 continue
      ghry8z = dlog(te4qac / (1.0d0 - te4qac))
23029 continue
23027 continue
      return
      end
      subroutine pjw1l(ur73jo, lq8reh, go0l1q, nfiumb4, lku8xq, vi231l, 
     &zxiwf1, dyt0pg, lir0o1, zvxw1l, h3mrfq, q121lc)
      implicit logical (a-z)
      integer nfiumb4, lku8xq, vi231l, zxiwf1, dyt0pg, lir0o1, zvxw1l, 
     &h3mrfq
      double precision ur73jo(vi231l,zxiwf1), lq8reh(zxiwf1), go0l1q(
     &lku8xq,nfiumb4), q121lc(nfiumb4)
      integer w3gohz, d9rjek, nd6mep, opf6cv, c3qxjo
      double precision nqvu3e
      if(.not.(dyt0pg .eq. 1))goto 23030
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23032
      c3qxjo = 2*lir0o1-1
      do 23034 w3gohz=1,nfiumb4 
      nqvu3e = 0.0d0
      do 23036 nd6mep=1,zxiwf1 
      nqvu3e = nqvu3e + ur73jo(2*w3gohz-1,nd6mep) * lq8reh(nd6mep)
23036 continue
      go0l1q(c3qxjo,w3gohz) = nqvu3e
23034 continue
      c3qxjo = 2*lir0o1
      do 23038 w3gohz=1,nfiumb4 
      nqvu3e = 0.0d0
      do 23040 nd6mep=1,zxiwf1 
      nqvu3e = nqvu3e + ur73jo(2*w3gohz ,nd6mep) * lq8reh(nd6mep)
23040 continue
      go0l1q(c3qxjo,w3gohz) = nqvu3e
23038 continue
      goto 23033
23032 continue
      do 23042 w3gohz=1,vi231l 
      nqvu3e = 0.0d0
      do 23044 nd6mep=1,zxiwf1 
      nqvu3e = nqvu3e + ur73jo(w3gohz,nd6mep) * lq8reh(nd6mep)
23044 continue
      go0l1q(lir0o1,w3gohz) = nqvu3e
23042 continue
23033 continue
      goto 23031
23030 continue
      opf6cv = 1
      do 23046 w3gohz=1,nfiumb4 
      do 23048 d9rjek=1,lku8xq 
      nqvu3e = 0.0d0
      do 23050 nd6mep=1,zxiwf1 
      nqvu3e = nqvu3e + ur73jo(opf6cv,nd6mep) * lq8reh(nd6mep)
23050 continue
      opf6cv = opf6cv + 1
      go0l1q(d9rjek,w3gohz) = nqvu3e
23048 continue
23046 continue
23031 continue
      if(.not.(h3mrfq .eq. 1))goto 23052
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23054
      do 23056 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = go0l1q(2*lir0o1-1,w3gohz) + q121lc(
     &w3gohz)
23056 continue
      goto 23055
23054 continue
      do 23058 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = go0l1q(lir0o1,w3gohz) + q121lc(w3gohz)
23058 continue
23055 continue
23052 continue
      return
      end
      subroutine o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l,
     & lir0o1)
      implicit logical (a-z)
      integer nfiumb4, lku8xq, aqk377, zvxw1l, lir0o1
      double precision go0l1q(lku8xq,nfiumb4), w5poyv(aqk377,nfiumb4)
      integer w3gohz, d9rjek
      double precision xkwp2m0
      if(.not.(lir0o1 .eq. 0))goto 23060
      if(.not.(zvxw1l .eq. 1))goto 23062
      do 23064 w3gohz=1,nfiumb4 
      do 23066 d9rjek=1,lku8xq 
      xkwp2m0 = dexp(go0l1q(d9rjek,w3gohz))
      w5poyv(d9rjek,w3gohz) = xkwp2m0 / (1.0d0 + xkwp2m0)
23066 continue
23064 continue
23062 continue
      if(.not.(zvxw1l .eq. 2))goto 23068
      do 23070 w3gohz=1,nfiumb4 
      do 23072 d9rjek=1,lku8xq 
      w5poyv(d9rjek,w3gohz) = dexp(go0l1q(d9rjek,w3gohz))
23072 continue
23070 continue
23068 continue
      if(.not.(zvxw1l .eq. 4))goto 23074
      do 23076 w3gohz=1,nfiumb4 
      do 23078 d9rjek=1,lku8xq 
      w5poyv(d9rjek,w3gohz) = 1.0d0-dexp(-dexp(go0l1q(d9rjek,w3gohz)))
23078 continue
23076 continue
23074 continue
      if(.not.(zvxw1l .eq. 5))goto 23080
      do 23082 w3gohz=1,nfiumb4 
      do 23084 d9rjek=1,aqk377 
      w5poyv(d9rjek,w3gohz) = dexp(go0l1q(2*d9rjek-1,w3gohz))
23084 continue
23082 continue
23080 continue
      if(.not.(zvxw1l .eq. 3))goto 23086
      do 23088 w3gohz=1,nfiumb4 
      do 23090 d9rjek=1,aqk377 
      w5poyv(d9rjek,w3gohz) = dexp(go0l1q(2*d9rjek-1,w3gohz))
23090 continue
23088 continue
23086 continue
      if(.not.(zvxw1l .eq. 8))goto 23092
      do 23094 w3gohz=1,nfiumb4 
      do 23096 d9rjek=1,lku8xq 
      w5poyv(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz)
23096 continue
23094 continue
23092 continue
      goto 23061
23060 continue
      if(.not.(zvxw1l .eq. 1))goto 23098
      do 23100 w3gohz=1,nfiumb4 
      xkwp2m0 = dexp(go0l1q(lir0o1,w3gohz))
      w5poyv(lir0o1,w3gohz) = xkwp2m0 / (1.0d0 + xkwp2m0)
23100 continue
23098 continue
      if(.not.(zvxw1l .eq. 2))goto 23102
      do 23104 w3gohz=1,nfiumb4 
      w5poyv(lir0o1,w3gohz) = dexp(go0l1q(lir0o1,w3gohz))
23104 continue
23102 continue
      if(.not.(zvxw1l .eq. 4))goto 23106
      do 23108 w3gohz=1,nfiumb4 
      w5poyv(lir0o1,w3gohz) = 1.0d0 - dexp(-dexp(go0l1q(lir0o1,w3gohz)))
23108 continue
23106 continue
      if(.not.(zvxw1l .eq. 5))goto 23110
      do 23112 w3gohz=1,nfiumb4 
      w5poyv(lir0o1,w3gohz) = dexp(go0l1q(2*lir0o1-1,w3gohz))
23112 continue
23110 continue
      if(.not.(zvxw1l .eq. 3))goto 23114
      do 23116 w3gohz=1,nfiumb4 
      w5poyv(lir0o1,w3gohz) = dexp(go0l1q(2*lir0o1-1,w3gohz))
23116 continue
23114 continue
      if(.not.(zvxw1l .eq. 8))goto 23118
      do 23120 w3gohz=1,nfiumb4 
      w5poyv(lir0o1,w3gohz) = go0l1q(lir0o1,w3gohz)
23120 continue
23118 continue
23061 continue
      return
      end
      subroutine kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq,
     & aqk377, xhe4cg, go0l1q, dev, lir0o1, fiumb4, wlkaa3, cll)
      implicit logical (a-z)
      integer zvxw1l, nfiumb4, lku8xq, aqk377, xhe4cg, lir0o1, cll
      double precision jmwo0z(nfiumb4, aqk377), w8xfic(nfiumb4, xhe4cg),
     & w5poyv(aqk377, nfiumb4), go0l1q(lku8xq,nfiumb4), dev, fiumb4, 
     &wlkaa3
      integer w3gohz, d9rjek
      double precision qe3jcd, ue1phr, mu4ygk, ig5cma, j0izmn, smu, 
     &mqs3rp, mdk7tp, oesul2
      double precision gq815b, aho01l, dh3rio
      logical xyiu19
      qe3jcd = 0.0d0
      if(.not.(lir0o1 .eq. 0))goto 23122
      if(.not.((zvxw1l .eq. 1) .or. (zvxw1l .eq. 4)))goto 23124
      do 23126 d9rjek=1,lku8xq 
      do 23128 w3gohz=1,nfiumb4 
      if(.not.(jmwo0z(w3gohz,d9rjek) .gt. 0.0d0))goto 23130
      mdk7tp = jmwo0z(w3gohz,d9rjek) * dlog(jmwo0z(w3gohz,d9rjek))
      goto 23131
23130 continue
      mdk7tp = 0.0d0
23131 continue
      if(.not.(jmwo0z(w3gohz,d9rjek) .lt. 1.0d0))goto 23132
      mdk7tp = mdk7tp + (1.0d0 - jmwo0z(w3gohz,d9rjek)) * dlog(1.0d0 - 
     &jmwo0z(w3gohz,d9rjek))
23132 continue
      mu4ygk = w5poyv(d9rjek,w3gohz) * (1.0d0 - w5poyv(d9rjek,w3gohz))
      if(.not.(mu4ygk .lt. fiumb4))goto 23134
      smu = w5poyv(d9rjek,w3gohz)
      if(.not.(smu .lt. fiumb4))goto 23136
      oesul2 = jmwo0z(w3gohz,d9rjek) * wlkaa3
      goto 23137
23136 continue
      oesul2 = jmwo0z(w3gohz,d9rjek) * dlog(smu)
23137 continue
      mqs3rp = 1.0d0 - smu
      if(.not.(mqs3rp .lt. fiumb4))goto 23138
      oesul2 = oesul2 + (1.0d0 - jmwo0z(w3gohz,d9rjek)) * wlkaa3
      goto 23139
23138 continue
      oesul2 = oesul2 + (1.0d0 - jmwo0z(w3gohz,d9rjek)) * dlog(mqs3rp)
23139 continue
      goto 23135
23134 continue
      oesul2 = (jmwo0z(w3gohz,d9rjek) * dlog(w5poyv(d9rjek,w3gohz)) + (
     &1.0d0 - jmwo0z(w3gohz,d9rjek)) * dlog(1.0d0 - w5poyv(d9rjek,
     &w3gohz)))
23135 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * (mdk7tp - oesul2)
23128 continue
23126 continue
23124 continue
      if(.not.(zvxw1l .eq. 2))goto 23140
      do 23142 d9rjek=1,lku8xq 
      do 23144 w3gohz=1,nfiumb4 
      if(.not.(jmwo0z(w3gohz,d9rjek) .gt. 0.0d0))goto 23146
      mu4ygk = w5poyv(d9rjek,w3gohz) - jmwo0z(w3gohz,d9rjek) + jmwo0z(
     &w3gohz,d9rjek) * dlog(jmwo0z(w3gohz,d9rjek) / w5poyv(d9rjek,
     &w3gohz))
      goto 23147
23146 continue
      mu4ygk = w5poyv(d9rjek,w3gohz) - jmwo0z(w3gohz,d9rjek)
23147 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk
23144 continue
23142 continue
23140 continue
      if(.not.(zvxw1l .eq. 5))goto 23148
      do 23150 d9rjek=1,aqk377 
      do 23152 w3gohz=1,nfiumb4 
      dh3rio = dexp(go0l1q(2*d9rjek,w3gohz))
      call atez9d(dh3rio, ig5cma)
      if(.not.(jmwo0z(w3gohz,d9rjek) .gt. 0.0d0))goto 23154
      mu4ygk = (dh3rio - 1.0d0) * dlog(jmwo0z(w3gohz,d9rjek)) + (dlog(
     &dh3rio)-jmwo0z(w3gohz,d9rjek) / w5poyv(d9rjek,w3gohz) - dlog(
     &w5poyv(d9rjek,w3gohz)) ) * dh3rio - ig5cma
      goto 23155
23154 continue
      mu4ygk = -1000.0d0
23155 continue
      mu4ygk = -mu4ygk
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk
23152 continue
23150 continue
23148 continue
      if(.not.(zvxw1l .eq. 3))goto 23156
      if(.not.(cll .eq. 0))goto 23158
      aho01l = 34.0d0
      do 23160 d9rjek=1,aqk377 
      do 23162 w3gohz=1,nfiumb4 
      if(.not.(go0l1q(2*d9rjek,w3gohz) .gt. aho01l))goto 23164
      gq815b = dexp(aho01l)
      xyiu19 = .true.
      goto 23165
23164 continue
      if(.not.(go0l1q(2*d9rjek,w3gohz) .lt. -aho01l))goto 23166
      gq815b = dexp(-aho01l)
      xyiu19 = .true.
      goto 23167
23166 continue
      gq815b = dexp(go0l1q(2*d9rjek,w3gohz))
      xyiu19 = .false.
23167 continue
23165 continue
      if(.not.(jmwo0z(w3gohz,d9rjek) .lt. 1.0d0))goto 23168
      mu4ygk = 1.0d0
      goto 23169
23168 continue
      mu4ygk = jmwo0z(w3gohz,d9rjek)
23169 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * (jmwo0z(w3gohz,d9rjek) * 
     &dlog(mu4ygk/w5poyv(d9rjek,w3gohz)) + (jmwo0z(w3gohz,d9rjek) + 
     &gq815b) * dlog((w5poyv(d9rjek,w3gohz)+gq815b) / (gq815b+ jmwo0z(
     &w3gohz,d9rjek))))
23162 continue
23160 continue
      goto 23159
23158 continue
      aho01l = 34.0d0
      do 23170 d9rjek=1,aqk377 
      do 23172 w3gohz=1,nfiumb4 
      if(.not.(go0l1q(2*d9rjek,w3gohz) .gt. aho01l))goto 23174
      gq815b = dexp(aho01l)
      xyiu19 = .true.
      goto 23175
23174 continue
      if(.not.(go0l1q(2*d9rjek,w3gohz) .lt. -aho01l))goto 23176
      gq815b = dexp(-aho01l)
      xyiu19 = .true.
      goto 23177
23176 continue
      gq815b = dexp(go0l1q(2*d9rjek,w3gohz))
      xyiu19 = .false.
23177 continue
23175 continue
      if(.not.( xyiu19 ))goto 23178
      ig5cma = 0.0d0
      j0izmn = 0.0d0
      goto 23179
23178 continue
      call atez9d(gq815b + jmwo0z(w3gohz,d9rjek), ig5cma)
      call atez9d(gq815b, j0izmn)
23179 continue
      call atez9d(1.0d0 + jmwo0z(w3gohz,d9rjek), ue1phr)
      mu4ygk = gq815b * dlog(gq815b / (gq815b + w5poyv(d9rjek,w3gohz))) 
     &+ ig5cma - j0izmn - ue1phr
      if(.not.(jmwo0z(w3gohz,d9rjek) .gt. 0.0d0))goto 23180
      mu4ygk = mu4ygk + jmwo0z(w3gohz,d9rjek) * dlog(w5poyv(d9rjek,
     &w3gohz) / (gq815b + w5poyv(d9rjek,w3gohz)))
23180 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk
23172 continue
23170 continue
      qe3jcd = -qe3jcd / 2.0d0
23159 continue
23156 continue
      if(.not.(zvxw1l .eq. 8))goto 23182
      do 23184 d9rjek=1,lku8xq 
      do 23186 w3gohz=1,nfiumb4 
      mu4ygk = jmwo0z(w3gohz,d9rjek) - w5poyv(d9rjek,w3gohz)
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk**2
23186 continue
23184 continue
23182 continue
      goto 23123
23122 continue
      if(.not.((zvxw1l .eq. 1) .or. (zvxw1l .eq. 4)))goto 23188
      do 23190 w3gohz=1,nfiumb4 
      if(.not.(jmwo0z(w3gohz,lir0o1) .gt. 0.0d0))goto 23192
      mdk7tp = jmwo0z(w3gohz,lir0o1) * dlog(jmwo0z(w3gohz,lir0o1))
      goto 23193
23192 continue
      mdk7tp = 0.0d0
23193 continue
      if(.not.(jmwo0z(w3gohz,lir0o1) .lt. 1.0d0))goto 23194
      mdk7tp = mdk7tp + (1.0d0 - jmwo0z(w3gohz,lir0o1)) * dlog(1.0d0 - 
     &jmwo0z(w3gohz,lir0o1))
23194 continue
      mu4ygk = w5poyv(lir0o1,w3gohz) * (1.0d0 - w5poyv(lir0o1,w3gohz))
      if(.not.(mu4ygk .lt. fiumb4))goto 23196
      smu = w5poyv(lir0o1,w3gohz)
      if(.not.(smu .lt. fiumb4))goto 23198
      oesul2 = jmwo0z(w3gohz,lir0o1) * wlkaa3
      goto 23199
23198 continue
      oesul2 = jmwo0z(w3gohz,lir0o1) * dlog(smu)
23199 continue
      mqs3rp = 1.0d0 - smu
      if(.not.(mqs3rp .lt. fiumb4))goto 23200
      oesul2 = oesul2 + (1.0d0-jmwo0z(w3gohz,lir0o1))*wlkaa3
      goto 23201
23200 continue
      oesul2 = oesul2 + (1.0d0-jmwo0z(w3gohz,lir0o1))*dlog(mqs3rp)
23201 continue
      goto 23197
23196 continue
      oesul2 = (jmwo0z(w3gohz,lir0o1) * dlog(w5poyv(lir0o1,w3gohz)) + (
     &1.0d0 - jmwo0z(w3gohz,lir0o1)) * dlog(1.0d0 - w5poyv(lir0o1,
     &w3gohz)))
23197 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * (mdk7tp - oesul2)
23190 continue
23188 continue
      if(.not.(zvxw1l .eq. 2))goto 23202
      do 23204 w3gohz=1,nfiumb4 
      if(.not.(jmwo0z(w3gohz,lir0o1) .gt. 0.0d0))goto 23206
      mu4ygk = w5poyv(lir0o1,w3gohz) - jmwo0z(w3gohz,lir0o1) + jmwo0z(
     &w3gohz,lir0o1) * dlog(jmwo0z(w3gohz,lir0o1) / w5poyv(lir0o1,
     &w3gohz))
      goto 23207
23206 continue
      mu4ygk = w5poyv(lir0o1,w3gohz) - jmwo0z(w3gohz,lir0o1)
23207 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk
23204 continue
23202 continue
      if(.not.(zvxw1l .eq. 5))goto 23208
      do 23210 w3gohz=1,nfiumb4 
      dh3rio = dexp(go0l1q(2*lir0o1,w3gohz))
      call atez9d(dh3rio, ig5cma)
      if(.not.(jmwo0z(w3gohz,lir0o1) .gt. 0.0d0))goto 23212
      mu4ygk = (dh3rio - 1.0d0) * dlog(jmwo0z(w3gohz,lir0o1)) + dh3rio *
     & (dlog(dh3rio) - jmwo0z(w3gohz,lir0o1) / w5poyv(lir0o1,w3gohz) - 
     &dlog(w5poyv(lir0o1,w3gohz))) - ig5cma
      goto 23213
23212 continue
      mu4ygk = -1000.0d0
23213 continue
      mu4ygk = -mu4ygk
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk
23210 continue
23208 continue
      if(.not.(zvxw1l .eq. 3))goto 23214
      if(.not.(cll .eq. 0))goto 23216
      aho01l = 34.0d0
      do 23218 w3gohz=1,nfiumb4 
      if(.not.(go0l1q(2*lir0o1,w3gohz) .gt. aho01l))goto 23220
      gq815b = dexp(aho01l)
      xyiu19 = .true.
      goto 23221
23220 continue
      if(.not.(go0l1q(2*lir0o1,w3gohz) .lt. -aho01l))goto 23222
      gq815b = dexp(-aho01l)
      xyiu19 = .true.
      goto 23223
23222 continue
      gq815b = dexp(go0l1q(2*lir0o1,w3gohz))
      xyiu19 = .false.
23223 continue
23221 continue
      if(.not.(jmwo0z(w3gohz,lir0o1) .lt. 1.0d0))goto 23224
      mu4ygk = 1.0d0
      goto 23225
23224 continue
      mu4ygk = jmwo0z(w3gohz,lir0o1)
23225 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * (jmwo0z(w3gohz,lir0o1) * 
     &dlog(mu4ygk/w5poyv(lir0o1,w3gohz)) + (jmwo0z(w3gohz,lir0o1)+
     &gq815b) * dlog((w5poyv(lir0o1,w3gohz) + gq815b) / ( gq815b+jmwo0z(
     &w3gohz,lir0o1))))
23218 continue
      goto 23217
23216 continue
      do 23226 w3gohz=1,nfiumb4 
      gq815b = dexp(go0l1q(2*lir0o1,w3gohz))
      call atez9d(gq815b + jmwo0z(w3gohz,lir0o1), ig5cma)
      call atez9d(gq815b, j0izmn)
      call atez9d(1.0d0 + jmwo0z(w3gohz,lir0o1), ue1phr)
      mu4ygk = gq815b * dlog(gq815b / (gq815b + w5poyv(lir0o1,w3gohz))) 
     &+ ig5cma - j0izmn - ue1phr
      if(.not.(jmwo0z(w3gohz,lir0o1) .gt. 0.0d0))goto 23228
      mu4ygk = mu4ygk + jmwo0z(w3gohz,lir0o1) * dlog(w5poyv(lir0o1,
     &w3gohz) / (gq815b + w5poyv(lir0o1,w3gohz)))
23228 continue
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk
23226 continue
      qe3jcd = -qe3jcd / 2.0d0
23217 continue
23214 continue
      if(.not.(zvxw1l .eq. 8))goto 23230
      do 23232 w3gohz=1,nfiumb4 
      mu4ygk = jmwo0z(w3gohz,lir0o1) - w5poyv(lir0o1,w3gohz)
      qe3jcd = qe3jcd + w8xfic(w3gohz,1) * mu4ygk**2
23232 continue
23230 continue
23123 continue
      dev = 2.0d0 * qe3jcd
      return
      end
      subroutine sptoq8(hft28, ur73jo, nfiumb4, vi231l, cqui1v, zvxw1l)
      implicit logical (a-z)
      integer nfiumb4, vi231l, cqui1v, zvxw1l
      double precision hft28(nfiumb4,cqui1v), ur73jo(vi231l,1)
      integer w3gohz, c3qxjo, pvnfr4
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq.5 )))goto 23234
      c3qxjo = 1
      do 23236 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = 1.0d0
      ur73jo(2*w3gohz, c3qxjo) = 0.0d0
23236 continue
      c3qxjo = c3qxjo + 1
      do 23238 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = 0.0d0
      ur73jo(2*w3gohz, c3qxjo) = 1.0d0
23238 continue
      c3qxjo = c3qxjo + 1
      do 23240 pvnfr4=1,cqui1v 
      do 23242 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = hft28(w3gohz,pvnfr4)
      ur73jo(2*w3gohz, c3qxjo) = 0.0d0
23242 continue
      c3qxjo = c3qxjo + 1
23240 continue
      goto 23235
23234 continue
      c3qxjo = 1
      do 23244 w3gohz=1,nfiumb4 
      ur73jo(w3gohz,c3qxjo) = 1.0d0
23244 continue
      c3qxjo = c3qxjo + 1
      do 23246 pvnfr4=1,cqui1v 
      do 23248 w3gohz=1,nfiumb4 
      ur73jo(w3gohz,c3qxjo)=hft28(w3gohz,pvnfr4)
23248 continue
      c3qxjo = c3qxjo + 1
23246 continue
23235 continue
      return
      end
      subroutine u16zxj(hft28, ur73jo, nfiumb4, cqui1v, zvxw1l, q121lc, 
     &vi231l, zxiwf1, i5uvkm, zqve1l, vvl1li, oju3yh, p1, h3mrfq)
      implicit logical (a-z)
      integer nfiumb4, cqui1v, zvxw1l, vi231l, zxiwf1, i5uvkm, zqve1l(
     &i5uvkm), vvl1li(i5uvkm), p1, h3mrfq
      double precision hft28(nfiumb4,cqui1v), ur73jo(vi231l,zxiwf1), 
     &oju3yh(nfiumb4,p1)
      double precision q121lc(nfiumb4)
      integer hv3wja, w3gohz, c3qxjo, pvnfr4
      double precision mw6reg, ig5cma
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23250
      do 23252 pvnfr4=1,cqui1v 
      do 23254 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,pvnfr4) = hft28(w3gohz,pvnfr4)
      ur73jo(2*w3gohz ,pvnfr4) = 0.0d0
23254 continue
23252 continue
      c3qxjo = cqui1v + 1
      if(.not.(h3mrfq .eq. 0))goto 23256
      do 23258 hv3wja=1,i5uvkm 
      do 23260 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = hft28(w3gohz,zqve1l(hv3wja)) * hft28(
     &w3gohz,vvl1li(hv3wja))
      ur73jo(2*w3gohz ,c3qxjo) = 0.0d0
23260 continue
      c3qxjo = c3qxjo + 1
23258 continue
      goto 23257
23256 continue
      do 23262 w3gohz=1,nfiumb4 
      mw6reg = 0.0d0
      do 23264 pvnfr4=1,cqui1v 
      ig5cma = hft28(w3gohz,pvnfr4)
      mw6reg = mw6reg + ig5cma * ig5cma
23264 continue
      q121lc(w3gohz) = -0.50d0 * mw6reg
23262 continue
23257 continue
      goto 23251
23250 continue
      do 23266 pvnfr4=1,cqui1v 
      do 23268 w3gohz=1,nfiumb4 
      ur73jo(w3gohz,pvnfr4) = hft28(w3gohz,pvnfr4)
23268 continue
23266 continue
      c3qxjo = cqui1v + 1
      if(.not.(h3mrfq .eq. 0))goto 23270
      do 23272 hv3wja=1,i5uvkm 
      do 23274 w3gohz=1,nfiumb4 
      ur73jo(w3gohz,c3qxjo) = hft28(w3gohz,zqve1l(hv3wja)) * hft28(
     &w3gohz,vvl1li(hv3wja))
23274 continue
      c3qxjo = c3qxjo + 1
23272 continue
      goto 23271
23270 continue
      do 23276 w3gohz=1,nfiumb4 
      mw6reg = 0.0d0
      do 23278 pvnfr4=1,cqui1v 
      ig5cma = hft28(w3gohz,pvnfr4)
      mw6reg = mw6reg + ig5cma * ig5cma
23278 continue
      q121lc(w3gohz) = -0.50d0 * mw6reg
23276 continue
23271 continue
23251 continue
      if(.not.(p1 .gt. 0))goto 23280
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23282
      do 23284 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = 1.0d0
      ur73jo(2*w3gohz, c3qxjo) = 0.0d0
23284 continue
      c3qxjo = c3qxjo + 1
      do 23286 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = 0.0d0
      ur73jo(2*w3gohz, c3qxjo) = 1.0d0
23286 continue
      c3qxjo = c3qxjo + 1
      if(.not.(p1 .gt. 1))goto 23288
      do 23290 hv3wja=2,p1 
      do 23292 w3gohz=1,nfiumb4 
      ur73jo(2*w3gohz-1,c3qxjo) = oju3yh(w3gohz,hv3wja)
      ur73jo(2*w3gohz, c3qxjo) = 0.0d0
23292 continue
      c3qxjo = c3qxjo + 1
23290 continue
23288 continue
      goto 23283
23282 continue
      do 23294 hv3wja=1,p1 
      do 23296 w3gohz=1,nfiumb4 
      ur73jo(w3gohz,c3qxjo) = oju3yh(w3gohz,hv3wja)
23296 continue
      c3qxjo = c3qxjo + 1
23294 continue
23283 continue
23280 continue
      return
      end
      subroutine p0lk40(hft28, ur73jo, nfiumb4, lku8xq, vi231l, cqui1v, 
     &zvxw1l, aqk377, w5tcfp, cr8hav, i5uvkm, zqve1l, vvl1li, h3mrfq, 
     &q121lc)
      implicit logical (a-z)
      integer nfiumb4, lku8xq, vi231l, cqui1v, zvxw1l, aqk377, w5tcfp, 
     &cr8hav, i5uvkm, zqve1l(i5uvkm), vvl1li(i5uvkm), h3mrfq
      double precision hft28(nfiumb4,cqui1v), ur73jo(vi231l,cr8hav), 
     &q121lc(nfiumb4)
      integer hv3wja, w3gohz, d9rjek, nd6mep, ptr, c3qxjo, pvnfr4
      double precision ig5cma, mw6reg
      do 23298 nd6mep=1,cr8hav 
      do 23300 w3gohz=1,vi231l 
      ur73jo(w3gohz,nd6mep) = 0.0d0
23300 continue
23298 continue
      c3qxjo = 0
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23302
      do 23304 pvnfr4=1,cqui1v 
      ptr = 1
      do 23306 w3gohz=1,nfiumb4 
      do 23308 d9rjek=1,aqk377 
      ur73jo(ptr,c3qxjo+d9rjek) = hft28(w3gohz,pvnfr4)
      ptr = ptr + 2
23308 continue
23306 continue
      c3qxjo = c3qxjo + aqk377
23304 continue
      goto 23303
23302 continue
      do 23310 pvnfr4=1,cqui1v 
      ptr = 0
      do 23312 w3gohz=1,nfiumb4 
      do 23314 d9rjek=1,lku8xq 
      ptr = ptr + 1
      ur73jo(ptr,c3qxjo+d9rjek) = hft28(w3gohz,pvnfr4)
23314 continue
23312 continue
      c3qxjo = c3qxjo + lku8xq
23310 continue
23303 continue
      if(.not.(w5tcfp .eq. 0))goto 23316
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23318
      do 23320 hv3wja=1,i5uvkm 
      ptr = 1
      do 23322 w3gohz=1,nfiumb4 
      ig5cma = hft28(w3gohz,zqve1l(hv3wja)) * hft28(w3gohz,vvl1li(
     &hv3wja))
      do 23324 d9rjek=1,aqk377 
      ur73jo(ptr,c3qxjo+d9rjek) = ig5cma
      ptr = ptr + 2
23324 continue
23322 continue
      c3qxjo = c3qxjo + aqk377
23320 continue
      goto 23319
23318 continue
      do 23326 hv3wja=1,i5uvkm 
      ptr = 0
      do 23328 w3gohz=1,nfiumb4 
      ig5cma = hft28(w3gohz,zqve1l(hv3wja)) * hft28(w3gohz,vvl1li(
     &hv3wja))
      do 23330 d9rjek=1,lku8xq 
      ptr = ptr + 1
      ur73jo(ptr,c3qxjo+d9rjek) = ig5cma
23330 continue
23328 continue
      c3qxjo = c3qxjo + lku8xq
23326 continue
23319 continue
      goto 23317
23316 continue
      if(.not.(h3mrfq .eq. 1))goto 23332
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23334
      do 23336 w3gohz=1,nfiumb4 
      mw6reg = 0.0d0
      do 23338 pvnfr4=1,cqui1v 
      ig5cma = hft28(w3gohz,pvnfr4)
      mw6reg = mw6reg + ig5cma * ig5cma
23338 continue
      q121lc(w3gohz) = -0.50d0 * mw6reg
23336 continue
      goto 23335
23334 continue
      do 23340 w3gohz=1,nfiumb4 
      mw6reg = 0.0d0
      do 23342 pvnfr4=1,cqui1v 
      ig5cma = hft28(w3gohz,pvnfr4)
      mw6reg = mw6reg + ig5cma * ig5cma
23342 continue
      q121lc(w3gohz) = -0.50d0 * mw6reg
23340 continue
23335 continue
      goto 23333
23332 continue
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23344
      do 23346 hv3wja=1,i5uvkm 
      ptr = 1
      do 23348 w3gohz=1,nfiumb4 
      ig5cma = hft28(w3gohz,zqve1l(hv3wja)) * hft28(w3gohz,vvl1li(
     &hv3wja))
      do 23350 d9rjek=1,aqk377 
      ur73jo(ptr,c3qxjo+hv3wja) = ig5cma
      ptr = ptr + 2
23350 continue
23348 continue
23346 continue
      c3qxjo = c3qxjo + i5uvkm
      goto 23345
23344 continue
      do 23352 hv3wja=1,i5uvkm 
      ptr = 0
      do 23354 w3gohz=1,nfiumb4 
      ig5cma = hft28(w3gohz,zqve1l(hv3wja)) * hft28(w3gohz,vvl1li(
     &hv3wja))
      do 23356 d9rjek=1,lku8xq 
      ptr = ptr + 1
      ur73jo(ptr,c3qxjo+hv3wja) = ig5cma
23356 continue
23354 continue
23352 continue
      c3qxjo = c3qxjo + i5uvkm
23345 continue
23333 continue
23317 continue
      return
      end
      subroutine nbq4ua(jmwo0z, go0l1q, l1zvxx, nfiumb4, lku8xq, aqk377,
     & zvxw1l, lir0o1, w8xfic, foej1u)
      implicit logical (a-z)
      integer nfiumb4, lku8xq, aqk377, zvxw1l, lir0o1, foej1u
      double precision jmwo0z(nfiumb4,aqk377), go0l1q(lku8xq,nfiumb4), 
     &l1zvxx(15)
      double precision w8xfic(nfiumb4,1)
      double precision nqvu3e, pg2aqx, ozpqa0, ghys4c, qg8fdc, bdgzx3, 
     &cy0nqs, reg6st, wo8cqk
      integer w3gohz
      if(.not.((zvxw1l .eq. 1) .or. (zvxw1l .eq. 4) .or.(zvxw1l .eq. 3) 
     &.or. (zvxw1l .eq. 5)))goto 23358
      nqvu3e = 0.0d0
      pg2aqx = 0.0d0
      do 23360 w3gohz=1,nfiumb4 
      nqvu3e = nqvu3e + jmwo0z(w3gohz,lir0o1) * w8xfic(w3gohz,1)
      pg2aqx = pg2aqx + w8xfic(w3gohz,1)
23360 continue
      ozpqa0 = nqvu3e / pg2aqx
23358 continue
      if(.not.(zvxw1l .eq. 1))goto 23362
      call u10e3o(ozpqa0, ghys4c)
      do 23364 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = ghys4c
23364 continue
23362 continue
      if(.not.(zvxw1l .eq. 2))goto 23366
      do 23368 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = dlog(jmwo0z(w3gohz,lir0o1) + 0.125d0)
23368 continue
23366 continue
      if(.not.(zvxw1l .eq. 4))goto 23370
      call wgf0al(ozpqa0, qg8fdc)
      do 23372 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = qg8fdc
23372 continue
23370 continue
      if(.not.(zvxw1l .eq. 5))goto 23374
      if(.not.(foej1u .eq. 1))goto 23376
      bdgzx3 = dlog(ozpqa0 + 0.03125d0)
      cy0nqs = dlog(l1zvxx(3+aqk377+lir0o1)+0.01d0)
      do 23378 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = bdgzx3
      go0l1q(2*lir0o1, w3gohz) = cy0nqs
23378 continue
      goto 23377
23376 continue
      if(.not.(foej1u .eq. 2))goto 23380
      bdgzx3 = dlog((6.0/8.0)*ozpqa0+0.000d0)
      cy0nqs = dlog(l1zvxx(3+aqk377+lir0o1)+0.01d0)
      do 23382 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = bdgzx3
      go0l1q(2*lir0o1 ,w3gohz) = cy0nqs
23382 continue
      goto 23381
23380 continue
      cy0nqs = dlog(l1zvxx(3+aqk377+lir0o1)+0.01d0)
      do 23384 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = dlog(jmwo0z(w3gohz,lir0o1) + 0.
     &03125d0)
      go0l1q(2*lir0o1, w3gohz) = cy0nqs
23384 continue
23381 continue
23377 continue
23374 continue
      if(.not.(zvxw1l .eq. 3))goto 23386
      if(.not.(foej1u .eq. 1))goto 23388
      bdgzx3 = dlog(ozpqa0 + 0.03125d0)
      cy0nqs = dlog(l1zvxx(3+lir0o1)+0.03125d0)
      do 23390 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = bdgzx3
      go0l1q(2*lir0o1,w3gohz) = cy0nqs
23390 continue
      goto 23389
23388 continue
      if(.not.(foej1u .eq. 2))goto 23392
      bdgzx3 = dlog(ozpqa0 + 0.03125d0)
      wo8cqk = l1zvxx(3+lir0o1)
      cy0nqs = dlog(wo8cqk)
      do 23394 w3gohz=1,nfiumb4 
      reg6st = jmwo0z(w3gohz,lir0o1) - ozpqa0
      if(.not.(reg6st .gt. 3.0 * ozpqa0))goto 23396
      go0l1q(2*lir0o1-1,w3gohz) = dlog(dsqrt(jmwo0z(w3gohz,lir0o1)))
      go0l1q(2*lir0o1 ,w3gohz) = cy0nqs
      goto 23397
23396 continue
      go0l1q(2*lir0o1-1,w3gohz) = bdgzx3
      go0l1q(2*lir0o1 ,w3gohz) = cy0nqs
23397 continue
23394 continue
      goto 23393
23392 continue
      if(.not.(foej1u .eq. 3))goto 23398
      bdgzx3 = dlog(ozpqa0 + 0.03125d0)
      wo8cqk = l1zvxx(3+lir0o1)
      cy0nqs = dlog(wo8cqk)
      do 23400 w3gohz=1,nfiumb4 
      reg6st = jmwo0z(w3gohz,lir0o1) - ozpqa0
      if(.not.(reg6st .gt. ozpqa0))goto 23402
      go0l1q(2*lir0o1-1,w3gohz) = dlog(0.5*(jmwo0z(w3gohz,lir0o1)+
     &ozpqa0))
      go0l1q(2*lir0o1 ,w3gohz) = dlog(wo8cqk / (reg6st / ozpqa0))
      goto 23403
23402 continue
      if(.not.(jmwo0z(w3gohz,lir0o1) .lt. (ozpqa0 / 4.0)))goto 23404
      go0l1q(2*lir0o1-1,w3gohz) = dlog(ozpqa0 / 4.0)
      go0l1q(2*lir0o1 ,w3gohz) = cy0nqs
      goto 23405
23404 continue
      go0l1q(2*lir0o1-1,w3gohz) = bdgzx3
      go0l1q(2*lir0o1 ,w3gohz) = cy0nqs
23405 continue
23403 continue
23400 continue
      goto 23399
23398 continue
      cy0nqs = dlog(l1zvxx(3+lir0o1))
      do 23406 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = dlog(jmwo0z(w3gohz,lir0o1) + 0.
     &03125d0)
      go0l1q(2*lir0o1, w3gohz) = cy0nqs
23406 continue
23399 continue
23393 continue
23389 continue
23386 continue
      if(.not.(zvxw1l .eq. 8))goto 23408
      do 23410 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = jmwo0z(w3gohz,lir0o1)
23410 continue
23408 continue
      return
      end
      subroutine kqsxz1(jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, lj4dph, 
     &jrxg6l, jftq1, fiumb4, zl11l0, nfiumb4, lku8xq, aqk377, vi231l, 
     &zkjqhi, lir0o1, zvxw1l, gqxvz8, h3mrfq, q121lc)
      implicit logical (a-z)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, lir0o1, gqxvz8, 
     &h3mrfq
      double precision jmwo0z(nfiumb4,aqk377), w8xfic(nfiumb4,1), 
     &go0l1q(lku8xq,nfiumb4), w5poyv(aqk377,nfiumb4), q121lc(nfiumb4), 
     &hr83e(nfiumb4,lku8xq), lj4dph(nfiumb4,lku8xq), jrxg6l(zkjqhi,
     &nfiumb4), jftq1, fiumb4, zl11l0
      integer w3gohz, zvxw1l
      double precision mu4ygka, zixm0o, mu4ygkc, aho01l
      logical xyiu19
      double precision gq815b, p4je8, da51l0o, hmr3dx, yqsco4, ogq67o, 
     &qpzx6l(1,1), qxvi5(1,1), xkcm3b(1,1)
      integer uxzze7, c4uxow, nh2qxl
      double precision dh3rio, ig5cmad, ig5v8gzsp, dldshape
      double precision l0zqm, q1znur
      integer lqhm2g
      vi231l = 1
      uxzze7 = 1
      ogq67o = 0.990d0
      ogq67o = 0.995d0
      if(.not.(zvxw1l .eq. 1))goto 23412
      do 23414 w3gohz=1,nfiumb4 
      mu4ygka = w5poyv(lir0o1,w3gohz) * (1.0d0 - w5poyv(lir0o1,w3gohz))
      zixm0o = mu4ygka * w8xfic(w3gohz,1)
      if(.not.(mu4ygka .lt. fiumb4))goto 23416
      mu4ygka = fiumb4
23416 continue
      if(.not.(zixm0o .lt. fiumb4))goto 23418
      zixm0o = fiumb4
      jrxg6l(lir0o1,w3gohz) = zl11l0
      goto 23419
23418 continue
      jrxg6l(lir0o1,w3gohz) = dsqrt(zixm0o)
23419 continue
      lj4dph(w3gohz,lir0o1) = zixm0o
      hr83e(w3gohz,lir0o1) = go0l1q(lir0o1,w3gohz) + (jmwo0z(w3gohz,
     &lir0o1)-w5poyv(lir0o1,w3gohz)) / mu4ygka
23414 continue
23412 continue
      if(.not.(zvxw1l .eq. 2))goto 23420
      do 23422 w3gohz=1,nfiumb4 
      mu4ygka = w5poyv(lir0o1,w3gohz)
      zixm0o = mu4ygka * w8xfic(w3gohz,1)
      if(.not.(zixm0o .lt. fiumb4))goto 23424
      zixm0o = fiumb4
      jrxg6l(lir0o1,w3gohz) = zl11l0
      goto 23425
23424 continue
      jrxg6l(lir0o1,w3gohz) = dsqrt(zixm0o)
23425 continue
      lj4dph(w3gohz,lir0o1) = zixm0o
      if(.not.(jmwo0z(w3gohz,lir0o1) .gt. 0.0d0))goto 23426
      mu4ygkc = mu4ygka
      if(.not.(mu4ygkc .lt. fiumb4))goto 23428
      mu4ygkc = fiumb4
23428 continue
      hr83e(w3gohz,lir0o1) = go0l1q(lir0o1,w3gohz) + (jmwo0z(w3gohz,
     &lir0o1)-mu4ygkc)/mu4ygkc
      goto 23427
23426 continue
      hr83e(w3gohz,lir0o1) = go0l1q(lir0o1,w3gohz) - 1.0d0
23427 continue
23422 continue
23420 continue
      if(.not.(zvxw1l .eq. 4))goto 23430
      do 23432 w3gohz=1,nfiumb4 
      if(.not.((w5poyv(lir0o1,w3gohz) .lt. fiumb4) .or.(w5poyv(lir0o1,
     &w3gohz) .gt. 1.0d0 - fiumb4)))goto 23434
      mu4ygka = fiumb4
      zixm0o = mu4ygka * w8xfic(w3gohz,1)
      if(.not.(zixm0o .lt. fiumb4))goto 23436
      zixm0o = fiumb4
      jrxg6l(lir0o1,w3gohz) = zl11l0
      goto 23437
23436 continue
      jrxg6l(lir0o1,w3gohz) = dsqrt(zixm0o)
23437 continue
      lj4dph(w3gohz,lir0o1) = zixm0o
      hr83e(w3gohz,lir0o1) = go0l1q(lir0o1,w3gohz) + (jmwo0z(w3gohz,
     &lir0o1)-w5poyv(lir0o1,w3gohz)) / mu4ygka
      goto 23435
23434 continue
      mu4ygka = -(1.0d0 - w5poyv(lir0o1,w3gohz)) * dlog(1.0d0 - w5poyv(
     &lir0o1,w3gohz))
      if(.not.(mu4ygka .lt. fiumb4))goto 23438
      mu4ygka = fiumb4
23438 continue
      zixm0o = -mu4ygka * w8xfic(w3gohz,1) * dlog(1.0d0 - w5poyv(lir0o1,
     &w3gohz)) / w5poyv(lir0o1,w3gohz)
      if(.not.(zixm0o .lt. fiumb4))goto 23440
      zixm0o = fiumb4
23440 continue
      lj4dph(w3gohz,lir0o1) = zixm0o
      jrxg6l(lir0o1,w3gohz) = dsqrt(zixm0o)
      hr83e(w3gohz,lir0o1) = go0l1q(lir0o1,w3gohz) + (jmwo0z(w3gohz,
     &lir0o1)-w5poyv(lir0o1,w3gohz)) / mu4ygka
23435 continue
23432 continue
23430 continue
      if(.not.(zvxw1l .eq. 5))goto 23442
      l0zqm = 1.0d-20
      aho01l = 34.0d0
      do 23444 w3gohz=1,nfiumb4 
      if(.not.(go0l1q(2*lir0o1,w3gohz) .gt. aho01l))goto 23446
      dh3rio = dexp(aho01l)
      xyiu19 = .true.
      goto 23447
23446 continue
      if(.not.(go0l1q(2*lir0o1,w3gohz) .lt. -aho01l))goto 23448
      dh3rio = dexp(-aho01l)
      xyiu19 = .true.
      goto 23449
23448 continue
      dh3rio = dexp(go0l1q(2*lir0o1,w3gohz))
      xyiu19 = .false.
23449 continue
23447 continue
      call vdgam1(dh3rio, ig5cmad, lqhm2g)
      if(.not.(lqhm2g .ne. 1))goto 23450
      call intpr("error in kqsxz1 lqhm2g 1: ",-1,lqhm2g,1)
23450 continue
      q1znur = w5poyv(lir0o1,w3gohz)
      if(.not.(q1znur .lt. l0zqm))goto 23452
      q1znur = l0zqm
23452 continue
      dldshape = dlog(jmwo0z(w3gohz,lir0o1)) + dlog(dh3rio) - dlog(
     &q1znur) + 1.0d0 - ig5cmad - jmwo0z(w3gohz,lir0o1) / q1znur
      call vtgam1(dh3rio, ig5v8gzsp, lqhm2g)
      if(.not.(lqhm2g .ne. 1))goto 23454
      call intpr("error in kqsxz1 lqhm2g 2: ",-1,lqhm2g,1)
23454 continue
      lj4dph(w3gohz,2*lir0o1-1) = w8xfic(w3gohz,1) * dh3rio
      mu4ygka = dh3rio * ig5v8gzsp - 1.0d0
      lj4dph(w3gohz,2*lir0o1 ) = w8xfic(w3gohz,1) * dh3rio * mu4ygka
      if(.not.(lj4dph(w3gohz,2*lir0o1-1) .lt. fiumb4))goto 23456
      lj4dph(w3gohz,2*lir0o1-1) = fiumb4
      jrxg6l(2*lir0o1-1,w3gohz) = zl11l0
      goto 23457
23456 continue
      jrxg6l(2*lir0o1-1,w3gohz) = dsqrt(lj4dph(w3gohz,2*lir0o1-1))
23457 continue
      if(.not.(lj4dph(w3gohz,2*lir0o1) .lt. fiumb4))goto 23458
      lj4dph(w3gohz,2*lir0o1) = fiumb4
      jrxg6l(2*lir0o1,w3gohz) = zl11l0
      goto 23459
23458 continue
      jrxg6l(2*lir0o1,w3gohz) = dsqrt(lj4dph(w3gohz,2*lir0o1))
23459 continue
      if(.not.(mu4ygka .lt. l0zqm))goto 23460
      mu4ygka = l0zqm
23460 continue
      hr83e(w3gohz,2*lir0o1-1) = go0l1q(2*lir0o1-1,w3gohz) + jmwo0z(
     &w3gohz,lir0o1) / q1znur - 1.0d0
      hr83e(w3gohz,2*lir0o1 ) = go0l1q(2*lir0o1 ,w3gohz) + dldshape / 
     &mu4ygka
23444 continue
23442 continue
      if(.not.(zvxw1l .eq. 3))goto 23462
      aho01l = 34.0d0
      l0zqm = 1.0d-20
      do 23464 w3gohz=1,nfiumb4 
      if(.not.(go0l1q(2*lir0o1,w3gohz) .gt. aho01l))goto 23466
      gq815b = dexp(aho01l)
      xyiu19 = .true.
      goto 23467
23466 continue
      if(.not.(go0l1q(2*lir0o1,w3gohz) .lt. -aho01l))goto 23468
      gq815b = dexp(-aho01l)
      xyiu19 = .true.
      goto 23469
23468 continue
      gq815b = dexp(go0l1q(2*lir0o1,w3gohz))
      xyiu19 = .false.
23469 continue
23467 continue
      q1znur = w5poyv(lir0o1,w3gohz)
      if(.not.(q1znur .lt. l0zqm))goto 23470
      q1znur = l0zqm
23470 continue
      call vdgam1(jmwo0z(w3gohz,lir0o1) + gq815b, mu4ygka, lqhm2g)
      if(.not.(lqhm2g .ne. 1))goto 23472
23472 continue
      call vdgam1(gq815b, zixm0o, lqhm2g)
      if(.not.(lqhm2g .ne. 1))goto 23474
23474 continue
      da51l0o = mu4ygka - zixm0o - (jmwo0z(w3gohz,lir0o1) + gq815b) / (
     &q1znur + gq815b) + 1.0d0 + dlog(gq815b / (q1znur + gq815b))
      p4je8 = gq815b
      qxvi5(1,1) = gq815b
      xkcm3b(1,1) = q1znur
      nh2qxl = 5000
      call enbin9(qpzx6l, qxvi5, xkcm3b, ogq67o, uxzze7, c4uxow, uxzze7,
     & hmr3dx, jftq1, nh2qxl)
      if(.not.(c4uxow .ne. 1))goto 23476
      gqxvz8 = 5
      return
23476 continue
      yqsco4 = -qpzx6l(1,1) - 1.0d0 / gq815b + 1.0d0 / (gq815b + q1znur)
      lj4dph(w3gohz,2*lir0o1-1) = w8xfic(w3gohz,1) * q1znur * gq815b / (
     &q1znur + gq815b)
      lj4dph(w3gohz,2*lir0o1 ) = w8xfic(w3gohz,1) * gq815b * (-qpzx6l(1,
     &1)*gq815b - 1.0d0 + gq815b / (gq815b + q1znur))
      if(.not.(lj4dph(w3gohz,2*lir0o1-1) .lt. fiumb4))goto 23478
      lj4dph(w3gohz,2*lir0o1-1) = fiumb4
      jrxg6l(2*lir0o1-1,w3gohz) = zl11l0
      goto 23479
23478 continue
      jrxg6l(2*lir0o1-1,w3gohz) = dsqrt(lj4dph(w3gohz,2*lir0o1-1))
23479 continue
      if(.not.(lj4dph(w3gohz,2*lir0o1) .lt. fiumb4))goto 23480
      lj4dph(w3gohz,2*lir0o1) = fiumb4
      jrxg6l(2*lir0o1,w3gohz) = zl11l0
      goto 23481
23480 continue
      jrxg6l(2*lir0o1,w3gohz) = dsqrt(lj4dph(w3gohz,2*lir0o1))
23481 continue
      hr83e(w3gohz,2*lir0o1-1) = go0l1q(2*lir0o1-1,w3gohz) + jmwo0z(
     &w3gohz,lir0o1) / q1znur - 1.0d0
      hr83e(w3gohz,2*lir0o1 ) = go0l1q(2*lir0o1 ,w3gohz) + da51l0o / (
     &p4je8 * yqsco4)
23464 continue
23462 continue
      if(.not.(zvxw1l .eq. 8))goto 23482
      do 23484 w3gohz=1,nfiumb4 
      lj4dph(w3gohz,lir0o1) = w8xfic(w3gohz,1)
      jrxg6l(lir0o1,w3gohz) = dsqrt(lj4dph(w3gohz,lir0o1))
      hr83e(w3gohz,lir0o1) = jmwo0z(w3gohz,lir0o1)
23484 continue
23482 continue
      if(.not.(h3mrfq .eq. 1))goto 23486
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23488
      do 23490 w3gohz=1,nfiumb4 
      hr83e(w3gohz,2*lir0o1-1) = hr83e(w3gohz,2*lir0o1-1) - q121lc(
     &w3gohz)
23490 continue
      goto 23489
23488 continue
      do 23492 w3gohz=1,nfiumb4 
      hr83e(w3gohz,lir0o1) = hr83e(w3gohz,lir0o1) - q121lc(w3gohz)
23492 continue
23489 continue
23486 continue
      return
      end
      subroutine cqo2f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, 
     &w5poyv, hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, 
     &nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, 
     &vvl1li, nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx)
      implicit logical (a-z)
      integer p1i8xz(18), zqve1l(1), vvl1li(1)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, i83h1(1)
      double precision hft28(nfiumb4,1), jmwo0z(nfiumb4,aqk377), oju3yh(
     &nfiumb4,9), w8xfic(nfiumb4,1), go0l1q(lku8xq,nfiumb4), q121lc(
     &nfiumb4), w5poyv(aqk377,nfiumb4)
      double precision hr83e(nfiumb4,lku8xq), lj4dph(nfiumb4,lku8xq), 
     &jrxg6l(zkjqhi,nfiumb4), ur73jo(vi231l,1)
      double precision ioqzvb(vi231l,1), i0qvzl(1), nx1bat, lq8reh(1), 
     &l1zvxx(4)
      double precision t5vlzq(lku8xq,nfiumb4,2), zxao0o(lku8xq*(lku8xq+
     &1))
      integer w3gohz, d9rjek, nd6mep, i5uvkm, ptr, opf6cv, i2, oht3ga, 
     &ucgi1r, w5tcfp, cqui1v, xhe4cg, ugsma5, zvxw1l, pga6nul
      integer tvyd2b, fjg0qv, zx1610, zxiwf1, cr8hav, dyt0pg, h3mrfq, 
     &uvnk0i
      integer uxzze7, foej1u
      double precision hq710, fiumb4, t7sbea, xmr7cj, elq2cs, qik6ym, 
     &zl11l0, wlkaa3, jftq1
      double precision epx9jf1, epx9jf2
      integer scvgce
      uxzze7 = 1
      scvgce = 0
      oju3yh(1,1) = 1
      zxao0o(1) = 0.0d0
      cqui1v = p1i8xz(1)
      w5tcfp = p1i8xz(2)
      zxiwf1 = p1i8xz(3)
      xhe4cg = p1i8xz(4)
      ugsma5 = p1i8xz(5)
      zvxw1l = p1i8xz(6)
      pga6nul = p1i8xz(7)
      p1i8xz(9) = 0
      cr8hav = p1i8xz(11)
      dyt0pg = p1i8xz(12)
      h3mrfq = p1i8xz(14)
      uvnk0i = p1i8xz(15)
      foej1u = p1i8xz(18)
      fiumb4 = l1zvxx(1)
      zl11l0 = dsqrt(fiumb4)
      if(.not.((zvxw1l .eq. 1) .or. (zvxw1l .eq. 4)))goto 23494
      wlkaa3 = dlog(fiumb4)
23494 continue
      qik6ym = l1zvxx(2)
      jftq1 = l1zvxx(3)
      elq2cs = 0.0d0
      oht3ga = 0
      gqxvz8 = 1
      call qh4ulb(zqve1l, vvl1li, cqui1v)
      i5uvkm = cqui1v * (cqui1v+1) / 2
      call p0lk40(hft28, ur73jo, nfiumb4, lku8xq, vi231l, cqui1v, 
     &zvxw1l, aqk377, w5tcfp, cr8hav, i5uvkm, zqve1l, vvl1li, h3mrfq, 
     &q121lc)
653   epx9jf2 = 1.0d0
      if(.not.(ugsma5 .eq. 0))goto 23496
      do 23498 d9rjek=1,aqk377 
      call nbq4ua(jmwo0z, go0l1q, l1zvxx, nfiumb4, lku8xq, aqk377, 
     &zvxw1l, d9rjek, w8xfic, foej1u)
23498 continue
      goto 23497
23496 continue
      if(.not.(ugsma5 .eq. 2))goto 23500
      call pjw1l(ur73jo, lq8reh, go0l1q, nfiumb4, lku8xq, vi231l, 
     &zxiwf1, dyt0pg, oht3ga, zvxw1l, h3mrfq, q121lc)
23500 continue
23497 continue
      call o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l, 
     &oht3ga)
      if(.not.(ugsma5 .eq. 2))goto 23502
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, hq710, oht3ga, fiumb4, wlkaa3, uxzze7)
      goto 23503
23502 continue
      hq710 = -1.0d0
23503 continue
      do 23504 ucgi1r=1,pga6nul 
      do 23506 d9rjek=1,aqk377 
      call kqsxz1(jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, lj4dph, jrxg6l,
     & jftq1, fiumb4, zl11l0, nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, 
     &d9rjek, zvxw1l, gqxvz8, h3mrfq, q121lc)
23506 continue
      do 23508 d9rjek=1,zxiwf1 
      do 23510 w3gohz=1,vi231l 
      ioqzvb(w3gohz,d9rjek) = ur73jo(w3gohz,d9rjek)
23510 continue
23508 continue
      do 23512 d9rjek=1,zxiwf1 
      ptr = 1
      do 23514 opf6cv=1,nfiumb4 
      do 23516 i2=1,lku8xq 
      ioqzvb(ptr,d9rjek) = jrxg6l(i2,opf6cv) * ioqzvb(ptr,d9rjek)
      ptr = ptr + 1
23516 continue
23514 continue
23512 continue
      do 23518 nd6mep=1,zxiwf1 
      i83h1(nd6mep) = nd6mep
23518 continue
      t7sbea = 1.0d-7
      call dhkt9w(ioqzvb,vi231l,vi231l,zxiwf1,i0qvzl,i83h1,t5vlzq,
     &zx1610,t7sbea)
      if(.not.(zx1610 .ne. zxiwf1))goto 23520
      gqxvz8 = 2
      return
23520 continue
      do 23522 w3gohz=1,nfiumb4 
      do 23524 d9rjek=1,lku8xq 
      t5vlzq(d9rjek,w3gohz,1) = jrxg6l(d9rjek,w3gohz) * hr83e(w3gohz,
     &d9rjek)
23524 continue
23522 continue
      tvyd2b = 101
      call vdqrsl(ioqzvb,vi231l,vi231l,zx1610,i0qvzl, t5vlzq, elq2cs, 
     &t5vlzq(1,1,2), lq8reh, elq2cs,go0l1q,tvyd2b,fjg0qv)
      do 23526 w3gohz=1,nfiumb4 
      do 23528 d9rjek=1,lku8xq 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) / jrxg6l(d9rjek,
     &w3gohz)
23528 continue
23526 continue
      if(.not.(h3mrfq .eq. 1))goto 23530
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23532
      do 23534 w3gohz=1,nfiumb4 
      do 23536 d9rjek=1,aqk377 
      go0l1q(2*d9rjek-1,w3gohz) = go0l1q(2*d9rjek-1,w3gohz) + q121lc(
     &w3gohz)
23536 continue
23534 continue
      goto 23533
23532 continue
      do 23538 w3gohz=1,nfiumb4 
      do 23540 d9rjek=1,lku8xq 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) + q121lc(w3gohz)
23540 continue
23538 continue
23533 continue
23530 continue
      call o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l, 
     &oht3ga)
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, nx1bat,oht3ga,fiumb4,wlkaa3, uxzze7)
      xmr7cj = dabs(nx1bat - hq710) / (1.0d0 + dabs(nx1bat))
      if(.not.(xmr7cj .lt. qik6ym))goto 23542
      gqxvz8 = 0
      p1i8xz(8) = ucgi1r
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23544
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, nx1bat,oht3ga,fiumb4,wlkaa3, oht3ga)
23544 continue
      scvgce = 1
      goto 20097
      goto 23543
23542 continue
      hq710 = nx1bat
      scvgce = 0
23543 continue
23504 continue
20097 epx9jf1 = 0.0d0
      if(.not.(scvgce .eq. 1))goto 23546
      return
23546 continue
      if(.not.(ugsma5 .eq. 1 .or. ugsma5 .eq. 2))goto 23548
      ugsma5 = 0
      p1i8xz(9) = 1
      goto 653
23548 continue
      gqxvz8 = 3
      return
      end
      subroutine cqo1f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, 
     &w5poyv, hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, 
     &nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, 
     &vvl1li, nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx)
      implicit logical (a-z)
      integer p1i8xz(18), zqve1l(1), vvl1li(1)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, i83h1(1)
      double precision hft28(nfiumb4,1), jmwo0z(nfiumb4,aqk377), w8xfic(
     &nfiumb4,1), go0l1q(lku8xq,nfiumb4), q121lc(nfiumb4), w5poyv(
     &aqk377,nfiumb4), oju3yh(nfiumb4,9), hr83e(nfiumb4,lku8xq), lj4dph(
     &nfiumb4,lku8xq), jrxg6l(zkjqhi,nfiumb4), ur73jo(vi231l,1)
      double precision ioqzvb(vi231l,1), i0qvzl(1), nx1bat, lq8reh(1), 
     &l1zvxx(4)
      double precision t5vlzq(vi231l,3), zxao0o(lku8xq*(lku8xq+1))
      integer w3gohz, lir0o1, i5uvkm, oht3ga, ucgi1r, w5tcfp, h3mrfq, 
     &cqui1v, xhe4cg, ugsma5, zvxw1l, pga6nul
      integer tvyd2b, fjg0qv, zx1610, zxiwf1, dyt0pg, uvnk0i
      integer uxzze7, p1, foej1u
      double precision hq710, fiumb4, t7sbea, xmr7cj, elq2cs, qik6ym, 
     &zl11l0, wlkaa3, jftq1
      integer nd6mep
      double precision ni1qfp, epx9jf
      ni1qfp = 0.0d0
      uxzze7 = 1
      zxao0o(1) = 1.0d0
      call intpr(
     &"entering cqo1f uxzze7 -------------------------------: ",-1,
     &uxzze7,1)
      call intpr("in cqo1f aqk377: ",-1,aqk377,1)
      cqui1v = p1i8xz(1)
      w5tcfp = p1i8xz(2)
      zxiwf1 = p1i8xz(3)
      xhe4cg = p1i8xz(4)
      ugsma5 = p1i8xz(5)
      zvxw1l = p1i8xz(6)
      pga6nul = p1i8xz(7)
      p1i8xz(9) = 0
      dyt0pg = p1i8xz(12)
      if(.not.(dyt0pg .ne. 1))goto 23550
      gqxvz8 = 4
      return
23550 continue
      h3mrfq = p1i8xz(14)
      uvnk0i = p1i8xz(15)
      p1 = p1i8xz(16)
      foej1u = p1i8xz(18)
      call intpr("Entry to cqo1f: ugsma5 ",-1,ugsma5,1)
      fiumb4 = l1zvxx(1)
      zl11l0 = dsqrt(fiumb4)
      if(.not.((zvxw1l .eq. 1) .or. (zvxw1l .eq. 4)))goto 23552
      wlkaa3 = dlog(fiumb4)
23552 continue
      qik6ym = l1zvxx(2)
      jftq1 = l1zvxx(3)
      elq2cs = 0.0d0
      oht3ga = 0
      gqxvz8 = 1
      call qh4ulb(zqve1l, vvl1li, cqui1v)
      i5uvkm = cqui1v * (cqui1v+1) / 2
      call u16zxj(hft28, ur73jo, nfiumb4, cqui1v, zvxw1l, q121lc, 
     &vi231l, zxiwf1, i5uvkm, zqve1l, vvl1li, oju3yh, p1, h3mrfq)
      call dblepr("cqo1f: q121lc()",-1,q121lc,nfiumb4)
      call dblepr("cqo1f: ur73jo(,)",-1,ur73jo,vi231l*zxiwf1)
      call dblepr("cqo1f: w8xfic(,1)",-1,w8xfic(1,1),nfiumb4)
      do 23554 lir0o1=1,aqk377 
      call intpr("cqo1f: lir0o1======================: ",-1,lir0o1,1)
653   epx9jf = 1.0d0
      if(.not.(ugsma5 .eq. 0))goto 23556
      call intpr("cqo1f: calling nbq4ua ",-1,lir0o1,1)
      call nbq4ua(jmwo0z, go0l1q, l1zvxx, nfiumb4, lku8xq, aqk377, 
     &zvxw1l, lir0o1, w8xfic, foej1u)
      goto 23557
23556 continue
      if(.not.(ugsma5 .eq. 2))goto 23558
      call intpr("cqo1f: calling pjw1l; dyt0pg== ",-1,dyt0pg,1)
      call pjw1l(ur73jo, lq8reh(1+(lir0o1-1)*zxiwf1), go0l1q, nfiumb4, 
     &lku8xq, vi231l, zxiwf1, dyt0pg, lir0o1, zvxw1l, h3mrfq, q121lc)
23558 continue
23557 continue
      call o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l, 
     &lir0o1)
      if(.not.(ugsma5 .eq. 2))goto 23560
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, hq710, lir0o1, fiumb4, wlkaa3, uxzze7)
      goto 23561
23560 continue
      hq710 = -1.0d0
23561 continue
      do 23562 ucgi1r=1,pga6nul 
      call intpr("ucgi1r: ",-1,ucgi1r,1)
      call intpr("posn 7: ",-1,uxzze7,1)
      call intpr("zvxw1l: ",-1,zvxw1l,1)
      call dblepr("hq710",-1,hq710,1)
      call kqsxz1(jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, lj4dph, jrxg6l,
     & jftq1, fiumb4, zl11l0, nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, 
     &lir0o1, zvxw1l, gqxvz8, h3mrfq, q121lc)
      call dblepr("cqo1f: go0l1q",-1,go0l1q,lku8xq*nfiumb4)
      call dblepr("cqo1f: jrxg6l",-1,jrxg6l,zkjqhi*nfiumb4)
      call dblepr("cqo1f: hr83e",-1,hr83e,nfiumb4*lku8xq)
      call dblepr("cqo1f: lj4dph",-1,lj4dph,nfiumb4*lku8xq)
      do 23564 nd6mep=1,zxiwf1 
      do 23566 w3gohz=1,vi231l 
      ioqzvb(w3gohz,nd6mep) = ur73jo(w3gohz,nd6mep)
23566 continue
23564 continue
      call intpr("posn 3: ",-1,uxzze7,1)
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23568
      do 23570 nd6mep=1,zxiwf1 
      do 23572 w3gohz=1,nfiumb4 
      ioqzvb(2*w3gohz-1,nd6mep) = jrxg6l(2*lir0o1-1,w3gohz) * ioqzvb(2*
     &w3gohz-1,nd6mep)
      ioqzvb(2*w3gohz ,nd6mep) = jrxg6l(2*lir0o1 ,w3gohz) * ioqzvb(2*
     &w3gohz ,nd6mep)
23572 continue
23570 continue
      goto 23569
23568 continue
      do 23574 nd6mep=1,zxiwf1 
      do 23576 w3gohz=1,nfiumb4 
      ioqzvb(w3gohz,nd6mep) = jrxg6l(lir0o1,w3gohz) * ioqzvb(w3gohz,
     &nd6mep)
23576 continue
23574 continue
23569 continue
      call intpr("posn 4: ",-1,uxzze7,1)
      do 23578 nd6mep=1,zxiwf1 
      i83h1(nd6mep) = nd6mep
23578 continue
      call dblepr("cqo1f: ioqzvb",-1,ioqzvb,vi231l*zxiwf1)
      call intpr("ucgi1r: ",-1,ucgi1r,1)
      t7sbea = 1.0d-7
      call dhkt9w(ioqzvb,vi231l,vi231l,zxiwf1,i0qvzl,i83h1,t5vlzq,
     &zx1610,t7sbea)
      call intpr("i83h1: ",-1,i83h1,zxiwf1)
      if(.not.(zx1610 .ne. zxiwf1))goto 23580
      gqxvz8 = 2
      return
23580 continue
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23582
      do 23584 w3gohz=1,nfiumb4 
      t5vlzq(2*w3gohz-1,1) = jrxg6l(2*lir0o1-1,w3gohz) * hr83e(w3gohz,2*
     &lir0o1-1)
      t5vlzq(2*w3gohz ,1) = jrxg6l(2*lir0o1 ,w3gohz) * hr83e(w3gohz,2*
     &lir0o1 )
23584 continue
      goto 23583
23582 continue
      do 23586 w3gohz=1,nfiumb4
      t5vlzq(w3gohz,1) = jrxg6l(lir0o1,w3gohz) * hr83e(w3gohz,lir0o1)
23586 continue
23583 continue
      call intpr("posn 5: ",-1,uxzze7,1)
      tvyd2b = 101
      call intpr("posn 6: ",-1,uxzze7,1)
      call vdqrsl(ioqzvb,vi231l,vi231l,zx1610,i0qvzl, t5vlzq(1,1), 
     &elq2cs, t5vlzq(1,2), lq8reh(1+(lir0o1-1)*zxiwf1), elq2cs,t5vlzq(1,
     &3),tvyd2b,fjg0qv)
      call dblepr("lq8reh(1+(lir0o1-1)*zxiwf1)",-1,lq8reh(1+(lir0o1-1)*
     &zxiwf1),zxiwf1)
      if(.not.(uvnk0i .gt. 1))goto 23588
23588 continue
      do 23590 nd6mep=1,zxiwf1 
      t5vlzq(nd6mep,1) = lq8reh((lir0o1-1)*zxiwf1 + nd6mep)
23590 continue
      do 23592 nd6mep=1,zxiwf1 
      lq8reh((lir0o1-1)*zxiwf1 + i83h1(nd6mep)) = t5vlzq(nd6mep,1)
23592 continue
      call intpr("posn 7: ",-1,uxzze7,1)
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23594
      do 23596 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = t5vlzq(2*w3gohz-1,3) / jrxg6l(2*
     &lir0o1-1,w3gohz)
      go0l1q(2*lir0o1 ,w3gohz) = t5vlzq(2*w3gohz ,3) / jrxg6l(2*lir0o1 ,
     &w3gohz)
23596 continue
      if(.not.(h3mrfq .eq. 1))goto 23598
      do 23600 w3gohz=1,nfiumb4 
      go0l1q(2*lir0o1-1,w3gohz) = go0l1q(2*lir0o1-1,w3gohz) + q121lc(
     &w3gohz)
23600 continue
23598 continue
      goto 23595
23594 continue
      do 23602 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = t5vlzq(w3gohz,3) / jrxg6l(lir0o1,w3gohz)
23602 continue
      if(.not.(h3mrfq .eq. 1))goto 23604
      do 23606 w3gohz=1,nfiumb4 
      go0l1q(lir0o1,w3gohz) = go0l1q(lir0o1,w3gohz) + q121lc(w3gohz)
23606 continue
23604 continue
23595 continue
      call intpr("posn 8: ",-1,uxzze7,1)
      call o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l, 
     &lir0o1)
      call intpr("posn 8b: ",-1,uxzze7,1)
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, nx1bat,lir0o1,fiumb4,wlkaa3,uxzze7)
      call intpr("posn 8c: ",-1,uxzze7,1)
      xmr7cj = dabs(nx1bat - hq710) / (1.0d0 + dabs(nx1bat))
      call intpr("cqo1f: ucgi1r -------------",-1,ucgi1r,1)
      call dblepr("cqo1f: xmr7cj",-1,xmr7cj,1)
      if(.not.(xmr7cj .lt. qik6ym))goto 23608
      gqxvz8 = 0
      p1i8xz(8)=ucgi1r
      call intpr("cqo1f p1i8xz(8): ",-1,p1i8xz(8),1)
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23610
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, nx1bat,lir0o1,fiumb4,wlkaa3, oht3ga)
23610 continue
      ni1qfp = ni1qfp + nx1bat
      goto 1011
      goto 23609
23608 continue
      hq710 = nx1bat
23609 continue
      call intpr("posn 9: ",-1,uxzze7,1)
23562 continue
      call intpr("cqo1f; unsuccessful convergence: ",-1,uxzze7,1)
      if(.not.(ugsma5 .eq. 1))goto 23612
      ugsma5 = 0
      p1i8xz(9) = 1
      goto 653
23612 continue
      gqxvz8 = 3
1011  epx9jf = 1.0d0
23554 continue
      call intpr(
     &"exiting cqo1f uxzze7 ============================ : ",-1,uxzze7,
     &1)
      nx1bat = ni1qfp
      return
      end
      subroutine vcao6f(hft28, jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, 
     &lj4dph, jrxg6l, ioqzvb, i0qvzl, i83h1, nfiumb4, lku8xq, aqk377, 
     &vi231l, zkjqhi, gqxvz8, p1i8xz, nx1bat, lq8reh, t5vlzq, zxao0o, 
     &l1zvxx, gqai81,h2mzlo, sq5cvf, ynk9ah, uxs1iq, vliac4, vfd2pw,
     &sazp9g,s0, zrcbl2, nyg3mt, e6tljz, ifo4ew, ozuw3p, hwi2tb, nbd5rl,
     & wj5shg, ykdc2t, wk2, wzxao0o, phqco4, vb81l0, bmb, rjcq9o, mwk, 
     &n1zwoi, j1l0o1, qc7zyb, vlni8d, jko0o1, mnh3up, fg3pxq)
      implicit logical (a-z)
      integer p1i8xz(19)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, i83h1(1)
      double precision hft28(nfiumb4,1), jmwo0z(nfiumb4,aqk377), w8xfic(
     &nfiumb4,1), go0l1q(lku8xq,nfiumb4), w5poyv(aqk377,nfiumb4)
      double precision hr83e(nfiumb4,lku8xq), lj4dph(nfiumb4,lku8xq), 
     &jrxg6l(zkjqhi,nfiumb4)
      double precision ioqzvb(vi231l,2), i0qvzl(1), nx1bat, lq8reh(1), 
     &l1zvxx(6)
      double precision t5vlzq(vi231l,3), zxao0o(lku8xq*(lku8xq+1))
      integer lir0o1, sglfr1, oht3ga, ucgi1r, cqui1v, xhe4cg, ugsma5, 
     &zvxw1l, pga6nul
      integer dyt0pg, uvnk0i, zxiwf1
      integer uxzze7, c3qxjo
      double precision hq710, fiumb4, xmr7cj, elq2cs, qik6ym, zl11l0, 
     &wlkaa3, jftq1
      double precision ni1qfp, epx9jf
      integer gqai81(15), h2mzlo, ynk9ah(1),uxs1iq(1),vliac4(1), ozuw3p(
     &1), hwi2tb(3), nbd5rl(1), wj5shg(1)
      integer foej1u, jko0o1(1), mnh3up(1), fg3pxq(2), vlni8d(2)
      double precision sq5cvf(aqk377)
      double precision vfd2pw(h2mzlo,nfiumb4), sazp9g(nfiumb4,1),s0(
     &lku8xq), zrcbl2(h2mzlo,nfiumb4), nyg3mt(h2mzlo,nfiumb4), e6tljz(
     &nfiumb4,2), ifo4ew(h2mzlo,1), ykdc2t(1), wk2(nfiumb4,h2mzlo), 
     &phqco4(1), vb81l0(1), bmb(1), rjcq9o(1), mwk(1), j1l0o1(1), 
     &qc7zyb(1)
      integer ymetu2
      integer w3gohz, myx3od, ibd3vc
      integer d8gwha, tiav4e
      double precision purf2k(2), x1boaf, ad3xzo, j6gbnx, rk3jet
      double precision h4fgoy, das4bx
      double precision q121lc(2)
      x1boaf=0.0d0
      ad3xzo=0.0d0
      j6gbnx=0.0d0
      rk3jet=0.0d0
      d8gwha = p1i8xz(19)
      ni1qfp = 0.0d0
      uxzze7 = 1
      zxao0o(1) = 1.0d0
      t5vlzq(1,1) = 1.0d0
      cqui1v = p1i8xz(1)
      zxiwf1 = p1i8xz(3)
      xhe4cg = p1i8xz(4)
      ugsma5 = p1i8xz(5)
      zvxw1l = p1i8xz(6)
      pga6nul = p1i8xz(7)
      p1i8xz(9) = 0
      tiav4e = p1i8xz(11)
      dyt0pg = p1i8xz(12)
      if(.not.((dyt0pg .ne. 1) .or. (tiav4e .ne. cqui1v)))goto 23614
      gqxvz8 = 4
      return
23614 continue
      uvnk0i = p1i8xz(15)
      foej1u = p1i8xz(18)
      h4fgoy = l1zvxx(3+aqk377+aqk377+2)
      fiumb4 = l1zvxx(1)
      zl11l0 = dsqrt(fiumb4)
      if(.not.((zvxw1l .eq. 1) .or. (zvxw1l .eq. 4)))goto 23616
      wlkaa3 = dlog(fiumb4)
23616 continue
      qik6ym = l1zvxx(2)
      jftq1 = l1zvxx(3)
      elq2cs = 0.0d0
      oht3ga = 0
      gqxvz8 = 1
      do 23618 lir0o1=1,aqk377 
653   epx9jf = 1.0d0
      if(.not.(ugsma5 .eq. 0))goto 23620
      call nbq4ua(jmwo0z, go0l1q, l1zvxx, nfiumb4, lku8xq, aqk377, 
     &zvxw1l, lir0o1, w8xfic, foej1u)
      goto 23621
23620 continue
      if(.not.(ugsma5 .ne. 1))goto 23622
      gqxvz8 = 6
      return
23622 continue
23621 continue
      call o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l, 
     &lir0o1)
      if(.not.(ugsma5 .eq. 2))goto 23624
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, hq710, lir0o1, fiumb4, wlkaa3, uxzze7)
      goto 23625
23624 continue
      hq710 = -1.0d0
23625 continue
      do 23626 ucgi1r=1,pga6nul 
      call sptoq8(hft28, ioqzvb, nfiumb4, vi231l, cqui1v, zvxw1l)
      gqai81(7) = 0
      call kqsxz1(jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, lj4dph, jrxg6l,
     & jftq1, fiumb4, zl11l0, nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, 
     &lir0o1, zvxw1l, gqxvz8, oht3ga, q121lc)
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23628
      ymetu2 = 2*lir0o1-1
      goto 23629
23628 continue
      ymetu2 = lir0o1
23629 continue
      do 23630 myx3od=1,h2mzlo 
      do 23632 w3gohz=1,nfiumb4 
      zrcbl2(myx3od,w3gohz) = jrxg6l(ymetu2-1+myx3od,w3gohz)
      nyg3mt(myx3od,w3gohz) = go0l1q(ymetu2-1+myx3od,w3gohz)
23632 continue
23630 continue
      c3qxjo = tiav4e * aqk377
      sglfr1 = cqui1v * (lir0o1-1)
      if(.not.(ucgi1r .eq. 1))goto 23634
      x1boaf = sq5cvf( sglfr1 + hwi2tb(1))
      ad3xzo = sq5cvf(c3qxjo + sglfr1 + hwi2tb(1))
      if(.not.(cqui1v .eq. 2))goto 23636
      j6gbnx = sq5cvf( sglfr1 + hwi2tb(2))
      rk3jet = sq5cvf(c3qxjo + sglfr1 + hwi2tb(2))
23636 continue
      do 23638 myx3od=1,tiav4e 
      do 23640 w3gohz=1,nfiumb4 
      sazp9g(w3gohz,sglfr1 + hwi2tb(myx3od)) = 0.0d0
23640 continue
23638 continue
      goto 23635
23634 continue
      sq5cvf( sglfr1 + hwi2tb(1)) = x1boaf
      sq5cvf(c3qxjo + sglfr1 + hwi2tb(1)) = ad3xzo
      if(.not.(cqui1v .eq. 2))goto 23642
      sq5cvf( sglfr1 + hwi2tb(2)) = j6gbnx
      sq5cvf(c3qxjo + sglfr1 + hwi2tb(2)) = rk3jet
23642 continue
23635 continue
      call vbfa(d8gwha,nfiumb4,h2mzlo,gqai81, e6tljz, hr83e(1,ymetu2), 
     &lj4dph(1,ymetu2), sq5cvf( sglfr1 + hwi2tb(1)), sq5cvf(c3qxjo + 
     &sglfr1 + hwi2tb(1)), ynk9ah,uxs1iq,vliac4, vfd2pw,sazp9g(1,sglfr1 
     &+ hwi2tb(1)), nyg3mt,s0, lq8reh(1+(lir0o1-1)*zxiwf1), zo5jyl,
     &h4fgoy, ioqzvb,i0qvzl, i83h1, purf2k, zrcbl2, ifo4ew, ozuw3p, 
     &hwi2tb, nbd5rl, wj5shg, ykdc2t, wk2, wzxao0o, phqco4, vb81l0, bmb,
     & rjcq9o, mwk, n1zwoi, j1l0o1(1+(lir0o1-1)*(vlni8d(1+cqui1v)-1)), 
     &qc7zyb, das4bx, vlni8d, jko0o1, mnh3up, fg3pxq)
      l1zvxx(3+aqk377+aqk377+1) = das4bx
      ibd3vc = gqai81(14)
      if(.not.(ibd3vc .ne. 0))goto 23644
      call intpr("vcao6f: exiting because of an error",-1,ibd3vc,1)
      gqxvz8 = 8
      return
23644 continue
      do 23646 myx3od=1,h2mzlo 
      do 23648 w3gohz=1,nfiumb4 
      go0l1q(ymetu2-1+myx3od,w3gohz) = nyg3mt(myx3od,w3gohz)
23648 continue
23646 continue
      call o47zxq(go0l1q, w5poyv, nfiumb4, lku8xq, aqk377, zvxw1l, 
     &lir0o1)
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, nx1bat, lir0o1, fiumb4, wlkaa3, uxzze7)
      xmr7cj = dabs(nx1bat - hq710) / (1.0d0 + dabs(nx1bat))
      if(.not.(xmr7cj .lt. qik6ym))goto 23650
      gqxvz8 = 0
      p1i8xz(8) = ucgi1r
      if(.not.((zvxw1l .eq. 3) .or. (zvxw1l .eq. 5)))goto 23652
      call kqx20o(zvxw1l, jmwo0z, w8xfic, w5poyv, nfiumb4, lku8xq, 
     &aqk377, xhe4cg, go0l1q, nx1bat,lir0o1,fiumb4,wlkaa3, oht3ga)
23652 continue
      ni1qfp = ni1qfp + nx1bat
      goto 1011
      goto 23651
23650 continue
      hq710 = nx1bat
23651 continue
23626 continue
      if(.not.(ugsma5 .eq. 1))goto 23654
      ugsma5 = 0
      p1i8xz(9) = 1
      goto 653
23654 continue
      gqxvz8 = 3
1011  epx9jf = 1.0d0
23618 continue
      nx1bat = ni1qfp
      return
      end
      subroutine dcqof(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, 
     &w5poyv, hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, 
     &nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, 
     &vvl1li, nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx, ize5km, ip0ox8, 
     &v8gzsp, p2, xt3fko, o4xmfj, zclbn2)
      implicit logical (a-z)
      integer p1i8xz(19), zqve1l(1), vvl1li(1)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, i83h1(1)
      integer dyt0pg
      double precision hft28(nfiumb4,1), jmwo0z(nfiumb4,aqk377), oju3yh(
     &nfiumb4,9), w8xfic(nfiumb4,1), go0l1q(lku8xq,nfiumb4), q121lc(
     &nfiumb4), w5poyv(aqk377,nfiumb4), hr83e(nfiumb4,lku8xq), lj4dph(
     &nfiumb4,lku8xq), jrxg6l(zkjqhi,nfiumb4), ur73jo(vi231l,1)
      double precision ioqzvb(vi231l,1), i0qvzl(1), nx1bat, lq8reh(1), 
     &l1zvxx(4)
      double precision t5vlzq(lku8xq,nfiumb4,2), zxao0o(lku8xq*(lku8xq+
     &1))
      integer p2
      double precision ize5km(nfiumb4,p2), ip0ox8(nfiumb4,1), v8gzsp(p2,
     &1), xt3fko(p2,1), o4xmfj, zclbn2(1)
      integer w3gohz, s9otpy, pvnfr4, cqui1v, fiy4lc, nd6mep, z2q1li, 
     &foej1u
      double precision wrbc3q, gibj6t
      cqui1v = p1i8xz(1)
      fiy4lc = p1i8xz(5)
      dyt0pg = p1i8xz(12)
      z2q1li = p1i8xz(13)
      foej1u = p1i8xz(18)
      do 23656 pvnfr4=1,cqui1v 
      do 23658 w3gohz=1,nfiumb4 
      wrbc3q = 0.0d0
      do 23660 s9otpy=1,p2 
      wrbc3q = wrbc3q + ize5km(w3gohz,s9otpy) * v8gzsp(s9otpy,pvnfr4)
23660 continue
      ip0ox8(w3gohz,pvnfr4) = wrbc3q
      hft28(w3gohz,pvnfr4) = wrbc3q
23658 continue
23656 continue
      if(.not.(dyt0pg.eq.1))goto 23662
      call cqo1f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &gibj6t, zclbn2, t5vlzq, zxao0o, l1zvxx)
      goto 23663
23662 continue
      call cqo2f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &gibj6t, zclbn2, t5vlzq, zxao0o, l1zvxx)
23663 continue
      do 23664 s9otpy=1,p2 
      do 23666 w3gohz=1,nfiumb4 
      ize5km(w3gohz,s9otpy) = o4xmfj * ize5km(w3gohz,s9otpy)
23666 continue
23664 continue
      do 23668 pvnfr4=1,cqui1v 
      do 23670 s9otpy=1,p2 
      do 23672 w3gohz=1,nfiumb4 
      hft28(w3gohz,pvnfr4)=ip0ox8(w3gohz,pvnfr4)+ize5km(w3gohz,s9otpy)
23672 continue
      p1i8xz(5) = 2
      do 23674 nd6mep=1,z2q1li 
      lq8reh(nd6mep) = zclbn2(nd6mep)
23674 continue
      if(.not.(dyt0pg.eq.1))goto 23676
      call cqo1f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx)
      goto 23677
23676 continue
      call cqo2f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx)
23677 continue
      if(.not.(gqxvz8 .ne. 0))goto 23678
      return
23678 continue
      xt3fko(s9otpy,pvnfr4) = (nx1bat - gibj6t) / o4xmfj
23670 continue
      if(.not.(cqui1v .gt. 1))goto 23680
      do 23682 w3gohz=1,nfiumb4 
      hft28(w3gohz,pvnfr4) = ip0ox8(w3gohz,pvnfr4)
23682 continue
23680 continue
23668 continue
      p1i8xz(5) = fiy4lc
      return
      end
      subroutine vdcaof(hft28, jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, 
     &lj4dph, jrxg6l, ioqzvb, i0qvzl, i83h1, nfiumb4, lku8xq, aqk377, 
     &vi231l, zkjqhi, gqxvz8, p1i8xz, nx1bat, lq8reh, t5vlzq, zxao0o, 
     &l1zvxx, ize5km, ip0ox8, v8gzsp, p2, xt3fko, zclbn2, gqai81,h2mzlo,
     & sq5cvf, ynk9ah, uxs1iq, vliac4, vfd2pw,sazp9g,s0, zrcbl2, nyg3mt,
     & e6tljz, ifo4ew, ozuw3p, hwi2tb, nbd5rl, wj5shg, ykdc2t, wk2, 
     &wzxao0o, phqco4, vb81l0, bmb, rjcq9o, mwk, n1zwoi, j1l0o1, qc7zyb,
     & vlni8d, jko0o1, mnh3up, fg3pxq)
      implicit logical (a-z)
      integer p1i8xz(19)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, i83h1(1)
      integer dyt0pg
      double precision hft28(nfiumb4,1), jmwo0z(nfiumb4,aqk377), w8xfic(
     &nfiumb4,1), go0l1q(lku8xq,nfiumb4), w5poyv(aqk377,nfiumb4), hr83e(
     &nfiumb4,lku8xq), lj4dph(nfiumb4,lku8xq), jrxg6l(zkjqhi,nfiumb4)
      double precision ioqzvb(vi231l,1), i0qvzl(1), nx1bat, lq8reh(1), 
     &l1zvxx(6)
      double precision t5vlzq(lku8xq,nfiumb4,2)
      double precision zxao0o(lku8xq*(lku8xq+1))
      integer p2
      double precision ize5km(nfiumb4,p2), ip0ox8(nfiumb4,1), v8gzsp(p2,
     &1), xt3fko(p2,1), o4xmfj, zclbn2(1)
      integer w3gohz, pp, pvnfr4, cqui1v, fiy4lc, z2q1li, foej1u
      double precision wrbc3q, gibj6t
      integer gqai81(15), h2mzlo, ynk9ah(1),uxs1iq(1),vliac4(1), ozuw3p(
     &1), hwi2tb(1), nbd5rl(1), wj5shg(1), vlni8d(2), jko0o1(1), mnh3up(
     &1), fg3pxq(2)
      double precision sq5cvf(aqk377)
      double precision vfd2pw(h2mzlo,nfiumb4), sazp9g(nfiumb4,1),s0(
     &lku8xq), zrcbl2(h2mzlo,nfiumb4)
      double precision nyg3mt(h2mzlo,nfiumb4), e6tljz(nfiumb4,1), 
     &ifo4ew(h2mzlo,1), ykdc2t(1), wk2(nfiumb4,h2mzlo), phqco4(1), 
     &vb81l0(1), bmb(1), rjcq9o(1), mwk(1), j1l0o1(1), qc7zyb(1), 
     &das4bx
      integer d8gwha
      double precision h4fgoy
      das4bx = 0.0d0
      d8gwha = 0
      cqui1v = p1i8xz(1)
      fiy4lc = p1i8xz(5)
      dyt0pg = p1i8xz(12)
      z2q1li = p1i8xz(13)
      foej1u = p1i8xz(18)
      h4fgoy = l1zvxx(3+aqk377+aqk377+2)
      o4xmfj = l1zvxx(3+aqk377+aqk377+3)
      do 23684 pvnfr4=1,cqui1v 
      do 23686 w3gohz=1,nfiumb4 
      wrbc3q = 0.0d0
      do 23688 pp=1,p2 
      wrbc3q = wrbc3q + ize5km(w3gohz,pp) * v8gzsp(pp,pvnfr4)
23688 continue
      ip0ox8(w3gohz,pvnfr4) = wrbc3q
      hft28(w3gohz,pvnfr4) = wrbc3q
23686 continue
23684 continue
      if(.not.(dyt0pg.eq.1))goto 23690
      call vcao6f(hft28, jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, lj4dph, 
     &jrxg6l, ioqzvb, i0qvzl, i83h1, nfiumb4, lku8xq, aqk377, vi231l, 
     &zkjqhi, gqxvz8, p1i8xz, gibj6t, zclbn2, t5vlzq, zxao0o, l1zvxx, 
     &gqai81,h2mzlo, sq5cvf, ynk9ah, uxs1iq, vliac4, vfd2pw,sazp9g,s0, 
     &zrcbl2, nyg3mt, e6tljz, ifo4ew, ozuw3p, hwi2tb, nbd5rl, wj5shg, 
     &ykdc2t, wk2, wzxao0o, phqco4, vb81l0, bmb, rjcq9o, mwk, n1zwoi, 
     &j1l0o1, qc7zyb, vlni8d, jko0o1, mnh3up, fg3pxq)
      l1zvxx(3+aqk377+aqk377+1) = das4bx
      goto 23691
23690 continue
23691 continue
      do 23692 pp=1,p2 
      do 23694 w3gohz=1,nfiumb4 
      ize5km(w3gohz,pp) = o4xmfj * ize5km(w3gohz,pp)
23694 continue
23692 continue
      do 23696 pvnfr4=1,cqui1v 
      do 23698 pp=1,p2 
      do 23700 w3gohz=1,nfiumb4 
      hft28(w3gohz,pvnfr4) = ip0ox8(w3gohz,pvnfr4) + ize5km(w3gohz,pp)
23700 continue
      p1i8xz(5) = 0
      if(.not.(dyt0pg.eq.1))goto 23702
      call vcao6f(hft28, jmwo0z, w8xfic, go0l1q, w5poyv, hr83e, lj4dph, 
     &jrxg6l, ioqzvb, i0qvzl, i83h1, nfiumb4, lku8xq, aqk377, vi231l, 
     &zkjqhi, gqxvz8, p1i8xz, nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx, 
     &gqai81,h2mzlo, sq5cvf, ynk9ah, uxs1iq, vliac4, vfd2pw,sazp9g,s0, 
     &zrcbl2, nyg3mt, e6tljz, ifo4ew, ozuw3p, hwi2tb, nbd5rl, wj5shg, 
     &ykdc2t, wk2, wzxao0o, phqco4, vb81l0, bmb, rjcq9o, mwk, n1zwoi, 
     &j1l0o1, qc7zyb, vlni8d, jko0o1, mnh3up, fg3pxq)
      l1zvxx(3+aqk377+aqk377+1) = das4bx
      goto 23703
23702 continue
23703 continue
      if(.not.(gqxvz8 .ne. 0))goto 23704
      return
23704 continue
      xt3fko(pp,pvnfr4) = (nx1bat - gibj6t) / o4xmfj
23698 continue
      if(.not.(cqui1v .gt. 1))goto 23706
      do 23708 w3gohz=1,nfiumb4 
      hft28(w3gohz,pvnfr4) = ip0ox8(w3gohz,pvnfr4)
23708 continue
23706 continue
23696 continue
      p1i8xz(5) = fiy4lc
      return
      end
      subroutine duqof(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, 
     &w5poyv, hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, 
     &nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, 
     &vvl1li, nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx, ip0ox8, xt3fko, 
     &o4xmfj, zclbn2)
      implicit logical (a-z)
      integer p1i8xz(19), zqve1l(1), vvl1li(1)
      integer nfiumb4, lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, i83h1(1)
      integer dyt0pg
      double precision hft28(nfiumb4,1), jmwo0z(nfiumb4,aqk377), oju3yh(
     &nfiumb4,9), w8xfic(nfiumb4,1), go0l1q(lku8xq,nfiumb4), q121lc(
     &nfiumb4), w5poyv(aqk377,nfiumb4), hr83e(nfiumb4,lku8xq), lj4dph(
     &nfiumb4,lku8xq), jrxg6l(zkjqhi,nfiumb4), ur73jo(vi231l,1)
      double precision ioqzvb(vi231l,1), i0qvzl(1), nx1bat, lq8reh(1), 
     &l1zvxx(4)
      double precision t5vlzq(lku8xq,nfiumb4,2), zxao0o(lku8xq*(lku8xq+
     &1))
      double precision ip0ox8(nfiumb4,1), xt3fko(nfiumb4,1), o4xmfj, 
     &zclbn2(1)
      integer w3gohz, pvnfr4, cqui1v, fiy4lc, nd6mep, z2q1li
      double precision gibj6t
      cqui1v = p1i8xz(1)
      fiy4lc = p1i8xz(5)
      dyt0pg = p1i8xz(12)
      z2q1li = p1i8xz(13)
      if(.not.(dyt0pg.eq.1))goto 23710
      call cqo1f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &gibj6t, zclbn2, t5vlzq, zxao0o, l1zvxx)
      goto 23711
23710 continue
      call cqo2f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &gibj6t, zclbn2, t5vlzq, zxao0o, l1zvxx)
23711 continue
      do 23712 pvnfr4=1,cqui1v 
      do 23714 w3gohz=1,nfiumb4 
      hft28(w3gohz,pvnfr4) = ip0ox8(w3gohz,pvnfr4) + o4xmfj
      p1i8xz(5) = 2
      do 23716 nd6mep=1,z2q1li 
      lq8reh(nd6mep) = zclbn2(nd6mep)
23716 continue
      if(.not.(dyt0pg.eq.1))goto 23718
      call cqo1f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx)
      goto 23719
23718 continue
      call cqo2f(hft28, jmwo0z, oju3yh, w8xfic, go0l1q, q121lc, w5poyv, 
     &hr83e, lj4dph, jrxg6l, ur73jo, ioqzvb, i0qvzl, i83h1, nfiumb4, 
     &lku8xq, aqk377, vi231l, zkjqhi, gqxvz8, p1i8xz, zqve1l, vvl1li, 
     &nx1bat, lq8reh, t5vlzq, zxao0o, l1zvxx)
23719 continue
      if(.not.(gqxvz8 .ne. 0))goto 23720
      return
23720 continue
      xt3fko(w3gohz,pvnfr4) = (nx1bat - gibj6t) / o4xmfj
      hft28(w3gohz,pvnfr4) = ip0ox8(w3gohz,pvnfr4)
23714 continue
23712 continue
      p1i8xz(5) = fiy4lc
      return
      end
