      subroutine vbvs(nfiumb4,knot,j1l0o1,nk,p3vlea,ocaxi1,ikscn4,
     &lku8xq)
      integer nfiumb4, nk, ikscn4, lku8xq
      double precision knot(nk+4), j1l0o1(nk,lku8xq), p3vlea(nfiumb4), 
     &ocaxi1(nfiumb4,lku8xq)
      double precision xvalue, bvalue
      integer w3gohz, d9rjek, def4wn
      def4wn = 4
      do 23000 d9rjek=1,lku8xq 
      do 23002 w3gohz=1,nfiumb4 
      xvalue = p3vlea(w3gohz)
      ocaxi1(w3gohz,d9rjek) = bvalue(knot, j1l0o1(1,d9rjek), nk, def4wn,
     & xvalue, ikscn4)
23002 continue
23000 continue
      return
      end
      subroutine j3navf(nkzg2p, nk, lku8xq, a51l0o, l6xrjt, nf8brk)
      implicit logical (a-z)
      integer nk, lku8xq, a51l0o
      double precision nkzg2p(a51l0o,nk*lku8xq), l6xrjt(lku8xq), nf8brk(
     &nk,4)
      integer w3gohz, d9rjek
      do 23004 w3gohz=1,nk 
      do 23006 d9rjek=1,lku8xq 
      nkzg2p(a51l0o,(w3gohz-1)*lku8xq+d9rjek) = nkzg2p(a51l0o,(w3gohz-1)
     &*lku8xq+d9rjek) + l6xrjt(d9rjek) * nf8brk(w3gohz,1)
23006 continue
23004 continue
      do 23008 w3gohz=1,(nk-1) 
      do 23010 d9rjek=1,lku8xq 
      nkzg2p(a51l0o-lku8xq,(w3gohz-0)*lku8xq+d9rjek) = nkzg2p(a51l0o-
     &lku8xq,(w3gohz-0)*lku8xq+d9rjek) + l6xrjt(d9rjek) * nf8brk(w3gohz,
     &2)
23010 continue
23008 continue
      do 23012 w3gohz=1,(nk-2) 
      do 23014 d9rjek=1,lku8xq 
      nkzg2p(a51l0o-2*lku8xq,(w3gohz+1)*lku8xq+d9rjek) = nkzg2p(a51l0o-
     &2*lku8xq,(w3gohz+1)*lku8xq+d9rjek) + l6xrjt(d9rjek) * nf8brk(
     &w3gohz,3)
23014 continue
23012 continue
      do 23016 w3gohz=1,(nk-3) 
      do 23018 d9rjek=1,lku8xq 
      nkzg2p(a51l0o-3*lku8xq,(w3gohz+2)*lku8xq+d9rjek) = nkzg2p(a51l0o-
     &3*lku8xq,(w3gohz+2)*lku8xq+d9rjek) + l6xrjt(d9rjek) * nf8brk(
     &w3gohz,4)
23018 continue
23016 continue
      return
      end
      subroutine wgy5ta(p1rifj, s17te9, nbbad, uq9jtc, nkzg2p, w8xfic, 
     &evgfu3, anke8p, lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, 
     &vvl1li)
      implicit logical (a-z)
      integer p1rifj, s17te9, nbbad, evgfu3, anke8p, lku8xq, a51l0o, 
     &xhe4cg, nfiumb4, nk, zqve1l(1), vvl1li(1)
      double precision uq9jtc(4,1), nkzg2p(a51l0o, nk*lku8xq), w8xfic(
     &nfiumb4,xhe4cg)
      double precision temp
      integer xi1mqb, j3ymns, bcol, brow, y9eilo, pazyk8
      bcol = s17te9 + nbbad
      brow = s17te9
      do 23020 xi1mqb=1,xhe4cg 
      temp = w8xfic(p1rifj,xi1mqb) * uq9jtc(evgfu3,1) * uq9jtc(anke8p,1)
      y9eilo = (brow-1)*lku8xq + zqve1l(xi1mqb)
      pazyk8 = (bcol-1)*lku8xq + vvl1li(xi1mqb)
      j3ymns = pazyk8 - y9eilo
      nkzg2p(a51l0o-j3ymns, pazyk8) = nkzg2p(a51l0o-j3ymns, pazyk8) + 
     &temp
      if(.not.(nbbad .gt. 0 .and. vvl1li(xi1mqb) .ne. zqve1l(xi1mqb)))
     &goto 23022
      y9eilo = (brow-1)*lku8xq + vvl1li(xi1mqb)
      pazyk8 = (bcol-1)*lku8xq + zqve1l(xi1mqb)
      j3ymns = pazyk8 - y9eilo
      nkzg2p(a51l0o-j3ymns, pazyk8) = nkzg2p(a51l0o-j3ymns, pazyk8) + 
     &temp
23022 continue
23020 continue
      return
      end
      subroutine vsplin(p3vlea,lj4dph,w8xfic,nfiumb4,onyz6j, nk,a51l0o,
     &lku8xq,xhe4cg, zqve1l,vvl1li, zxao0o, l6xrjt, fjg0qv, w5poyv, 
     &j1l0o1, nkzg2p, cy3dhl, vb81l0, l8dgox, y6jcvk, bmb, rjcq9o, dof, 
     &sz6ohy, la5dcf, e5jrsg)
      implicit logical (a-z)
      integer nfiumb4, nk, a51l0o, lku8xq, xhe4cg, zqve1l(1), vvl1li(1),
     & fjg0qv, la5dcf, e5jrsg
      integer y6jcvk
      double precision p3vlea(nfiumb4), lj4dph(nfiumb4,lku8xq), w8xfic(
     &nfiumb4,xhe4cg), onyz6j(nk+4), zxao0o(lku8xq,lku8xq,16), l6xrjt(
     &lku8xq), w5poyv(nfiumb4,lku8xq), j1l0o1(nk,lku8xq), nkzg2p(a51l0o,
     &nk*lku8xq), cy3dhl(lku8xq,nk)
      double precision vb81l0(nk,lku8xq), l8dgox(e5jrsg,lku8xq), bmb(
     &lku8xq,lku8xq), rjcq9o(nfiumb4,lku8xq), dof(lku8xq), sz6ohy(1)
      integer d9rjek, w3gohz, m5xudf, i6ndbu, xi1mqb, rlhz2a
      integer yc1ezl, mk2vyr, thfyl1, hts1gp(3), ispar, opf6cv
      double precision kqoy6w, uq9jtc(4,1), z2djpt(16), egoxa3, n9peut, 
     &bt9lgm, v2isnf, fjo2dy, fpcb2n(3)
      do 23024 d9rjek=1,lku8xq 
      if(.not.(l6xrjt(d9rjek) .eq. 0.0d0))goto 23026
      ispar=0
      rlhz2a=3
      goto 23027
23026 continue
      ispar=1
      rlhz2a=1
23027 continue
      if(.not.((lku8xq .eq. 1) .or. (xhe4cg.eq.lku8xq) .or. (ispar .eq. 
     &0)))goto 23028
      mk2vyr = 4
      bt9lgm = 1.50d0
      v2isnf = -0.25d0
      thfyl1 = 1
      fjo2dy=0.001d0
      hts1gp(1) = rlhz2a
      hts1gp(2) = ispar
      hts1gp(3) = 300
      fpcb2n(1) = v2isnf
      fpcb2n(2) = bt9lgm
      fpcb2n(3) = fjo2dy
      yc1ezl=0
      if(.not.((lku8xq .eq. 1) .or. (xhe4cg.eq.lku8xq)))goto 23030
      do 23032 w3gohz=1,nfiumb4 
      lj4dph(w3gohz,d9rjek) = lj4dph(w3gohz,d9rjek) / w8xfic(w3gohz,
     &d9rjek)
23032 continue
      call nvhb7f(egoxa3, dof(d9rjek), p3vlea, lj4dph(1,d9rjek), w8xfic(
     &1,d9rjek), nfiumb4,nk, onyz6j,j1l0o1(1,d9rjek), w5poyv(1,d9rjek), 
     &rjcq9o(1,d9rjek), n9peut,l6xrjt(d9rjek),fpcb2n, sz6ohy, yc1ezl,
     &hts1gp, mk2vyr,thfyl1,la5dcf)
      if(.not.(la5dcf .ne. 0))goto 23034
      return
23034 continue
      do 23036 w3gohz=1,nfiumb4 
      w8xfic(w3gohz,d9rjek) = w8xfic(w3gohz,d9rjek) * w8xfic(w3gohz,
     &d9rjek)
23036 continue
      if(.not.(y6jcvk .ne. 0))goto 23038
      do 23040 w3gohz=1,nfiumb4 
      l8dgox(w3gohz,d9rjek) = rjcq9o(w3gohz,d9rjek) / w8xfic(w3gohz,
     &d9rjek)
23040 continue
23038 continue
      goto 23031
23030 continue
      call nvhb7f(egoxa3, dof(d9rjek), p3vlea, cy3dhl(1,d9rjek), w8xfic(
     &1,d9rjek), nfiumb4,nk, onyz6j,j1l0o1(1,d9rjek),w5poyv(1,d9rjek), 
     &rjcq9o(1,d9rjek), n9peut,l6xrjt(d9rjek),fpcb2n, sz6ohy, yc1ezl,
     &hts1gp, mk2vyr,thfyl1,la5dcf)
      if(.not.(la5dcf .ne. 0))goto 23042
      return
23042 continue
      do 23044 w3gohz=1,nfiumb4 
      w8xfic(w3gohz,d9rjek) = w8xfic(w3gohz,d9rjek) * w8xfic(w3gohz,
     &d9rjek)
23044 continue
23031 continue
      if(.not.(la5dcf .ne. 0))goto 23046
      return
23046 continue
23028 continue
23024 continue
      if(.not.((lku8xq .eq. 1) .or. (xhe4cg .eq. lku8xq)))goto 23048
      return
23048 continue
      do 23050 w3gohz=1,nk 
      do 23052 d9rjek=1,lku8xq 
      cy3dhl(d9rjek,w3gohz)=0.0d0
23052 continue
23050 continue
      do 23054 w3gohz=1,(nk*lku8xq) 
      do 23056 d9rjek=1,a51l0o 
      nkzg2p(d9rjek,w3gohz) = 0.0d0
23056 continue
23054 continue
      kqoy6w = 0.1d-9
      do 23058 w3gohz=1,nfiumb4 
      call vinterv(onyz6j(1),(nk+1),p3vlea(w3gohz),m5xudf,i6ndbu)
      if(.not.(i6ndbu .eq. 1))goto 23060
      if(.not.(p3vlea(w3gohz) .le. (onyz6j(m5xudf)+kqoy6w)))goto 23062
      m5xudf=m5xudf-1
      goto 23063
23062 continue
      return
23063 continue
23060 continue
      call vbsplvd(onyz6j,4,p3vlea(w3gohz),m5xudf,z2djpt,uq9jtc,1)
      d9rjek= m5xudf-4+1
      do 23064 xi1mqb=1,lku8xq 
      cy3dhl(xi1mqb,d9rjek)=cy3dhl(xi1mqb,d9rjek) + lj4dph(w3gohz,
     &xi1mqb) * uq9jtc(1,1)
23064 continue
      call wgy5ta(w3gohz, d9rjek, 0, uq9jtc, nkzg2p, w8xfic, 1, 1, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      call wgy5ta(w3gohz, d9rjek, 1, uq9jtc, nkzg2p, w8xfic, 1, 2, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      call wgy5ta(w3gohz, d9rjek, 2, uq9jtc, nkzg2p, w8xfic, 1, 3, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      call wgy5ta(w3gohz, d9rjek, 3, uq9jtc, nkzg2p, w8xfic, 1, 4, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      d9rjek= m5xudf-4+2
      do 23066 xi1mqb=1,lku8xq 
      cy3dhl(xi1mqb,d9rjek)=cy3dhl(xi1mqb,d9rjek) + lj4dph(w3gohz,
     &xi1mqb) * uq9jtc(2,1)
23066 continue
      call wgy5ta(w3gohz, d9rjek, 0, uq9jtc, nkzg2p, w8xfic, 2, 2, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      call wgy5ta(w3gohz, d9rjek, 1, uq9jtc, nkzg2p, w8xfic, 2, 3, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      call wgy5ta(w3gohz, d9rjek, 2, uq9jtc, nkzg2p, w8xfic, 2, 4, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      d9rjek= m5xudf-4+3
      do 23068 xi1mqb=1,lku8xq 
      cy3dhl(xi1mqb,d9rjek)=cy3dhl(xi1mqb,d9rjek) + lj4dph(w3gohz,
     &xi1mqb) * uq9jtc(3,1)
23068 continue
      call wgy5ta(w3gohz, d9rjek, 0, uq9jtc, nkzg2p, w8xfic, 3, 3, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      call wgy5ta(w3gohz, d9rjek, 1, uq9jtc, nkzg2p, w8xfic, 3, 4, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
      d9rjek= m5xudf-4+4
      do 23070 xi1mqb=1,lku8xq 
      cy3dhl(xi1mqb,d9rjek)=cy3dhl(xi1mqb,d9rjek) + lj4dph(w3gohz,
     &xi1mqb) * uq9jtc(4,1)
23070 continue
      call wgy5ta(w3gohz, d9rjek, 0, uq9jtc, nkzg2p, w8xfic, 4, 4, 
     &lku8xq, a51l0o, xhe4cg, nfiumb4, nk, zqve1l, vvl1li)
23058 continue
      call poqy8c(vb81l0(1,1), vb81l0(1,2), vb81l0(1,3), vb81l0(1,4), 
     &onyz6j, nk)
      call j3navf(nkzg2p, nk, lku8xq, a51l0o, l6xrjt, vb81l0)
      call vdpbfa7(nkzg2p, a51l0o, nk*lku8xq, a51l0o-1, fjg0qv, vb81l0)
      if(.not.(fjg0qv .ne. 0))goto 23072
      return
23072 continue
      call vdpbsl7(nkzg2p, a51l0o, nk*lku8xq, a51l0o-1, cy3dhl, vb81l0)
      opf6cv = 0
      do 23074 w3gohz=1,nk 
      do 23076 d9rjek=1,lku8xq 
      opf6cv = opf6cv + 1
      j1l0o1(w3gohz,d9rjek) = cy3dhl(d9rjek,w3gohz)
23076 continue
23074 continue
      call ye3zvn(onyz6j, p3vlea, j1l0o1, nfiumb4, nk, lku8xq, w5poyv)
      call gyzcj5(nkzg2p, nkzg2p, vb81l0, zxao0o, a51l0o-1, nk*lku8xq)
      call jiyw4z(nkzg2p, p3vlea, onyz6j, l8dgox, a51l0o, nfiumb4, nk, 
     &lku8xq, y6jcvk, bmb, zxao0o, w8xfic, rjcq9o, xhe4cg, zqve1l, 
     &vvl1li, e5jrsg)
      return
      end
      subroutine ye3zvn(knot, p3vlea, j1l0o1, nfiumb4, nk, lku8xq, 
     &w5poyv)
      implicit logical (a-z)
      integer nfiumb4, nk, lku8xq
      double precision knot(nk+4), p3vlea(nfiumb4), j1l0o1(nk,lku8xq), 
     &w5poyv(nfiumb4,lku8xq)
      double precision xvalue, bvalue
      integer w3gohz, d9rjek
      do 23078 w3gohz=1,nfiumb4 
      xvalue = p3vlea(w3gohz)
      do 23080 d9rjek=1,lku8xq 
      w5poyv(w3gohz,d9rjek) = bvalue(knot, j1l0o1(1,d9rjek), nk, 4, 
     &xvalue, 0)
23080 continue
23078 continue
      return
      end
      subroutine vsuff9(nfiumb4,uxs1iq,ynk9ah, p3vlea,jmwo0z,w8xfic, 
     &qxy6aj,bz3pyo,ax1cdp,f0pzmy,lg3zhr, lku8xq, xhe4cg, zkjqhi, 
     &zqve1l, vvl1li, bgu6fw, ve2mqu, ifo4ew, du8jbv, wj5shg, x6rito, 
     &c4uxow)
      implicit logical (a-z)
      integer nfiumb4, uxs1iq, ynk9ah(nfiumb4), lku8xq, xhe4cg, zkjqhi, 
     &zqve1l(1),vvl1li(1), du8jbv, wj5shg, x6rito, c4uxow
      double precision p3vlea(nfiumb4), jmwo0z(nfiumb4,lku8xq), w8xfic(
     &nfiumb4,xhe4cg), qxy6aj(uxs1iq), bz3pyo(uxs1iq,lku8xq), ax1cdp(
     &uxs1iq,1), f0pzmy(zkjqhi,uxs1iq), lg3zhr(uxs1iq,lku8xq), bgu6fw(
     &lku8xq,lku8xq+1), ve2mqu(du8jbv,du8jbv+1), ifo4ew(lku8xq,du8jbv)
      integer w3gohz, d9rjek, nd6mep, xi1mqb, i1nkrb, j0qwtz
      integer uxzze7
      uxzze7 = 1
      if(.not.(wj5shg .eq. 1))goto 23082
      if(.not.((zkjqhi .ne. xhe4cg) .or. (du8jbv .ne. lku8xq)))goto 2308
     &4
      c4uxow = 0
      return
23084 continue
23082 continue
      j0qwtz = lku8xq * (lku8xq+1) / 2
      if(.not.(xhe4cg .gt. j0qwtz))goto 23086
23086 continue
      call qh4ulb(zqve1l, vvl1li, lku8xq)
      do 23088 w3gohz=1,nfiumb4 
      qxy6aj(ynk9ah(w3gohz))=p3vlea(w3gohz)
23088 continue
      do 23090 d9rjek=1,lku8xq 
      do 23092 w3gohz=1,uxs1iq 
      lg3zhr(w3gohz,d9rjek) = 0.0d0
23092 continue
23090 continue
      do 23094 d9rjek=1,xhe4cg 
      do 23096 w3gohz=1,uxs1iq 
      ax1cdp(w3gohz,d9rjek) = 0.0d0
23096 continue
23094 continue
      if(.not.(xhe4cg .ne. j0qwtz))goto 23098
      do 23100 nd6mep=1,lku8xq 
      do 23102 d9rjek=1,lku8xq 
      bgu6fw(d9rjek,nd6mep) = 0.0d0
23102 continue
23100 continue
23098 continue
      do 23104 w3gohz=1,nfiumb4 
      do 23106 d9rjek=1,xhe4cg 
      bgu6fw(zqve1l(d9rjek),vvl1li(d9rjek)) = w8xfic(w3gohz,d9rjek)
      bgu6fw(vvl1li(d9rjek),zqve1l(d9rjek)) = bgu6fw(zqve1l(d9rjek),
     &vvl1li(d9rjek))
23106 continue
      do 23108 d9rjek=1,lku8xq 
      do 23110 nd6mep=1,lku8xq 
      lg3zhr(ynk9ah(w3gohz),d9rjek) = lg3zhr(ynk9ah(w3gohz),d9rjek) + 
     &bgu6fw(d9rjek,nd6mep)*jmwo0z(w3gohz,nd6mep)
23110 continue
23108 continue
      do 23112 d9rjek=1,xhe4cg 
      ax1cdp(ynk9ah(w3gohz),d9rjek) = ax1cdp(ynk9ah(w3gohz),d9rjek) + 
     &w8xfic(w3gohz,d9rjek)
23112 continue
23104 continue
      c4uxow = 1
      if(.not.(wj5shg .eq. 1))goto 23114
      do 23116 w3gohz=1,uxs1iq 
      do 23118 d9rjek=1,xhe4cg 
      bgu6fw(zqve1l(d9rjek),vvl1li(d9rjek)) = ax1cdp(w3gohz,d9rjek)
      bgu6fw(vvl1li(d9rjek),zqve1l(d9rjek)) = bgu6fw(zqve1l(d9rjek),
     &vvl1li(d9rjek))
23118 continue
      do 23120 d9rjek=1,lku8xq 
      bgu6fw(d9rjek,lku8xq+1)=lg3zhr(w3gohz,d9rjek)
23120 continue
      call vcholf(bgu6fw, bgu6fw(1,lku8xq+1), lku8xq, c4uxow, uxzze7)
      if(.not.(c4uxow .ne. 1))goto 23122
      return
23122 continue
      if(.not.(x6rito .ne. 0))goto 23124
      do 23126 d9rjek=1,xhe4cg 
      f0pzmy(d9rjek,w3gohz) = bgu6fw(zqve1l(d9rjek),vvl1li(d9rjek))
23126 continue
23124 continue
      do 23128 d9rjek=1,lku8xq 
      bz3pyo(w3gohz,d9rjek)=bgu6fw(d9rjek,lku8xq+1)
23128 continue
23116 continue
      goto 23115
23114 continue
      if(.not.(xhe4cg .ne. j0qwtz))goto 23130
      do 23132 d9rjek=1,lku8xq 
      do 23134 nd6mep=1,lku8xq 
      bgu6fw(d9rjek,nd6mep) = 0.0d0
23134 continue
23132 continue
23130 continue
      do 23136 w3gohz=1,uxs1iq 
      call qh4ulb(zqve1l, vvl1li, lku8xq)
      do 23138 d9rjek=1,xhe4cg 
      bgu6fw(zqve1l(d9rjek),vvl1li(d9rjek)) = ax1cdp(w3gohz,d9rjek)
      bgu6fw(vvl1li(d9rjek),zqve1l(d9rjek)) = bgu6fw(zqve1l(d9rjek),
     &vvl1li(d9rjek))
23138 continue
      do 23140 d9rjek=1,lku8xq 
      bgu6fw(d9rjek,lku8xq+1)=lg3zhr(w3gohz,d9rjek)
23140 continue
      do 23142 d9rjek=1,du8jbv 
      do 23144 nd6mep=d9rjek,du8jbv 
      ve2mqu(d9rjek,nd6mep) = 0.0d0
      do 23146 xi1mqb=1,lku8xq 
      do 23148 i1nkrb=1,lku8xq 
      ve2mqu(d9rjek,nd6mep) = ve2mqu(d9rjek,nd6mep) + ifo4ew(xi1mqb,
     &d9rjek) * bgu6fw(xi1mqb,i1nkrb) * ifo4ew(i1nkrb,nd6mep)
23148 continue
23146 continue
23144 continue
23142 continue
      call qh4ulb(zqve1l, vvl1li, du8jbv)
      do 23150 d9rjek=1,zkjqhi 
      ax1cdp(w3gohz,d9rjek) = ve2mqu(zqve1l(d9rjek),vvl1li(d9rjek))
23150 continue
      do 23152 d9rjek=1,du8jbv 
      ve2mqu(d9rjek,du8jbv+1) = 0.0d0
      do 23154 xi1mqb=1,lku8xq 
      ve2mqu(d9rjek,du8jbv+1) = ve2mqu(d9rjek,du8jbv+1) + ifo4ew(xi1mqb,
     &d9rjek) * bgu6fw(xi1mqb,lku8xq+1)
23154 continue
23152 continue
      do 23156 d9rjek=1,du8jbv 
      lg3zhr(w3gohz,d9rjek) = ve2mqu(d9rjek,du8jbv+1)
23156 continue
      call vcholf(ve2mqu, ve2mqu(1,du8jbv+1), du8jbv, c4uxow, uxzze7)
      if(.not.(c4uxow .ne. 1))goto 23158
      return
23158 continue
      if(.not.(x6rito .ne. 0))goto 23160
      do 23162 d9rjek=1,zkjqhi 
      f0pzmy(d9rjek,w3gohz) = ve2mqu(zqve1l(d9rjek),vvl1li(d9rjek))
23162 continue
23160 continue
      do 23164 d9rjek=1,du8jbv 
      bz3pyo(w3gohz,d9rjek) = ve2mqu(d9rjek,du8jbv+1)
23164 continue
23136 continue
23115 continue
      return
      end
      subroutine jiyw4z(n5fkml, p3vlea, onyz6j, svpr1i, a51l0o, nfiumb4,
     & nk, lku8xq, ifvar, bmb, bgu6fw, w8xfic, rjcq9o, xhe4cg, zqve1l, 
     &vvl1li, e5jrsg)
      implicit logical (a-z)
      integer a51l0o, nfiumb4, nk, lku8xq, ifvar, xhe4cg, zqve1l(1), 
     &vvl1li(1), e5jrsg
      double precision n5fkml(a51l0o,nk*lku8xq), p3vlea(nfiumb4), 
     &onyz6j(nk+4), svpr1i(e5jrsg,lku8xq), bmb(lku8xq,lku8xq), bgu6fw(
     &lku8xq,lku8xq), w8xfic(nfiumb4,xhe4cg), rjcq9o(nfiumb4,lku8xq)
      integer w3gohz, d9rjek, nd6mep, m5xudf, i6ndbu, xi1mqb, i1nkrb
      double precision kqoy6w, z2djpt(16), uq9jtc(4,1)
      if(.not.(ifvar .ne. 0))goto 23166
      do 23168 nd6mep=1,lku8xq 
      do 23170 w3gohz=1,nfiumb4 
      svpr1i(w3gohz,nd6mep) = 0.0d0
23170 continue
23168 continue
23166 continue
      kqoy6w = 0.10d-9
      call qh4ulb(zqve1l, vvl1li, lku8xq)
      do 23172 w3gohz=1,nfiumb4 
      do 23174 d9rjek=1,lku8xq 
      do 23176 nd6mep=1,lku8xq 
      bmb(d9rjek,nd6mep)=0.0d0
23176 continue
23174 continue
      call vinterv(onyz6j(1), (nk+1), p3vlea(w3gohz), m5xudf, i6ndbu)
      if(.not.(i6ndbu.eq. 1))goto 23178
      if(.not.(p3vlea(w3gohz) .le. (onyz6j(m5xudf)+kqoy6w)))goto 23180
      m5xudf=m5xudf-1
      goto 23181
23180 continue
      return
23181 continue
23178 continue
      call vbsplvd(onyz6j, 4, p3vlea(w3gohz), m5xudf, z2djpt, uq9jtc, 1)
      d9rjek= m5xudf-4+1
      do 23182 xi1mqb=d9rjek,d9rjek+3 
      call vsel(xi1mqb, xi1mqb, lku8xq, nk, a51l0o, n5fkml, bgu6fw)
      call bf7qci(lku8xq, uq9jtc(xi1mqb-d9rjek+1,1) * uq9jtc(xi1mqb-
     &d9rjek+1,1), bgu6fw, bmb)
23182 continue
      do 23184 xi1mqb=d9rjek,d9rjek+3 
      do 23186 i1nkrb=xi1mqb+1,d9rjek+3 
      call vsel(xi1mqb, i1nkrb, lku8xq, nk, a51l0o, n5fkml, bgu6fw)
      call bf7qci(lku8xq, 2.0d0 * uq9jtc(xi1mqb-d9rjek+1,1) * uq9jtc(
     &i1nkrb-d9rjek+1,1), bgu6fw, bmb)
23186 continue
23184 continue
      if(.not.(ifvar .ne. 0))goto 23188
      do 23190 d9rjek=1,lku8xq 
      svpr1i(w3gohz,d9rjek) = bmb(d9rjek,d9rjek)
23190 continue
23188 continue
      call dp2zwv(bmb, w8xfic, bgu6fw, rjcq9o, lku8xq, nfiumb4, xhe4cg, 
     &zqve1l, vvl1li, w3gohz)
23172 continue
      return
      end
      subroutine bf7qci(lku8xq, uq9jtc, bgu6fw, bmb)
      implicit logical (a-z)
      integer lku8xq
      double precision uq9jtc, bgu6fw(lku8xq,lku8xq), bmb(lku8xq,lku8xq)
      integer d9rjek, nd6mep
      do 23192 d9rjek=1,lku8xq 
      do 23194 nd6mep=1,lku8xq 
      bgu6fw(d9rjek,nd6mep) = bgu6fw(d9rjek,nd6mep) * uq9jtc
23194 continue
23192 continue
      do 23196 d9rjek=1,lku8xq 
      do 23198 nd6mep=1,lku8xq 
      bmb(nd6mep,d9rjek) = bmb(nd6mep,d9rjek) + bgu6fw(nd6mep,d9rjek)
23198 continue
23196 continue
      return
      end
      subroutine vsel(s, t, lku8xq, nk, a51l0o, minv, bgu6fw)
      implicit logical (a-z)
      integer s, t, lku8xq, nk, a51l0o
      double precision minv(a51l0o,nk*lku8xq), bgu6fw(lku8xq,lku8xq)
      integer w3gohz, d9rjek, y9eilo, pazyk8
      do 23200 w3gohz=1,lku8xq 
      do 23202 d9rjek=1,lku8xq 
      bgu6fw(w3gohz,d9rjek) = 0.0d0
23202 continue
23200 continue
      if(.not.(s .ne. t))goto 23204
      do 23206 w3gohz=1,lku8xq 
      y9eilo = (s-1)*lku8xq + w3gohz
      do 23208 d9rjek=1,lku8xq 
      pazyk8 = (t-1)*lku8xq + d9rjek
      bgu6fw(w3gohz,d9rjek) = minv(a51l0o-(pazyk8-y9eilo), pazyk8)
23208 continue
23206 continue
      goto 23205
23204 continue
      do 23210 w3gohz=1,lku8xq 
      y9eilo = (s-1)*lku8xq + w3gohz
      do 23212 d9rjek=w3gohz,lku8xq 
      pazyk8 = (t-1)*lku8xq + d9rjek
      bgu6fw(w3gohz,d9rjek) = minv(a51l0o-(pazyk8-y9eilo), pazyk8)
23212 continue
23210 continue
      do 23214 w3gohz=1,lku8xq 
      do 23216 d9rjek=w3gohz+1,lku8xq 
      bgu6fw(d9rjek,w3gohz) = bgu6fw(w3gohz,d9rjek)
23216 continue
23214 continue
23205 continue
      return
      end
      subroutine dp2zwv(bmb, w8xfic, bgu6fw, rjcq9o, lku8xq, nfiumb4, 
     &xhe4cg, zqve1l, vvl1li, p1rifj)
      implicit logical (a-z)
      integer lku8xq, nfiumb4, xhe4cg, zqve1l(1), vvl1li(1), p1rifj
      double precision bmb(lku8xq,lku8xq), w8xfic(nfiumb4,xhe4cg), 
     &bgu6fw(lku8xq,lku8xq), rjcq9o(nfiumb4,lku8xq)
      double precision qnk4zf, temp
      integer d9rjek, nd6mep, xi1mqb, i1nkrb
      do 23218 i1nkrb=1,lku8xq 
      do 23220 d9rjek=1,lku8xq 
      do 23222 nd6mep=1,lku8xq 
      bgu6fw(nd6mep,d9rjek) = 0.0d0
23222 continue
23220 continue
      do 23224 xi1mqb=1,xhe4cg 
      temp = w8xfic(p1rifj,xi1mqb)
      bgu6fw(zqve1l(xi1mqb),vvl1li(xi1mqb)) = temp
      bgu6fw(vvl1li(xi1mqb),zqve1l(xi1mqb)) = temp
23224 continue
      qnk4zf = 0.0d0
      do 23226 d9rjek=1,lku8xq 
      qnk4zf = qnk4zf + bmb(i1nkrb,d9rjek) * bgu6fw(d9rjek,i1nkrb)
23226 continue
      rjcq9o(p1rifj,i1nkrb) = qnk4zf
23218 continue
      return
      end
      subroutine gyzcj5(n5fkml, jrxg6l, d, uu, lku8xq, nfiumb4)
      implicit logical (a-z)
      integer lku8xq, nfiumb4
      double precision n5fkml(lku8xq+1,nfiumb4), jrxg6l(lku8xq+1,
     &nfiumb4), d(nfiumb4), uu(lku8xq+1,lku8xq+1)
      integer w3gohz, nd6mep, dibm1x, p4gdax, c3qxjo, j0izmn, myx3od
      n5fkml(lku8xq+1,nfiumb4) = 1.0d0 / d(nfiumb4)
      j0izmn = lku8xq+1
      c3qxjo = nfiumb4+1 - j0izmn
      do 23228 myx3od=c3qxjo,nfiumb4 
      do 23230 w3gohz=1,j0izmn 
      uu(w3gohz, myx3od-c3qxjo+1) = jrxg6l(w3gohz, myx3od)
23230 continue
23228 continue
      w3gohz = nfiumb4-1 
23232 if(.not.(w3gohz.ge.1))goto 23234
      if(.not.(lku8xq .lt. nfiumb4-w3gohz))goto 23235
      p4gdax = lku8xq
      goto 23236
23235 continue
      p4gdax = nfiumb4-w3gohz
23236 continue
      dibm1x=1
23237 if(.not.(dibm1x.le.p4gdax))goto 23239
      n5fkml(-dibm1x+lku8xq+1,w3gohz+dibm1x) = 0.0d0
      nd6mep=1
23240 if(.not.(nd6mep.le.dibm1x))goto 23242
      n5fkml(-dibm1x+lku8xq+1,w3gohz+dibm1x) = n5fkml(-dibm1x+lku8xq+1,
     &w3gohz+dibm1x) - uu(-nd6mep+lku8xq+1,w3gohz+nd6mep -c3qxjo+1) * 
     &n5fkml(nd6mep-dibm1x+lku8xq+1,w3gohz+dibm1x)
       nd6mep=nd6mep+1
      goto 23240
23242 continue
23243 if(.not.(nd6mep.le.p4gdax))goto 23245
      n5fkml(-dibm1x+lku8xq+1,w3gohz+dibm1x) = n5fkml(-dibm1x+lku8xq+1,
     &w3gohz+dibm1x) - uu(-nd6mep+lku8xq+1,w3gohz+nd6mep -c3qxjo+1) * 
     &n5fkml(dibm1x-nd6mep+lku8xq+1,w3gohz+nd6mep)
       nd6mep=nd6mep+1
      goto 23243
23245 continue
       dibm1x=dibm1x+1
      goto 23237
23239 continue
      n5fkml(lku8xq+1,w3gohz) = 1.0d0 / d(w3gohz)
      dibm1x = 1
23246 if(.not.(dibm1x.le.p4gdax))goto 23248
      n5fkml(lku8xq+1,w3gohz) = n5fkml(lku8xq+1,w3gohz) - uu(-dibm1x+
     &lku8xq+1,w3gohz+dibm1x -c3qxjo+1) * n5fkml(-dibm1x+lku8xq+1,
     &w3gohz+dibm1x)
       dibm1x=dibm1x+1
      goto 23246
23248 continue
      if(.not.(w3gohz .eq. c3qxjo))goto 23249
      c3qxjo = c3qxjo-1
      if(.not.(c3qxjo .lt. 1))goto 23251
      c3qxjo = 1
      goto 23252
23251 continue
      myx3od=j0izmn-1
23253 if(.not.(myx3od.ge.1))goto 23255
      nd6mep=1
23256 if(.not.(nd6mep.le.j0izmn))goto 23258
      uu(nd6mep,myx3od+1) = uu(nd6mep,myx3od)
       nd6mep=nd6mep+1
      goto 23256
23258 continue
       myx3od=myx3od-1
      goto 23253
23255 continue
      nd6mep=1
23259 if(.not.(nd6mep.le.j0izmn))goto 23261
      uu(nd6mep,1) = jrxg6l(nd6mep,c3qxjo)
       nd6mep=nd6mep+1
      goto 23259
23261 continue
23252 continue
23249 continue
       w3gohz = w3gohz-1
      goto 23232
23234 continue
      return
      end
      subroutine ntju9b(bz4guf,jmwo0z,w8xfic, nfiumb4,lku8xq,ynk9ah,
     &uxs1iq, l6xrjt,dof,smo,zo5jyl, s0, vy5hmo,yin,lj4dph,win, ykdc2t,
     &phqco4, xhe4cg, la5dcf, a51l0o, fjg0qv, y6jcvk, vb81l0, j1l0o1, 
     &qc7zyb, jko0o1,zqve1l,vvl1li, bmb, rjcq9o, zxao0o, wj5shg,du8jbv,
     &i83h1, ifo4ew, lq8reh, i0qvzl, jq8lra, kn7qya, vfd2pw, blq5vu, 
     &dfsom3)
      implicit logical (a-z)
      integer nfiumb4,lku8xq,ynk9ah(nfiumb4),uxs1iq, xhe4cg, la5dcf, 
     &a51l0o, fjg0qv, y6jcvk, jko0o1,zqve1l(1),vvl1li(1), wj5shg, 
     &du8jbv, i83h1(du8jbv*2)
      double precision bz4guf(nfiumb4), jmwo0z(nfiumb4,lku8xq), w8xfic(
     &nfiumb4,xhe4cg), l6xrjt(du8jbv), dof(du8jbv), smo(nfiumb4,du8jbv),
     & zo5jyl(nfiumb4,du8jbv)
      double precision s0(2*du8jbv, 2*du8jbv,2)
      double precision ykdc2t(1), phqco4(1), vb81l0(1), j1l0o1(1), 
     &qc7zyb(jko0o1+4)
      double precision vy5hmo(uxs1iq), yin(uxs1iq,lku8xq), lj4dph(
     &uxs1iq,lku8xq), win(uxs1iq,1), bmb(1), rjcq9o(uxs1iq,du8jbv), 
     &zxao0o(lku8xq,lku8xq,16), ifo4ew(lku8xq,du8jbv)
      double precision lq8reh(2*du8jbv), i0qvzl(2*du8jbv), jq8lra(
     &uxs1iq,du8jbv), kn7qya(du8jbv,uxs1iq), vfd2pw(du8jbv,uxs1iq), 
     &blq5vu(uxs1iq*du8jbv), dfsom3(1)
      integer ybfr6z
      integer w3gohz, d9rjek, nd6mep, c4bdmu, o9ljyn, tvyd2b, zx1610, 
     &c4uxow
      integer uxzze7
      double precision kogeb2, tap0km, t7sbea
      uxzze7 = 1
      if(.not.(wj5shg .eq. 1))goto 23262
      ybfr6z = xhe4cg
      goto 23263
23262 continue
      ybfr6z = du8jbv*(du8jbv+1)/2
23263 continue
      call qh4ulb(zqve1l, vvl1li, lku8xq)
      call vsuff9(nfiumb4,uxs1iq,ynk9ah, bz4guf,jmwo0z,w8xfic, vy5hmo,
     &yin,win,dfsom3,lj4dph, lku8xq, xhe4cg, ybfr6z, zqve1l, vvl1li, 
     &zxao0o, zxao0o(1,1,3), ifo4ew, du8jbv, wj5shg, uxzze7, c4uxow)
      if(.not.(c4uxow .ne. 1))goto 23264
      return
23264 continue
      kogeb2 = vy5hmo(1)
      tap0km = vy5hmo(uxs1iq)-vy5hmo(1)
      do 23266 w3gohz=1,uxs1iq 
      vy5hmo(w3gohz) = (vy5hmo(w3gohz)-kogeb2)/tap0km
23266 continue
      a51l0o = 4*du8jbv
      la5dcf = 0
      do 23268 d9rjek=1,du8jbv 
      if(.not.(l6xrjt(d9rjek) .eq. 0.0d0))goto 23270
      dof(d9rjek) = dof(d9rjek) + 1.0d0
23270 continue
23268 continue
      call qh4ulb(zqve1l, vvl1li, du8jbv)
      call vsplin(vy5hmo,lj4dph,win,uxs1iq,qc7zyb, jko0o1,a51l0o,du8jbv,
     &ybfr6z, zqve1l,vvl1li, zxao0o, l6xrjt, fjg0qv, jq8lra, j1l0o1, 
     &phqco4(1), phqco4(1+jko0o1*du8jbv*a51l0o), vb81l0, zo5jyl, y6jcvk,
     & bmb, rjcq9o, dof, ykdc2t, la5dcf, nfiumb4)
      do 23272 d9rjek=1,du8jbv 
      dof(d9rjek) = -1.0d0
      do 23274 w3gohz=1,uxs1iq 
      dof(d9rjek)=dof(d9rjek)+rjcq9o(w3gohz,d9rjek)
23274 continue
23272 continue
      if(.not.(du8jbv .ge. 1))goto 23276
      t7sbea = 1.0d-7
      c4bdmu = uxs1iq*du8jbv
      o9ljyn = 2*du8jbv
      tvyd2b = 101
      fjg0qv = 1
      call kgevo5(vy5hmo, phqco4, uxs1iq, du8jbv)
      call qh4ulb(zqve1l, vvl1li, du8jbv)
      call mux17f(dfsom3, phqco4, du8jbv, o9ljyn, uxs1iq, zxao0o(1,1,1),
     & zxao0o(1,1,2), zqve1l, vvl1li, ybfr6z, c4bdmu)
      do 23278 nd6mep=1,o9ljyn 
      i83h1(nd6mep) = nd6mep
23278 continue
      call dhkt9w(phqco4,c4bdmu,c4bdmu,o9ljyn,i0qvzl,i83h1,ykdc2t,
     &zx1610,t7sbea)
      call qh4ulb(zqve1l, vvl1li, du8jbv)
      call mux22f(dfsom3,jq8lra,kn7qya,ybfr6z,zqve1l,vvl1li,uxs1iq,
     &du8jbv,zxao0o)
      call vdqrsl(phqco4,c4bdmu,c4bdmu,zx1610,i0qvzl,kn7qya,ykdc2t(1),
     &blq5vu,lq8reh, ykdc2t(1),vfd2pw,tvyd2b,fjg0qv)
      call vbksf(dfsom3,vfd2pw,du8jbv,uxs1iq,zxao0o,zqve1l,vvl1li,
     &ybfr6z)
      if(.not.(y6jcvk .ne. 0))goto 23280
      call vrinvf9(phqco4, c4bdmu, o9ljyn, c4uxow, s0(1,1,1), s0(1,1,2))
      if(.not.(c4uxow .ne. 1))goto 23282
      return
23282 continue
      do 23284 d9rjek=1,du8jbv 
      do 23286 w3gohz=1,uxs1iq 
      zo5jyl(w3gohz,d9rjek) = zo5jyl(w3gohz,d9rjek) - s0(d9rjek,d9rjek,
     &1) - vy5hmo(w3gohz) * (2.0d0 * s0(d9rjek,d9rjek+du8jbv,1) + 
     &vy5hmo(w3gohz) * s0(d9rjek+du8jbv,d9rjek+du8jbv,1))
23286 continue
23284 continue
23280 continue
      goto 23277
23276 continue
      call rpfnk6(uxs1iq, vy5hmo, win, jq8lra, vfd2pw, zo5jyl, y6jcvk)
23277 continue
      do 23288 w3gohz=1,uxs1iq 
      do 23290 d9rjek=1,du8jbv 
      jq8lra(w3gohz,d9rjek) = jq8lra(w3gohz,d9rjek) - vfd2pw(d9rjek,
     &w3gohz)
23290 continue
23288 continue
      do 23292 d9rjek=1,du8jbv 
      call uwye7d(nfiumb4, uxs1iq, ynk9ah, jq8lra(1,d9rjek), smo(1,
     &d9rjek))
23292 continue
      return
      end
      subroutine vbfa(d8gwha,n,lku8xq,gqai81, p3vlea,jmwo0z,w8xfic,
     &l6xrjt,dof, ynk9ah,uxs1iq,vliac4, vfd2pw,sazp9g,go0l1q,s0, lq8reh,
     &zo5jyl,h4fgoy, ioqzvb,i0qvzl, i83h1, xbig, jrxg6l, ifo4ew, ozuw3p,
     & hwi2tb, nbd5rl, wj5shg, ykdc2t, wk2, zxao0o, phqco4, vb81l0, bmb,
     & rjcq9o, mwk, t5vlzq, j1l0o1, qc7zyb, das4bx, vlni8d, jko0o1, 
     &mnh3up, fg3pxq)
      implicit logical (a-z)
      integer d8gwha, n, lku8xq, gqai81(15), ynk9ah(1),uxs1iq(1),vliac4(
     &1), i83h1(1)
      integer ozuw3p(1), hwi2tb(1), nbd5rl(1), wj5shg(1), vlni8d(1), 
     &jko0o1(1), mnh3up(1), fg3pxq(1)
      double precision p3vlea(1),jmwo0z(1),w8xfic(1),l6xrjt(1),dof(1), 
     &vfd2pw(1),sazp9g(1), go0l1q(1), s0(lku8xq), lq8reh(1),zo5jyl(1),
     &h4fgoy, ioqzvb(1),i0qvzl(1)
      double precision xbig(1), jrxg6l(1), ifo4ew(1), ykdc2t(1), wk2(n,
     &lku8xq,3), zxao0o(lku8xq,lku8xq,16), phqco4(1), vb81l0(1), bmb(1),
     & rjcq9o(1), mwk(1), t5vlzq(1), j1l0o1(1), qc7zyb(1), das4bx
      integer p,q,y6jcvk,nucgi1r,no2fik, c4bdmu, o9ljyn, tiav4e, xhe4cg,
     & zkjqhi, la5dcf,a51l0o
      integer ucgi1r
      integer sehz7y
      integer w3gohz, j0qwtz, zx1610
      j0qwtz = lku8xq*(lku8xq+1)/2
      p=gqai81(2)
      q=gqai81(3)
      y6jcvk= 0
      if(.not.(gqai81(4) .eq. 1))goto 23294
      y6jcvk = 1
23294 continue
      no2fik=gqai81(6)
      zx1610=gqai81(7)
      c4bdmu=gqai81(9)
      o9ljyn=gqai81(10)
      tiav4e=gqai81(11)
      xhe4cg=gqai81(12)
      zkjqhi=gqai81(13)
      la5dcf = 0
      a51l0o=gqai81(15)
      sehz7y = 1
      if(.not.(tiav4e .gt. 0))goto 23296
      do 23298 w3gohz=1,tiav4e 
      ykdc2t(w3gohz) = dof(w3gohz)
      ykdc2t(w3gohz+tiav4e) = l6xrjt(w3gohz)
      ykdc2t(w3gohz+2*tiav4e) = dof(w3gohz)
23298 continue
23296 continue
      ucgi1r = 0
23300 if(.not.(sehz7y .ne. 0))goto 23301
      ucgi1r = ucgi1r+1
      if(.not.(ucgi1r .gt. 1))goto 23302
      if(.not.(tiav4e .gt. 0))goto 23304
      do 23306 w3gohz=1,tiav4e 
      if(.not.(ykdc2t(w3gohz+tiav4e).eq.0.0d0 .and.(dabs(ykdc2t(w3gohz+
     &2*tiav4e)-dof(w3gohz))/dof(w3gohz).gt.0.05d0)))goto 23308
      ykdc2t(w3gohz+2*tiav4e) = dof(w3gohz)
      dof(w3gohz)=ykdc2t(w3gohz)
      l6xrjt(w3gohz)=0.0d0
      goto 23309
23308 continue
      ykdc2t(w3gohz+2*tiav4e) = dof(w3gohz)
23309 continue
23306 continue
23304 continue
23302 continue
      call xqasw0(d8gwha,n,lku8xq, p3vlea,jmwo0z,w8xfic,l6xrjt,dof, 
     &ynk9ah,uxs1iq,vliac4, vfd2pw,sazp9g,go0l1q,s0, lq8reh,zo5jyl,
     &h4fgoy, ioqzvb,i0qvzl, zx1610,i83h1, xbig, jrxg6l, ifo4ew, ozuw3p,
     & hwi2tb, nbd5rl(1), nbd5rl(1 + j0qwtz), wj5shg, ykdc2t(1+3*tiav4e)
     &, zxao0o, phqco4, vb81l0, bmb, rjcq9o, mwk, t5vlzq, j1l0o1, 
     &qc7zyb, das4bx, vlni8d, jko0o1, mnh3up, fg3pxq, p,q,y6jcvk,
     &nucgi1r,no2fik, wk2(1,1,1), wk2(1,1,2), wk2(1,1,3), c4bdmu, 
     &o9ljyn, tiav4e, xhe4cg, zkjqhi, la5dcf, a51l0o)
      if(.not.(d8gwha .ne. 0))goto 23310
      call vcall2(sehz7y,w,y,go0l1q,lq8reh,jrxg6l)
      goto 23311
23310 continue
      sehz7y = 0
23311 continue
      if(.not.(sehz7y .ne. 0))goto 23312
      zx1610=0
23312 continue
      goto 23300
23301 continue
      gqai81(7) = zx1610
      gqai81(5) = nucgi1r
      gqai81(14) = la5dcf
      return
      end
      subroutine xqasw0(d8gwha,nfiumb4,lku8xq, p3vlea,jmwo0z,w8xfic,
     &l6xrjt,dof, ynk9ah,uxs1iq,vliac4, vfd2pw,sazp9g,go0l1q,s0, lq8reh,
     &zo5jyl,h4fgoy, ioqzvb,i0qvzl, zx1610,i83h1, xbig, jrxg6l, ifo4ew, 
     &ozuw3p, hwi2tb, zqve1l, vvl1li, wj5shg, ykdc2t, zxao0o, phqco4, 
     &vb81l0, bmb, rjcq9o, mwk, t5vlzq, j1l0o1, qc7zyb, das4bx, vlni8d, 
     &jko0o1, mnh3up, fg3pxq, p, q, y6jcvk, nucgi1r, no2fik, hr83e, 
     &x7aort, wk2, c4bdmu, o9ljyn, tiav4e, xhe4cg, zkjqhi, la5dcf, 
     &a51l0o)
      implicit logical (a-z)
      integer zx1610
      integer vvl1li(1), zqve1l(1)
      integer p, q, y6jcvk, nucgi1r, no2fik, c4bdmu, o9ljyn, tiav4e, 
     &xhe4cg, zkjqhi, la5dcf, a51l0o
      integer d8gwha, nfiumb4, lku8xq, ynk9ah(nfiumb4,q),uxs1iq(q),
     &vliac4(q), i83h1(o9ljyn)
      integer ozuw3p(q), hwi2tb(q), wj5shg(q), vlni8d(q+1), jko0o1(q), 
     &mnh3up(1), fg3pxq(q+1)
      double precision p3vlea(nfiumb4,p), jmwo0z(nfiumb4,lku8xq), 
     &w8xfic(nfiumb4,xhe4cg), l6xrjt(tiav4e), dof(tiav4e)
      double precision vfd2pw(lku8xq,nfiumb4), sazp9g(nfiumb4,tiav4e), 
     &go0l1q(lku8xq,nfiumb4), s0(lku8xq), lq8reh(o9ljyn), zo5jyl(
     &nfiumb4,tiav4e), h4fgoy, ioqzvb(c4bdmu,o9ljyn), i0qvzl(o9ljyn)
      double precision xbig(c4bdmu,o9ljyn), jrxg6l(zkjqhi,nfiumb4), 
     &ifo4ew(lku8xq,tiav4e), ykdc2t(1), wk2(nfiumb4,lku8xq), zxao0o(
     &lku8xq,lku8xq,16), phqco4(1), vb81l0(1), bmb(1), rjcq9o(1), mwk(1)
     &, t5vlzq(1), j1l0o1(1), qc7zyb(1), das4bx
      double precision hr83e(nfiumb4,lku8xq), x7aort(nfiumb4,lku8xq)
      integer tvyd2b,fjg0qv,rbwx6v
      integer w3gohz, d9rjek, nd6mep, jbyv3q
      double precision gwu72m, jcp1ti,dyb3po, njdgw8, gcjn3k,t7sbea
      t7sbea = 1.0d-7
      tvyd2b = 101
      fjg0qv = 1
      if(.not.(q .eq. 0))goto 23314
      no2fik = 1
23314 continue
      if(.not.(d8gwha .ne. 0))goto 23316
      do 23318 d9rjek=1,o9ljyn 
      do 23320 w3gohz=1,c4bdmu 
      ioqzvb(w3gohz,d9rjek)=xbig(w3gohz,d9rjek)
23320 continue
23318 continue
23316 continue
      if(.not.(zx1610.eq.0))goto 23322
      call qh4ulb(zqve1l,vvl1li,lku8xq)
      call mux17f(jrxg6l, ioqzvb, lku8xq, o9ljyn, nfiumb4, zxao0o(1,1,1)
     &, zxao0o(1,1,2), zqve1l, vvl1li, zkjqhi, c4bdmu)
      do 23324 nd6mep=1,o9ljyn 
      i83h1(nd6mep) = nd6mep
23324 continue
      call dhkt9w(ioqzvb,c4bdmu,c4bdmu,o9ljyn,i0qvzl,i83h1,t5vlzq,
     &zx1610,t7sbea)
23322 continue
      do 23326 d9rjek=1,lku8xq 
      do 23328 w3gohz=1,nfiumb4 
      go0l1q(d9rjek,w3gohz)=0.0d0
23328 continue
      if(.not.(q .gt. 0))goto 23330
      do 23332 nd6mep=1,q 
      if(.not.(wj5shg(nd6mep).eq.1))goto 23334
      do 23336 w3gohz=1,nfiumb4 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) + sazp9g(w3gohz,
     &hwi2tb(nd6mep)+d9rjek-1)
23336 continue
      goto 23335
23334 continue
      do 23338 jbyv3q=1,ozuw3p(nd6mep) 
      do 23340 w3gohz=1,nfiumb4 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) + ifo4ew(d9rjek,
     &hwi2tb(nd6mep)+jbyv3q-1) * sazp9g(w3gohz,hwi2tb(nd6mep)+jbyv3q-1)
23340 continue
23338 continue
23335 continue
23332 continue
23330 continue
23326 continue
      nucgi1r = 0
      dyb3po = 1.0d0
23342 if(.not.((dyb3po .gt. h4fgoy ) .and. (nucgi1r .lt. no2fik)))goto 2
     &3343
      nucgi1r = nucgi1r + 1
      njdgw8 = 0.0d0
      do 23344 d9rjek=1,lku8xq 
      do 23346 w3gohz=1,nfiumb4 
      hr83e(w3gohz,d9rjek)=jmwo0z(w3gohz,d9rjek)-go0l1q(d9rjek,w3gohz)
23346 continue
23344 continue
      call qh4ulb(zqve1l,vvl1li,lku8xq)
      call mux22f(jrxg6l,hr83e, t5vlzq, zkjqhi,zqve1l,vvl1li,nfiumb4,
     &lku8xq,zxao0o)
      call vdqrsl(ioqzvb,c4bdmu,c4bdmu,zx1610,i0qvzl, t5vlzq, wk2,wk2, 
     &lq8reh, wk2,vfd2pw,tvyd2b,fjg0qv)
      das4bx=0.0d0
      do 23348 w3gohz=1,nfiumb4 
      do 23350 d9rjek=1,lku8xq 
      gwu72m = t5vlzq((w3gohz-1)*lku8xq+d9rjek) - vfd2pw(d9rjek,w3gohz)
      das4bx = das4bx + gwu72m * gwu72m
23350 continue
23348 continue
      call vbksf(jrxg6l,vfd2pw,lku8xq,nfiumb4,zxao0o,zqve1l,vvl1li,
     &zkjqhi)
      if(.not.(q .gt. 0))goto 23352
      do 23354 nd6mep=1,q 
      do 23356 d9rjek=1,lku8xq 
      if(.not.(wj5shg(nd6mep).eq.1))goto 23358
      do 23360 w3gohz=1,nfiumb4 
      x7aort(w3gohz,d9rjek)=sazp9g(w3gohz,hwi2tb(nd6mep)+d9rjek-1)
      hr83e(w3gohz,d9rjek) = jmwo0z(w3gohz,d9rjek) - vfd2pw(d9rjek,
     &w3gohz) - go0l1q(d9rjek,w3gohz) + x7aort(w3gohz,d9rjek)
23360 continue
      goto 23359
23358 continue
      do 23362 w3gohz=1,nfiumb4 
      x7aort(w3gohz,d9rjek)=0.0d0
      do 23364 jbyv3q=1,ozuw3p(nd6mep) 
      x7aort(w3gohz,d9rjek)=x7aort(w3gohz,d9rjek) + ifo4ew(d9rjek,
     &hwi2tb(nd6mep)+jbyv3q-1) * sazp9g(w3gohz,hwi2tb(nd6mep)+jbyv3q-1)
23364 continue
      hr83e(w3gohz,d9rjek) = jmwo0z(w3gohz,d9rjek) - vfd2pw(d9rjek,
     &w3gohz) - go0l1q(d9rjek,w3gohz) + x7aort(w3gohz,d9rjek)
23362 continue
23359 continue
23356 continue
      rbwx6v = uxs1iq(nd6mep)
      call ntju9b(p3vlea(1,vliac4(nd6mep)),hr83e,w8xfic, nfiumb4,lku8xq,
     &ynk9ah(1,nd6mep),rbwx6v, l6xrjt(hwi2tb(nd6mep)), dof(hwi2tb(
     &nd6mep)), sazp9g(1,hwi2tb(nd6mep)), zo5jyl(1,hwi2tb(nd6mep)), s0, 
     &mwk(1), mwk(1+rbwx6v), mwk(1+rbwx6v*(lku8xq+1)), mwk(1+rbwx6v*(2*
     &lku8xq+1)), ykdc2t, phqco4, xhe4cg, la5dcf, a51l0o, fjg0qv, 
     &y6jcvk, vb81l0, j1l0o1(vlni8d(nd6mep)), qc7zyb(fg3pxq(nd6mep)), 
     &jko0o1(nd6mep),zqve1l, vvl1li, bmb, rjcq9o, zxao0o, wj5shg(nd6mep)
     &,ozuw3p(nd6mep),mnh3up, ifo4ew(1,hwi2tb(nd6mep)), t5vlzq(1), 
     &t5vlzq(1+2*ozuw3p(nd6mep)), t5vlzq(1+4*ozuw3p(nd6mep)), t5vlzq(1+(
     &4+rbwx6v)*ozuw3p(nd6mep)), t5vlzq(1+(4+2*rbwx6v)*ozuw3p(nd6mep)), 
     &t5vlzq(1+(4+3*rbwx6v)*ozuw3p(nd6mep)), t5vlzq(1+(4+4*rbwx6v)*
     &ozuw3p(nd6mep)))
      do 23366 d9rjek=1,lku8xq 
      if(.not.(wj5shg(nd6mep).eq.1))goto 23368
      do 23370 w3gohz=1,nfiumb4 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) + sazp9g(w3gohz,
     &hwi2tb(nd6mep)+d9rjek-1)
23370 continue
      goto 23369
23368 continue
      do 23372 jbyv3q=1,ozuw3p(nd6mep) 
      do 23374 w3gohz=1,nfiumb4 
      go0l1q(d9rjek,w3gohz)=go0l1q(d9rjek,w3gohz) + ifo4ew(d9rjek,
     &hwi2tb(nd6mep)+jbyv3q-1) * sazp9g(w3gohz,hwi2tb(nd6mep)+jbyv3q-1)
23374 continue
23372 continue
23369 continue
      do 23376 w3gohz=1,nfiumb4 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) - x7aort(w3gohz,
     &d9rjek)
23376 continue
23366 continue
      do 23378 d9rjek=1,lku8xq 
      if(.not.(wj5shg(nd6mep) .eq. 1))goto 23380
      njdgw8 = njdgw8 + jcp1ti(nfiumb4,x7aort(1,d9rjek),sazp9g(1,hwi2tb(
     &nd6mep)+d9rjek-1), w8xfic(1,d9rjek))
      goto 23381
23380 continue
      do 23382 w3gohz=1,nfiumb4 
      t5vlzq(w3gohz) = 0.0d0
      do 23384 jbyv3q=1,ozuw3p(nd6mep) 
      t5vlzq(w3gohz) = t5vlzq(w3gohz) + ifo4ew(d9rjek,hwi2tb(nd6mep)+
     &jbyv3q-1) * sazp9g(w3gohz,hwi2tb(nd6mep)+jbyv3q-1)
23384 continue
23382 continue
      njdgw8 = njdgw8 + jcp1ti(nfiumb4, x7aort(1,d9rjek), t5vlzq, 
     &w8xfic(1,d9rjek))
23381 continue
23378 continue
      do 23386 d9rjek=1,lku8xq 
      do 23388 w3gohz=1,nfiumb4 
      hr83e(w3gohz,d9rjek)=jmwo0z(w3gohz,d9rjek)-go0l1q(d9rjek,w3gohz)
23388 continue
23386 continue
      call qh4ulb(zqve1l,vvl1li,lku8xq)
      call mux22f(jrxg6l,hr83e, t5vlzq, zkjqhi,zqve1l,vvl1li,nfiumb4,
     &lku8xq,zxao0o)
      call vdqrsl(ioqzvb,c4bdmu,c4bdmu,zx1610,i0qvzl, t5vlzq, wk2,wk2, 
     &lq8reh, wk2,vfd2pw,tvyd2b,fjg0qv)
      call vbksf(jrxg6l,vfd2pw,lku8xq,nfiumb4,zxao0o,zqve1l,vvl1li,
     &zkjqhi)
23354 continue
23352 continue
      if(.not.(q .gt. 0))goto 23390
      gcjn3k=0.0d0
      do 23392 d9rjek=1,lku8xq 
      do 23394 w3gohz=1,nfiumb4 
      gcjn3k = gcjn3k + w8xfic(w3gohz,d9rjek) * go0l1q(d9rjek,w3gohz)**
     &2
23394 continue
23392 continue
      if(.not.(gcjn3k .gt. 0.0d0))goto 23396
      dyb3po = dsqrt(njdgw8/gcjn3k)
      goto 23397
23396 continue
      dyb3po = 0.0d0
23397 continue
23390 continue
      if(.not.(nucgi1r .eq. 1))goto 23398
      dyb3po = 1.0d0
23398 continue
      goto 23342
23343 continue
      do 23400 d9rjek=1,o9ljyn 
      t5vlzq(d9rjek)=lq8reh(d9rjek)
23400 continue
      do 23402 d9rjek=1,o9ljyn 
      lq8reh(i83h1(d9rjek))=t5vlzq(d9rjek)
23402 continue
      do 23404 w3gohz=1,nfiumb4 
      do 23406 d9rjek=1,lku8xq 
      go0l1q(d9rjek,w3gohz) = go0l1q(d9rjek,w3gohz) + vfd2pw(d9rjek,
     &w3gohz)
23406 continue
23404 continue
      if(.not.((y6jcvk .ne. 0) .and. (q .gt. 0)))goto 23408
      do 23410 nd6mep=1,q 
      do 23412 jbyv3q=1,ozuw3p(nd6mep) 
      call uwye7d(nfiumb4,uxs1iq(nd6mep),ynk9ah(1,nd6mep), zo5jyl(1,
     &hwi2tb(nd6mep)+jbyv3q-1),x7aort)
      do 23414 w3gohz=1,nfiumb4 
      zo5jyl(w3gohz,hwi2tb(nd6mep)+jbyv3q-1) = x7aort(w3gohz,1)
23414 continue
23412 continue
23410 continue
23408 continue
      return
      end
      subroutine kgevo5(p3vlea, xout, nfiumb4, lku8xq)
      implicit logical (a-z)
      integer nfiumb4, lku8xq
      double precision p3vlea(nfiumb4), xout(1)
      integer w3gohz, d9rjek, nd6mep, xtiel4
      xtiel4=1
      do 23416 d9rjek=1,lku8xq 
      do 23418 w3gohz=1,nfiumb4 
      do 23420 nd6mep=1,lku8xq 
      if(.not.(d9rjek .eq. nd6mep))goto 23422
      xout(xtiel4) = 1.0d0
      goto 23423
23422 continue
      xout(xtiel4) = 0.0d0
23423 continue
      xtiel4=xtiel4+1
23420 continue
23418 continue
23416 continue
      do 23424 d9rjek=1,lku8xq 
      do 23426 w3gohz=1,nfiumb4 
      do 23428 nd6mep=1,lku8xq 
      if(.not.(d9rjek .eq. nd6mep))goto 23430
      xout(xtiel4) = p3vlea(w3gohz)
      goto 23431
23430 continue
      xout(xtiel4) = 0.0d0
23431 continue
      xtiel4=xtiel4+1
23428 continue
23426 continue
23424 continue
      return
      end
      double precision function jcp1ti(nfiumb4, yvec, go0l1q, wvec)
      integer nfiumb4
      double precision yvec(nfiumb4), go0l1q(nfiumb4), wvec(nfiumb4)
      integer w3gohz
      double precision wtot, risyv0, bgu6fw
      risyv0 = 0.0d0
      wtot = 0.0d0
      do 23432 w3gohz=1,nfiumb4 
      bgu6fw = yvec(w3gohz) - go0l1q(w3gohz)
      risyv0 = risyv0 + wvec(w3gohz)*bgu6fw*bgu6fw
      wtot = wtot + wvec(w3gohz)
23432 continue
      if(.not.(wtot .gt. 0.0d0))goto 23434
      jcp1ti=risyv0/wtot
      goto 23435
23434 continue
      jcp1ti=0.0d0
23435 continue
      return
      end
      subroutine usytl1(nfiumb4, yvec, wvec, ghry8z, wtot)
      implicit logical (a-z)
      integer nfiumb4
      double precision yvec(nfiumb4), wvec(nfiumb4), ghry8z, wtot
      double precision risyv0
      integer w3gohz
      wtot = 0.0d0
      risyv0 = 0.0d0
      do 23436 w3gohz=1,nfiumb4 
      risyv0 = risyv0 + yvec(w3gohz) * wvec(w3gohz)
      wtot = wtot + wvec(w3gohz)
23436 continue
      if(.not.(wtot .gt. 0.0d0))goto 23438
      ghry8z = risyv0 / wtot
      goto 23439
23438 continue
      ghry8z = 0.0d0
23439 continue
      return
      end
      subroutine rpfnk6(nfiumb4, x, w, yvec, vfd2pw, zo5jyl, y6jcvk)
      implicit logical (a-z)
      integer nfiumb4
      integer y6jcvk
      double precision x(nfiumb4), w(nfiumb4), yvec(nfiumb4), vfd2pw(
     &nfiumb4)
      double precision zo5jyl(nfiumb4,1)
      integer w3gohz
      double precision bz3pyo, qxy6aj, qnk4zf, vgh4cp, u7hbqo, agfy3b, 
     &qe3jcd, j0izmn, wtot
      call usytl1(nfiumb4,yvec,w,bz3pyo, wtot)
      call usytl1(nfiumb4,x,w,qxy6aj, wtot)
      vgh4cp = 0.0d0
      qnk4zf = 0.0d0
      do 23440 w3gohz=1,nfiumb4 
      j0izmn = x(w3gohz)-qxy6aj
      vgh4cp = vgh4cp + j0izmn * (yvec(w3gohz)-bz3pyo) * w(w3gohz)
      j0izmn = j0izmn * j0izmn
      qnk4zf = qnk4zf + j0izmn * w(w3gohz)
23440 continue
      u7hbqo = vgh4cp/qnk4zf
      agfy3b = bz3pyo - u7hbqo * qxy6aj
      do 23442 w3gohz=1,nfiumb4 
      vfd2pw(w3gohz) = agfy3b + u7hbqo * x(w3gohz)
23442 continue
      qe3jcd = agfy3b + u7hbqo * x(1)
      if(.not.(y6jcvk .ne. 0))goto 23444
      do 23446 w3gohz=1,nfiumb4 
      j0izmn = x(w3gohz)-qxy6aj
      if(.not.(w(w3gohz) .gt. 0.0d0))goto 23448
      zo5jyl(w3gohz,1) = zo5jyl(w3gohz,1) - 1.0d0/wtot - j0izmn * 
     &j0izmn / qnk4zf
      goto 23449
23448 continue
      zo5jyl(w3gohz,1) = 0.0d0
23449 continue
23446 continue
23444 continue
      return
      end
      subroutine uwye7d(nfiumb4, p, ynk9ah, qxy6aj, x)
      implicit logical (a-z)
      integer nfiumb4, p, ynk9ah(nfiumb4)
      double precision qxy6aj(p), x(nfiumb4)
      integer w3gohz
      do 23450 w3gohz=1,nfiumb4 
      x(w3gohz) = qxy6aj(ynk9ah(w3gohz))
23450 continue
      return
      end
      subroutine vknotl2(x, nfiumb4, knot, xl6qgm, q9wyop)
      implicit logical (a-z)
      integer nfiumb4, xl6qgm, q9wyop
      double precision x(nfiumb4), knot(nfiumb4)
      integer ndk, d9rjek
      if(.not.(q9wyop .eq. 0))goto 23452
      if(.not.(nfiumb4 .le. 40))goto 23454
      ndk = nfiumb4
      goto 23455
23454 continue
      ndk = 40 + dexp(0.25d0 * dlog(nfiumb4-40.0d0))
23455 continue
      goto 23453
23452 continue
      ndk = xl6qgm - 6
23453 continue
      xl6qgm = ndk + 6
      do 23456 d9rjek = 1,3 
      knot(d9rjek) = x(1) 
23456 continue
      do 23458 d9rjek = 1,ndk 
      knot(d9rjek+3) = x( 1 + (d9rjek-1)*(nfiumb4-1)/(ndk-1) ) 
23458 continue
      do 23460 d9rjek = 1,3 
      knot(ndk+3+d9rjek) = x(nfiumb4) 
23460 continue
      return
      end
      subroutine pknotl2(knot, nfiumb4, keep, fjo2dy)
      implicit logical (a-z)
      integer nfiumb4, keep(nfiumb4)
      double precision knot(nfiumb4), fjo2dy
      integer w3gohz, ilower
      do 23462 w3gohz=1,4 
      keep(w3gohz) = 1
23462 continue
      ilower = 4
      do 23464 w3gohz=5,(nfiumb4-4) 
      if(.not.((knot(w3gohz) - knot(ilower) .ge. fjo2dy) .and.(knot(
     &nfiumb4) - knot(w3gohz) .ge. fjo2dy)))goto 23466
      keep(w3gohz) = 1
      ilower = w3gohz
      goto 23467
23466 continue
      keep(w3gohz) = 0
23467 continue
23464 continue
      do 23468 w3gohz=(nfiumb4-3),nfiumb4 
      keep(w3gohz) = 1
23468 continue
      return
      end
      subroutine vglmf(xbig,c4bdmu,o9ljyn,d8gwha,nfiumb4, jmwo0z,lq8reh,
     &go0l1q,blq5vu, ioqzvb,i0qvzl, zx1610,i83h1, bgu6fw,zxao0o,jrxg6l, 
     &lku8xq,zkjqhi,xhe4cg, zqve1l, vvl1li, cpxbig, das4bx)
      implicit logical (a-z)
      integer c4bdmu,o9ljyn,d8gwha,nfiumb4, zx1610,i83h1(o9ljyn), 
     &lku8xq,zkjqhi,xhe4cg, zqve1l(1), vvl1li(1), cpxbig
      double precision xbig(c4bdmu,o9ljyn), jmwo0z(nfiumb4,lku8xq),
     &lq8reh(o9ljyn),go0l1q(lku8xq,nfiumb4),blq5vu(c4bdmu), ioqzvb(
     &c4bdmu,o9ljyn),i0qvzl(o9ljyn), bgu6fw(1), zxao0o(lku8xq,lku8xq,5),
     & jrxg6l(1), das4bx
      integer sehz7y
      call qh4ulb(zqve1l,vvl1li,lku8xq)
      sehz7y = 1
23470 if(.not.(sehz7y .ne. 0))goto 23471
      call vfit(lku8xq,c4bdmu,o9ljyn,nfiumb4, xbig,jmwo0z,lq8reh,go0l1q,
     &blq5vu, jrxg6l,ioqzvb,i0qvzl,das4bx, zx1610,i83h1, bgu6fw,zxao0o, 
     &xhe4cg,zkjqhi,zqve1l,vvl1li)
      if(.not.(d8gwha .ne. 0))goto 23472
      call vcall1(sehz7y,jmwo0z,go0l1q,lq8reh,jrxg6l,xbig,cpxbig)
      goto 23473
23472 continue
      sehz7y= 0
23473 continue
      if(.not.(sehz7y .ne. 0))goto 23474
      zx1610=0
23474 continue
      goto 23470
23471 continue
      return
      end
      subroutine vfit(lku8xq,c4bdmu,o9ljyn,nfiumb4, xbig,jmwo0z,lq8reh,
     &go0l1q,blq5vu, jrxg6l,ioqzvb,i0qvzl,das4bx, zx1610,i83h1, bgu6fw,
     &zxao0o, xhe4cg, zkjqhi, zqve1l, vvl1li)
      implicit logical (a-z)
      integer lku8xq, c4bdmu, o9ljyn, nfiumb4, zx1610, i83h1(o9ljyn), 
     &xhe4cg, zkjqhi, zqve1l(1), vvl1li(1)
      double precision xbig(c4bdmu,o9ljyn), jmwo0z(nfiumb4,lku8xq), 
     &lq8reh(o9ljyn), go0l1q(lku8xq,nfiumb4), blq5vu(c4bdmu), jrxg6l(
     &zkjqhi,nfiumb4), ioqzvb(c4bdmu,o9ljyn), i0qvzl(o9ljyn), das4bx, 
     &bgu6fw(c4bdmu), zxao0o(lku8xq,lku8xq,5)
      integer w3gohz, d9rjek, nd6mep, xi1mqb, hv3wja, tvyd2b, fjg0qv
      double precision gwu72m, dyb3po, t7sbea
      t7sbea=1.0d-7
      tvyd2b=101
      fjg0qv=1
      dyb3po=1.0d0
      if(.not.(zx1610 .eq. 0))goto 23476
      do 23478 d9rjek=1,o9ljyn 
      do 23480 w3gohz=1,c4bdmu 
      ioqzvb(w3gohz,d9rjek) = xbig(w3gohz,d9rjek)
23480 continue
23478 continue
      do 23482 nd6mep=1,o9ljyn 
      i83h1(nd6mep) = nd6mep
23482 continue
      call mux17f(jrxg6l, ioqzvb, lku8xq, o9ljyn, nfiumb4, zxao0o(1,1,1)
     &, zxao0o(1,1,2), zqve1l, vvl1li, zkjqhi, c4bdmu)
      call dhkt9w(ioqzvb,c4bdmu,c4bdmu,o9ljyn,i0qvzl,i83h1,bgu6fw,
     &zx1610,t7sbea)
23476 continue
      call mux22f(jrxg6l,jmwo0z,bgu6fw,zkjqhi,zqve1l,vvl1li,nfiumb4,
     &lku8xq,zxao0o)
      nd6mep=1
      do 23484 d9rjek=1,lku8xq 
      do 23486 w3gohz=1,nfiumb4 
      jmwo0z(w3gohz,d9rjek)=bgu6fw(nd6mep)
      nd6mep=nd6mep+1
23486 continue
23484 continue
      call vdqrsl(ioqzvb,c4bdmu,c4bdmu,zx1610,i0qvzl,jmwo0z,bgu6fw(1),
     &blq5vu,lq8reh, bgu6fw(1),go0l1q,tvyd2b,fjg0qv)
      das4bx=0.0d0
      xi1mqb=0
      hv3wja=1
      do 23488 w3gohz=1,nfiumb4 
      do 23490 d9rjek=1,lku8xq 
      xi1mqb = xi1mqb + 1
      if(.not.(xi1mqb .gt. nfiumb4))goto 23492
      xi1mqb = 1
      hv3wja = hv3wja + 1
23492 continue
      gwu72m = jmwo0z(xi1mqb,hv3wja) - go0l1q(d9rjek,w3gohz)
      das4bx = das4bx + gwu72m * gwu72m
23490 continue
23488 continue
      call vbksf(jrxg6l,go0l1q,lku8xq,nfiumb4,zxao0o,zqve1l,vvl1li,
     &xhe4cg)
      do 23494 d9rjek=1,o9ljyn 
      bgu6fw(d9rjek) = lq8reh(d9rjek)
23494 continue
      do 23496 d9rjek=1,o9ljyn 
      lq8reh(i83h1(d9rjek)) = bgu6fw(d9rjek)
23496 continue
      return
      end
