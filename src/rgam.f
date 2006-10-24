      subroutine nvhb7f(egoxa3,atqh0o,xs,ys,ws, nfiumb4,nk, knot,coef,
     &sz,rjcq9o, n9peut,l6xrjt,fpcb2n, sz6ohy, yc1ezl,hts1gp,mk2vyr,
     &thfyl1,la5dcf)
      implicit logical (a-z)
      integer nfiumb4, nk, yc1ezl, hts1gp(3), mk2vyr, thfyl1, la5dcf
      double precision egoxa3, atqh0o, xs(nfiumb4), ys(nfiumb4), ws(
     &nfiumb4), knot(nk+4), coef(nk), sz(nfiumb4), rjcq9o(nfiumb4), 
     &n9peut, l6xrjt, fpcb2n(3), sz6ohy(1)
      call gzmfi3(egoxa3,atqh0o,xs,ys,ws, nfiumb4,nk, knot,coef,sz,
     &rjcq9o, n9peut,hts1gp(1),l6xrjt,hts1gp(2), hts1gp(3), fpcb2n(1),
     &fpcb2n(2),fpcb2n(3), yc1ezl, sz6ohy(1), sz6ohy(nk+1),sz6ohy(2*nk+
     &1),sz6ohy(3*nk+1),sz6ohy(4*nk+1), sz6ohy(5*nk+1),sz6ohy(6*nk+1),
     &sz6ohy(7*nk+1),sz6ohy(8*nk+1), sz6ohy(9*nk+1),sz6ohy(9*nk+mk2vyr*
     &nk+1),sz6ohy(9*nk+2*mk2vyr*nk+1), mk2vyr,thfyl1,la5dcf)
      return
      end
      subroutine gzmfi3(egoxa3,atqh0o,xs,ys,ws, nfiumb4,nk, knot,coef,
     &sz,rjcq9o, n9peut,rlhz2a,dwgkz6,ispar, pga6nu, lspar,uspar,fjo2dy,
     & yc1ezl, mheq6i, n7cuql,dvpc8x,hdv8br,cbg5ys, vf1jtn,eh6nly,
     &mvx9at,vbxpg4, rlep7v,lunah2,p2ip, mk2vyr,thfyl1,la5dcf)
      implicit logical (a-z)
      integer nfiumb4,nk, rlhz2a,ispar, yc1ezl, mk2vyr,thfyl1,la5dcf
      integer pga6nu
      double precision egoxa3,atqh0o,xs(nfiumb4),ys(nfiumb4),ws(nfiumb4)
     &, knot(nk+4), coef(nk),sz(nfiumb4),rjcq9o(nfiumb4), n9peut,dwgkz6,
     &lspar,uspar,fjo2dy, mheq6i(nk), n7cuql(nk),dvpc8x(nk),hdv8br(nk),
     &cbg5ys(nk), vf1jtn(nk),eh6nly(nk),mvx9at(nk),vbxpg4(nk), rlep7v(
     &mk2vyr,nk),lunah2(mk2vyr,nk),p2ip(thfyl1,nk)
      double precision t1,t2,dyb3po, a,b,c,d,e,kqoy6w,xm,p,q,r,fjo2dy1,
     &fjo2dy2,u,v,w, fu,fv,fw,fx,x, ax,bx
      integer w3gohz, vucgi1r
      double precision hz0fmy, epx9jf
      hz0fmy = 8.0d88
      epx9jf = 0.0d0
      d = 0.5d0
      u = 0.5d0
      dyb3po = 0.5d0
      w3gohz = 1
23000 if(.not.(w3gohz.le.nfiumb4))goto 23002
      if(.not.(ws(w3gohz).gt.0.0d0))goto 23003
      ws(w3gohz) = dsqrt(ws(w3gohz))
23003 continue
       w3gohz = w3gohz+1
      goto 23000
23002 continue
      if(.not.(yc1ezl .eq. 0))goto 23005
      call poqy8c(vf1jtn,eh6nly,mvx9at,vbxpg4,knot,nk)
      call ak9vxi(xs,ys,ws,knot, nfiumb4,nk, mheq6i,n7cuql,dvpc8x,
     &hdv8br,cbg5ys)
      t1 = 0.0d0 
      t2 = 0.0d0
      do 23007 w3gohz = 3,nk-3 
      t1 = t1 + n7cuql(w3gohz) 
23007 continue
      do 23009 w3gohz = 3,nk-3 
      t2 = t2 + vf1jtn(w3gohz) 
23009 continue
      dyb3po = t1/t2
      yc1ezl = 1
23005 continue
      if(.not.(ispar .eq. 1))goto 23011
      call oipu6h(egoxa3,atqh0o,xs,ys,ws, nfiumb4,nk,rlhz2a, knot,coef,
     &sz,rjcq9o,n9peut, dwgkz6, mheq6i, n7cuql,dvpc8x,hdv8br,cbg5ys, 
     &vf1jtn,eh6nly,mvx9at,vbxpg4, rlep7v,lunah2,p2ip,mk2vyr,thfyl1,
     &la5dcf)
      return
23011 continue
      ax = lspar 
      bx = uspar
      c = 0.381966011250105097d0
      kqoy6w = 2.0d-7
      vucgi1r = 0
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0d0
      dwgkz6 = dyb3po * dexp((-2.0d0 + x*6.0d0) * dlog(16.0d0))
      call oipu6h(egoxa3,atqh0o,xs,ys,ws, nfiumb4,nk,rlhz2a, knot,coef,
     &sz,rjcq9o,n9peut, dwgkz6, mheq6i, n7cuql,dvpc8x,hdv8br,cbg5ys, 
     &vf1jtn,eh6nly,mvx9at,vbxpg4, rlep7v,lunah2,p2ip,mk2vyr,thfyl1,
     &la5dcf)
      fx = n9peut
      fv = fx
      fw = fx
23013 if(.not.(la5dcf .eq. 0))goto 23014
      vucgi1r = vucgi1r + 1
      xm = 0.5d0*(a + b)
      fjo2dy1 = kqoy6w*dabs(x) + fjo2dy/3.0d0
      fjo2dy2 = 2.0d0*fjo2dy1
      if(.not.((dabs(x - xm) .le. (fjo2dy2 - 0.5d0*(b - a))) .or.(
     &vucgi1r .gt. pga6nu)))goto 23015
      go to 90
23015 continue
      if(.not.((dabs(e) .le. fjo2dy1) .or.(fx .ge. hz0fmy) .or.(fv .ge. 
     &hz0fmy) .or.(fw .ge. hz0fmy)))goto 23017
      go to 40
23017 continue
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.0d0 * (q - r)
      if(.not.(q .gt. 0.0d0))goto 23019
      p = -p
23019 continue
      q = dabs(q)
      r = e
      e = d
30    if(.not.((dabs(p) .ge. dabs(0.5d0*q*r)) .or.(q .eq. 0.0d0)))goto 2
     &3021
      go to 40
23021 continue
      if(.not.((p .le. q*(a - x)) .or. (p .ge. q*(b - x))))goto 23023
      go to 40
23023 continue
      d = p/q
      u = x + d
      if(.not.((u - a) .lt. fjo2dy2))goto 23025
      d = dsign(fjo2dy1, xm - x)
23025 continue
      if(.not.((b - u) .lt. fjo2dy2))goto 23027
      d = dsign(fjo2dy1, xm - x)
23027 continue
      go to 50
40    if(.not.(x .ge. xm))goto 23029
      e = a - x
      goto 23030
23029 continue
      e = b - x
23030 continue
      d = c*e
50    if(.not.(dabs(d) .ge. fjo2dy1))goto 23031
      u = x + d
      goto 23032
23031 continue
      u = x + dsign(fjo2dy1, d)
23032 continue
      dwgkz6 = dyb3po * dexp((-2.0d0 + u*6.0) * dlog(16.0d0))
      call oipu6h(egoxa3,atqh0o,xs,ys,ws, nfiumb4,nk,rlhz2a, knot,coef,
     &sz,rjcq9o,n9peut, dwgkz6, mheq6i, n7cuql,dvpc8x,hdv8br,cbg5ys, 
     &vf1jtn,eh6nly,mvx9at,vbxpg4, rlep7v,lunah2,p2ip,mk2vyr,thfyl1,
     &la5dcf)
      fu = n9peut
      if(.not.(fu .gt. hz0fmy))goto 23033
      fu = 2.0d0 * hz0fmy
23033 continue
      if(.not.(fu .le. fx))goto 23035
      if(.not.(u .ge. x))goto 23037
      a = x
      goto 23038
23037 continue
      b = x
23038 continue
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      goto 23036
23035 continue
      if(.not.(u .lt. x))goto 23039
      a = u
      goto 23040
23039 continue
      b = u
23040 continue
      if(.not.((fu .le. fw) .or. (w .eq. x)))goto 23041
      v = w
      fv = fw
      w = u
      fw = fu
      goto 23042
23041 continue
      if(.not.((fu .le. fv) .or. (v .eq. x) .or. (v .eq. w)))goto 23043
      v = u
      fv = fu
23043 continue
23042 continue
23036 continue
      goto 23013
23014 continue
90    epx9jf = 0.0d0
      dwgkz6 = dyb3po * dexp((-2.0d0 + x*6.0d0) * dlog(16.0d0))
      n9peut = fx
      return
23012 continue
      return
      end
      subroutine poqy8c(vf1jtn,eh6nly,mvx9at,vbxpg4,tb,nb)
      implicit logical (a-z)
      integer nb
      double precision vf1jtn(nb),eh6nly(nb),mvx9at(nb),vbxpg4(nb),tb(
     &nb+4)
      integer m5xudf,ilo,i6ndbu, ynmzp6, def4wn
      integer w3gohz,p1rifj,d9rjek
      double precision uq9jtc(4,3),bgu6fw(16),avoe4y(4),yw2(4), wpt
      double precision uoqx2m
      uoqx2m = 1.0d0 / 3.0d0
      ynmzp6 = 3
      def4wn = 4
      do 23045 w3gohz = 1,nb 
      vf1jtn(w3gohz) = 0.0d0
      eh6nly(w3gohz) = 0.0d0
      mvx9at(w3gohz) = 0.0d0
      vbxpg4(w3gohz) = 0.0d0 
23045 continue
      ilo = 1
      do 23047 w3gohz = 1,nb 
      call vinterv(tb(1),(nb+1),tb(w3gohz),m5xudf,i6ndbu)
      call vbsplvd(tb,def4wn,tb(w3gohz),m5xudf,bgu6fw,uq9jtc,ynmzp6)
      do 23049 p1rifj = 1,4 
      avoe4y(p1rifj) = uq9jtc(p1rifj,3) 
23049 continue
      call vbsplvd(tb,def4wn,tb(w3gohz+1),m5xudf,bgu6fw,uq9jtc,ynmzp6)
      do 23051 p1rifj = 1,4 
      yw2(p1rifj) = uq9jtc(p1rifj,3) - avoe4y(p1rifj) 
23051 continue
      wpt = tb(w3gohz+1) - tb(w3gohz)
      if(.not.(m5xudf .ge. 4))goto 23053
      do 23055 p1rifj = 1,4 
      d9rjek = p1rifj
      vf1jtn(m5xudf-4+p1rifj) = vf1jtn(m5xudf-4+p1rifj) + wpt * (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
      d9rjek = p1rifj+1
      if(.not.(d9rjek .le. 4))goto 23057
      eh6nly(m5xudf+p1rifj-4) = eh6nly(m5xudf+p1rifj-4) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23057 continue
      d9rjek = p1rifj+2
      if(.not.(d9rjek .le. 4))goto 23059
      mvx9at(m5xudf+p1rifj-4) = mvx9at(m5xudf+p1rifj-4) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23059 continue
      d9rjek = p1rifj+3
      if(.not.(d9rjek .le. 4))goto 23061
      vbxpg4(m5xudf+p1rifj-4) = vbxpg4(m5xudf+p1rifj-4) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23061 continue
23055 continue
      goto 23054
23053 continue
      if(.not.(m5xudf .eq. 3))goto 23063
      do 23065 p1rifj = 1,3 
      d9rjek = p1rifj
      vf1jtn(m5xudf-3+p1rifj) = vf1jtn(m5xudf-3+p1rifj) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
      d9rjek = p1rifj+1
      if(.not.(d9rjek .le. 3))goto 23067
      eh6nly(m5xudf+p1rifj-3) = eh6nly(m5xudf+p1rifj-3) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23067 continue
      d9rjek = p1rifj+2
      if(.not.(d9rjek .le. 3))goto 23069
      mvx9at(m5xudf+p1rifj-3) = mvx9at(m5xudf+p1rifj-3) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23069 continue
23065 continue
      goto 23064
23063 continue
      if(.not.(m5xudf .eq. 2))goto 23071
      do 23073 p1rifj = 1,2 
      d9rjek = p1rifj
      vf1jtn(m5xudf-2+p1rifj) = vf1jtn(m5xudf-2+p1rifj) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
      d9rjek = p1rifj+1
      if(.not.(d9rjek .le. 2))goto 23075
      eh6nly(m5xudf+p1rifj-2) = eh6nly(m5xudf+p1rifj-2) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23075 continue
23073 continue
      goto 23072
23071 continue
      if(.not.(m5xudf .eq. 1))goto 23077
      do 23079 p1rifj = 1,1 
      d9rjek = p1rifj
      vf1jtn(m5xudf-1+p1rifj) = vf1jtn(m5xudf-1+p1rifj) + wpt* (avoe4y(
     &p1rifj)*avoe4y(d9rjek) + (yw2(p1rifj)*avoe4y(d9rjek) + yw2(d9rjek)
     &*avoe4y(p1rifj))*0.50 + yw2(p1rifj)*yw2(d9rjek)*uoqx2m)
23079 continue
23077 continue
23072 continue
23064 continue
23054 continue
23047 continue
      return
      end
      subroutine gayot2(rlep7v,lunah2,p2ip, mk2vyr,nk,thfyl1,sbkvx6)
      implicit logical (a-z)
      integer mk2vyr,nk,thfyl1,sbkvx6
      double precision rlep7v(mk2vyr,nk), lunah2(mk2vyr,nk), p2ip(
     &thfyl1,nk)
      integer w3gohz, d9rjek, nd6mep
      double precision yrbij3(3),vef2gk(2),cfko0l(1),c0,c1,c2,c3
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      yrbij3(1) = 0.0d0
      yrbij3(2) = 0.0d0
      yrbij3(1) = 0.0d0
      vef2gk(1) = 0.0d0
      vef2gk(2) = 0.0d0
      cfko0l(1) = 0.0d0
      do 23081 w3gohz = 1,nk 
      d9rjek = nk-w3gohz+1
      c0 = 1.0d0 / rlep7v(4,d9rjek)
      if(.not.(d9rjek .le. nk-3))goto 23083
      c1 = rlep7v(1,d9rjek+3)*c0
      c2 = rlep7v(2,d9rjek+2)*c0
      c3 = rlep7v(3,d9rjek+1)*c0
      goto 23084
23083 continue
      if(.not.(d9rjek .eq. nk-2))goto 23085
      c1 = 0.0d0
      c2 = rlep7v(2,d9rjek+2)*c0
      c3 = rlep7v(3,d9rjek+1)*c0
      goto 23086
23085 continue
      if(.not.(d9rjek .eq. nk-1))goto 23087
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = rlep7v(3,d9rjek+1)*c0
      goto 23088
23087 continue
      if(.not.(d9rjek .eq. nk))goto 23089
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
23089 continue
23088 continue
23086 continue
23084 continue
      lunah2(1,d9rjek) = 0.0d0 - (c1*yrbij3(1)+c2*yrbij3(2)+c3*yrbij3(3)
     &)
      lunah2(2,d9rjek) = 0.0d0 - (c1*yrbij3(2)+c2*vef2gk(1)+c3*vef2gk(2)
     &)
      lunah2(3,d9rjek) = 0.0d0 - (c1*yrbij3(3)+c2*vef2gk(2)+c3*cfko0l(1)
     &)
      lunah2(4,d9rjek) = c0**2 + c1**2 * yrbij3(1) + 2.0d0*c1*c2*yrbij3(
     &2)+2.0d0*c1*c3*yrbij3(3) + c2**2 * vef2gk(1) + 2.0d0*c2*c3*vef2gk(
     &2) + c3**2 * cfko0l(1)
      yrbij3(1) = vef2gk(1)
      yrbij3(2) = vef2gk(2)
      yrbij3(3) = lunah2(2,d9rjek)
      vef2gk(1) = cfko0l(1)
      vef2gk(2) = lunah2(3,d9rjek)
      cfko0l(1) = lunah2(4,d9rjek)
23081 continue
      if(.not.(sbkvx6 .eq. 0))goto 23091
      return
23091 continue
      do 23093 w3gohz = 1,nk 
      d9rjek = nk-w3gohz+1
      nd6mep = 1
23095 if(.not.(nd6mep.le.4.and.d9rjek+nd6mep-1.le.nk))goto 23097
      p2ip(d9rjek,d9rjek+nd6mep-1) = lunah2(5-nd6mep,d9rjek)
       nd6mep = nd6mep+1
      goto 23095
23097 continue
23093 continue
      do 23098 w3gohz = 1,nk 
      d9rjek = nk-w3gohz+1
      nd6mep = d9rjek-4
23100 if(.not.(nd6mep.ge.1))goto 23102
      c0 = 1.0 / rlep7v(4,nd6mep) 
      c1 = rlep7v(1,nd6mep+3)*c0
      c2 = rlep7v(2,nd6mep+2)*c0 
      c3 = rlep7v(3,nd6mep+1)*c0
      p2ip(nd6mep,d9rjek) = 0.0d0- ( c1*p2ip(nd6mep+3,d9rjek) + c2*p2ip(
     &nd6mep+2,d9rjek) + c3*p2ip(nd6mep+1,d9rjek) )
       nd6mep = nd6mep-1
      goto 23100
23102 continue
23098 continue
      return
23092 continue
      end
      subroutine oipu6h(egoxa3,atqh0o,x,y,w, nfiumb4,nk,rlhz2a, knot,
     &coef,sz,rjcq9o, n9peut, dwgkz6, mheq6i, n7cuql,dvpc8x,hdv8br,
     &cbg5ys, vf1jtn,eh6nly,mvx9at,vbxpg4, rlep7v,lunah2,p2ip, mk2vyr,
     &thfyl1,fjg0qv)
      implicit logical (a-z)
      integer nfiumb4,nk,rlhz2a, mk2vyr,thfyl1,fjg0qv
      double precision egoxa3,atqh0o,x(nfiumb4),y(nfiumb4),w(nfiumb4), 
     &knot(nk+4), coef(nk),sz(nfiumb4),rjcq9o(nfiumb4), n9peut, dwgkz6, 
     &mheq6i(nk), n7cuql(nk),dvpc8x(nk),hdv8br(nk),cbg5ys(nk), vf1jtn(
     &nk),eh6nly(nk),mvx9at(nk),vbxpg4(nk), rlep7v(mk2vyr,nk),lunah2(
     &mk2vyr,nk),p2ip(thfyl1,nk)
      double precision das4bx, bgu6fw(16), b0,b1,b2,b3,kqoy6w, uq9jtc(4,
     &1), xv,bvalue,df
      double precision risyv0
      integer oht3ga, ynmzp6, ilo, i6ndbu, d9rjek, w3gohz, px1yhr, 
     &m5xudf, def4wn
      ilo = 1
      kqoy6w = 0.1d-10
      oht3ga = 0
      ynmzp6 = 3
      def4wn = 4
      do 23103 w3gohz = 1,nk 
      coef(w3gohz) = mheq6i(w3gohz) 
23103 continue
      do 23105 w3gohz = 1,nk 
      rlep7v(4,w3gohz) = n7cuql(w3gohz)+dwgkz6*vf1jtn(w3gohz) 
23105 continue
      do 23107 w3gohz = 1,(nk-1) 
      rlep7v(3,w3gohz+1) = dvpc8x(w3gohz)+dwgkz6*eh6nly(w3gohz) 
23107 continue
      do 23109 w3gohz = 1,(nk-2) 
      rlep7v(2,w3gohz+2) = hdv8br(w3gohz)+dwgkz6*mvx9at(w3gohz) 
23109 continue
      do 23111 w3gohz = 1,(nk-3) 
      rlep7v(1,w3gohz+3) = cbg5ys(w3gohz)+dwgkz6*vbxpg4(w3gohz) 
23111 continue
      call dpbfa8(rlep7v,mk2vyr,nk,ynmzp6,fjg0qv)
      if(.not.(fjg0qv .ne. 0))goto 23113
      return
23113 continue
      call dpbsl8(rlep7v,mk2vyr,nk,ynmzp6,coef)
      px1yhr = 1
      do 23115 w3gohz = 1,nfiumb4 
      xv = x(w3gohz)
      sz(w3gohz) = bvalue(knot,coef, nk,def4wn,xv,oht3ga)
23115 continue
      if(.not.(rlhz2a .eq. 0))goto 23117
      return
23117 continue
      call gayot2(rlep7v,lunah2,p2ip, mk2vyr,nk,thfyl1,oht3ga)
      do 23119 w3gohz = 1,nfiumb4 
      xv = x(w3gohz)
      call vinterv(knot(1),(nk+1),xv,m5xudf,i6ndbu)
      if(.not.(i6ndbu .eq. -1))goto 23121
      m5xudf = 4 
      xv = knot(4) + kqoy6w 
23121 continue
      if(.not.(i6ndbu .eq. 1))goto 23123
      m5xudf = nk 
      xv = knot(nk+1) - kqoy6w 
23123 continue
      d9rjek = m5xudf-3
      call vbsplvd(knot,4,xv,m5xudf,bgu6fw,uq9jtc,1)
      b0 = uq9jtc(1,1)
      b1 = uq9jtc(2,1)
      b2 = uq9jtc(3,1)
      b3 = uq9jtc(4,1)
      rjcq9o(w3gohz) = (lunah2(4,d9rjek)*b0**2 + 2.0d0*lunah2(3,d9rjek)*
     &b0*b1 + 2.0d0*lunah2(2,d9rjek)*b0*b2 + 2.0d0*lunah2(1,d9rjek)*b0*
     &b3 + lunah2(4,d9rjek+1)*b1**2 + 2.0d0*lunah2(3,d9rjek+1)*b1*b2 + 
     &2.0d0*lunah2(2,d9rjek+1)*b1*b3 + lunah2(4,d9rjek+2)*b2**2 + 2.0d0*
     &lunah2(3,d9rjek+2)*b2*b3 + lunah2(4,d9rjek+3)*b3**2 ) * w(w3gohz)*
     &*2
23119 continue
      if(.not.(rlhz2a .eq. 1))goto 23125
      das4bx = 0.0d0 
      df = 0.0d0 
      risyv0 = 0.0d0
      do 23127 w3gohz = 1,nfiumb4 
      das4bx = das4bx + ((y(w3gohz)-sz(w3gohz))*w(w3gohz))**2
      df = df + rjcq9o(w3gohz)
      risyv0 = risyv0 + w(w3gohz)*w(w3gohz)
23127 continue
      n9peut = (das4bx/risyv0)/((1.0d0-(atqh0o+egoxa3*df)/risyv0)**2)
      goto 23126
23125 continue
      if(.not.(rlhz2a .eq. 2))goto 23129
      n9peut = 0.0d0
      risyv0 = 0.0d0
      do 23131 w3gohz = 1,nfiumb4 
      n9peut = n9peut + (((y(w3gohz)-sz(w3gohz))*w(w3gohz))/(1.0d0-
     &rjcq9o(w3gohz)))**2
      risyv0 = risyv0 + w(w3gohz)*w(w3gohz)
23131 continue
      n9peut = n9peut / risyv0
      goto 23130
23129 continue
      n9peut = 0.0d0
      do 23133 w3gohz = 1,nfiumb4 
      n9peut = n9peut+rjcq9o(w3gohz)
23133 continue
      n9peut = 3.0d0 + (atqh0o-n9peut)**2
23130 continue
23126 continue
      return
23118 continue
      end
      subroutine ak9vxi(p3vlea,hr83e,w,onyz6j, xl6qgm,nfiumb4, wevr5o,
     &n7cuql,dvpc8x,hdv8br,cbg5ys)
      implicit logical (a-z)
      integer xl6qgm,nfiumb4
      double precision p3vlea(xl6qgm),hr83e(xl6qgm),w(xl6qgm),onyz6j(
     &nfiumb4+4), wevr5o(nfiumb4), n7cuql(nfiumb4),dvpc8x(nfiumb4),
     &hdv8br(nfiumb4),cbg5ys(nfiumb4)
      double precision kqoy6w,uq9jtc(4,1),bgu6fw(16)
      integer d9rjek,w3gohz,ilo,m5xudf,i6ndbu
      do 23135 w3gohz = 1,nfiumb4 
      wevr5o(w3gohz) = 0.0d0 
      n7cuql(w3gohz) = 0.0d0 
      dvpc8x(w3gohz) = 0.0d0
      hdv8br(w3gohz) = 0.0d0 
      cbg5ys(w3gohz) = 0.0d0
23135 continue
      ilo = 1
      kqoy6w = 0.1d-9
      do 23137 w3gohz = 1,xl6qgm 
      call vinterv(onyz6j(1),(nfiumb4+1),p3vlea(w3gohz),m5xudf,i6ndbu)
      if(.not.(i6ndbu .eq. 1))goto 23139
      if(.not.(p3vlea(w3gohz) .le. (onyz6j(m5xudf)+kqoy6w)))goto 23141
      m5xudf = m5xudf-1
      goto 23142
23141 continue
      return
23142 continue
23139 continue
      call vbsplvd(onyz6j,4,p3vlea(w3gohz),m5xudf,bgu6fw,uq9jtc,1)
      d9rjek = m5xudf-4+1
      wevr5o(d9rjek) = wevr5o(d9rjek)+w(w3gohz)**2*hr83e(w3gohz)*uq9jtc(
     &1,1)
      n7cuql(d9rjek) = n7cuql(d9rjek)+w(w3gohz)**2*uq9jtc(1,1)**2
      dvpc8x(d9rjek) = dvpc8x(d9rjek)+w(w3gohz)**2*uq9jtc(1,1)*uq9jtc(2,
     &1)
      hdv8br(d9rjek) = hdv8br(d9rjek)+w(w3gohz)**2*uq9jtc(1,1)*uq9jtc(3,
     &1)
      cbg5ys(d9rjek) = cbg5ys(d9rjek)+w(w3gohz)**2*uq9jtc(1,1)*uq9jtc(4,
     &1)
      d9rjek = m5xudf-4+2
      wevr5o(d9rjek) = wevr5o(d9rjek)+w(w3gohz)**2*hr83e(w3gohz)*uq9jtc(
     &2,1)
      n7cuql(d9rjek) = n7cuql(d9rjek)+w(w3gohz)**2*uq9jtc(2,1)**2
      dvpc8x(d9rjek) = dvpc8x(d9rjek)+w(w3gohz)**2*uq9jtc(2,1)*uq9jtc(3,
     &1)
      hdv8br(d9rjek) = hdv8br(d9rjek)+w(w3gohz)**2*uq9jtc(2,1)*uq9jtc(4,
     &1)
      d9rjek = m5xudf-4+3
      wevr5o(d9rjek) = wevr5o(d9rjek)+w(w3gohz)**2*hr83e(w3gohz)*uq9jtc(
     &3,1)
      n7cuql(d9rjek) = n7cuql(d9rjek)+w(w3gohz)**2*uq9jtc(3,1)**2
      dvpc8x(d9rjek) = dvpc8x(d9rjek)+w(w3gohz)**2*uq9jtc(3,1)*uq9jtc(4,
     &1)
      d9rjek = m5xudf-4+4
      wevr5o(d9rjek) = wevr5o(d9rjek)+w(w3gohz)**2*hr83e(w3gohz)*uq9jtc(
     &4,1)
      n7cuql(d9rjek) = n7cuql(d9rjek)+w(w3gohz)**2*uq9jtc(4,1)**2
23137 continue
      return
      end
