      subroutine dnaoqj0l(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk, 
     &ankcghz2,coef,sz,ifys6woa, qcpiaj7f,wbkq9zyi,parms, scrtch, 
     &gp0xjetb,l3zpbstu,e5knafcg,wep0oibc,fbd5yktj)
      implicit logical (a-z)
      integer kuzxj1lo, nk, gp0xjetb, l3zpbstu(3), e5knafcg, wep0oibc, 
     &fbd5yktj
      double precision penalt, pjb6wfoq, xs(kuzxj1lo), ys(kuzxj1lo), ws(
     &kuzxj1lo), ankcghz2(nk+4), coef(nk), sz(kuzxj1lo), ifys6woa(
     &kuzxj1lo), qcpiaj7f, wbkq9zyi, parms(3), scrtch(*)
      call hbzuprs6(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk, ankcghz2,
     &coef,sz,ifys6woa, qcpiaj7f,l3zpbstu(1),wbkq9zyi,l3zpbstu(2), 
     &l3zpbstu(3), parms(1),parms(2),parms(3), gp0xjetb, scrtch(1), 
     &scrtch(nk+1),scrtch(2*nk+1),scrtch(3*nk+1),scrtch(4*nk+1), scrtch(
     &5*nk+1),scrtch(6*nk+1),scrtch(7*nk+1),scrtch(8*nk+1), scrtch(9*nk+
     &1),scrtch(9*nk+e5knafcg*nk+1),scrtch(9*nk+2*e5knafcg*nk+1), 
     &e5knafcg,wep0oibc,fbd5yktj)
      return
      end
      subroutine hbzuprs6(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk, 
     &ankcghz2,coef,sz,ifys6woa, qcpiaj7f,icrit,i9mwnvqt,ispar, 
     &c5aesxku, mynl7uaq,zustx4fw,tol, gp0xjetb, xwy, zvau2lct,f6lsuzax,
     &fvh2rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, buhyalv4,
     &fulcp8wa,plj0trqx, e5knafcg,wep0oibc,fbd5yktj)
      implicit logical (a-z)
      integer kuzxj1lo,nk, icrit,ispar, gp0xjetb, e5knafcg,wep0oibc,
     &fbd5yktj
      integer c5aesxku
      double precision penalt,pjb6wfoq,xs(kuzxj1lo),ys(kuzxj1lo),ws(
     &kuzxj1lo), ankcghz2(nk+4), coef(nk),sz(kuzxj1lo),ifys6woa(
     &kuzxj1lo), qcpiaj7f,i9mwnvqt,mynl7uaq,zustx4fw,tol, xwy(nk), 
     &zvau2lct(nk),f6lsuzax(nk),fvh2rwtc(nk),dcfir2no(nk), xecbg0pf(nk),
     &z4grbpiq(nk),d7glzhbj(nk),v2eydbxs(nk), buhyalv4(e5knafcg,nk),
     &fulcp8wa(e5knafcg,nk),plj0trqx(wep0oibc,nk)
      double precision t1,t2,ratio, a,b,c,d,e,qaltf0nz,xm,p,q,r,tol1,
     &tol2,u,v,w, fu,fv,fw,fx,x, ax,bx
      integer ayfnwr1v, viter
      double precision yjpnro8d, hmayv1xt
      yjpnro8d = 8.0d88
      hmayv1xt = 0.0d0
      d = 0.5d0
      u = 0.5d0
      ratio = 0.5d0
      ayfnwr1v = 1
23000 if(.not.(ayfnwr1v.le.kuzxj1lo))goto 23002
      if(.not.(ws(ayfnwr1v).gt.0.0d0))goto 23003
      ws(ayfnwr1v) = dsqrt(ws(ayfnwr1v))
23003 continue
       ayfnwr1v = ayfnwr1v+1
      goto 23000
23002 continue
      if(.not.(gp0xjetb .eq. 0))goto 23005
      call zosq7hub(xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs,ankcghz2,nk)
      call gt9iulbf(xs,ys,ws,ankcghz2, kuzxj1lo,nk, xwy,zvau2lct,
     &f6lsuzax,fvh2rwtc,dcfir2no)
      t1 = 0.0d0 
      t2 = 0.0d0
      do 23007 ayfnwr1v = 3,nk-3 
      t1 = t1 + zvau2lct(ayfnwr1v) 
23007 continue
      do 23009 ayfnwr1v = 3,nk-3 
      t2 = t2 + xecbg0pf(ayfnwr1v) 
23009 continue
      ratio = t1/t2
      gp0xjetb = 1
23005 continue
      if(.not.(ispar .eq. 1))goto 23011
      call wmhctl9x(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk,icrit, 
     &ankcghz2,coef,sz,ifys6woa,qcpiaj7f, i9mwnvqt, xwy, zvau2lct,
     &f6lsuzax,fvh2rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, 
     &buhyalv4,fulcp8wa,plj0trqx,e5knafcg,wep0oibc,fbd5yktj)
      return
23011 continue
      ax = mynl7uaq 
      bx = zustx4fw
      c = 0.381966011250105097d0
      qaltf0nz = 2.0d-5
      viter = 0
      a = ax
      b = bx
      v = a + c*(b - a)
      w = v
      x = v
      e = 0.0d0
      i9mwnvqt = ratio * dexp((-2.0d0 + x*6.0d0) * dlog(16.0d0))
      call wmhctl9x(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk,icrit, 
     &ankcghz2,coef,sz,ifys6woa,qcpiaj7f, i9mwnvqt, xwy, zvau2lct,
     &f6lsuzax,fvh2rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, 
     &buhyalv4,fulcp8wa,plj0trqx,e5knafcg,wep0oibc,fbd5yktj)
      fx = qcpiaj7f
      fv = fx
      fw = fx
23013 if(.not.(fbd5yktj .eq. 0))goto 23014
      viter = viter + 1
      xm = 0.5d0*(a + b)
      tol1 = qaltf0nz*dabs(x) + tol/3.0d0
      tol2 = 2.0d0*tol1
      if(.not.((dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) .or.(viter 
     &.gt. c5aesxku)))goto 23015
      go to 90
23015 continue
      if(.not.((dabs(e) .le. tol1) .or.(fx .ge. yjpnro8d) .or.(fv .ge. 
     &yjpnro8d) .or.(fw .ge. yjpnro8d)))goto 23017
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
      if(.not.((dabs(p) .ge. dabs(0.5d0*q*r)) .or.(q .eq. 0.0d0)))goto 2
     &3021
      go to 40
23021 continue
      if(.not.((p .le. q*(a - x)) .or. (p .ge. q*(b - x))))goto 23023
      go to 40
23023 continue
      d = p/q
      u = x + d
      if(.not.((u - a) .lt. tol2))goto 23025
      d = dsign(tol1, xm - x)
23025 continue
      if(.not.((b - u) .lt. tol2))goto 23027
      d = dsign(tol1, xm - x)
23027 continue
      go to 50
40    if(.not.(x .ge. xm))goto 23029
      e = a - x
      goto 23030
23029 continue
      e = b - x
23030 continue
      d = c*e
50    if(.not.(dabs(d) .ge. tol1))goto 23031
      u = x + d
      goto 23032
23031 continue
      u = x + dsign(tol1, d)
23032 continue
      i9mwnvqt = ratio * dexp((-2.0d0 + u*6.0) * dlog(16.0d0))
      call wmhctl9x(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk,icrit, 
     &ankcghz2,coef,sz,ifys6woa,qcpiaj7f, i9mwnvqt, xwy, zvau2lct,
     &f6lsuzax,fvh2rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, 
     &buhyalv4,fulcp8wa,plj0trqx,e5knafcg,wep0oibc,fbd5yktj)
      fu = qcpiaj7f
      if(.not.(fu .gt. yjpnro8d))goto 23033
      fu = 2.0d0 * yjpnro8d
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
90    hmayv1xt = 0.0d0
      i9mwnvqt = ratio * dexp((-2.0d0 + x*6.0d0) * dlog(16.0d0))
      qcpiaj7f = fx
      return
      return
      end
      subroutine zosq7hub(xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs,tb,nb)
      implicit logical (a-z)
      integer nb
      double precision xecbg0pf(nb),z4grbpiq(nb),d7glzhbj(nb),v2eydbxs(
     &nb),tb(nb+4)
      integer dqlr5bse,ilo,pqzfxw4i, three3, ifour4, nbp1
      integer ayfnwr1v,iii,yq6lorbx
      integer i2svdbx3tk
      double precision g9fvdrbw(4,3),work(16),yw1(4),yw2(4), wpt
      double precision othird
      othird = 1.0d0 / 3.0d0
      three3 = 3
      ifour4 = 4
      nbp1 = nb + 1
      do 23045 ayfnwr1v = 1,nb 
      xecbg0pf(ayfnwr1v) = 0.0d0
      z4grbpiq(ayfnwr1v) = 0.0d0
      d7glzhbj(ayfnwr1v) = 0.0d0
      v2eydbxs(ayfnwr1v) = 0.0d0 
23045 continue
      ilo = 1
      do 23047 ayfnwr1v = 1,nb 
      call vinterv(tb(1), nbp1 ,tb(ayfnwr1v),dqlr5bse,pqzfxw4i)
      call vbsplvd(tb,ifour4,tb(ayfnwr1v),dqlr5bse,work,g9fvdrbw,three3)
      do 23049 iii = 1,4 
      yw1(iii) = g9fvdrbw(iii,3) 
23049 continue
      call vbsplvd(tb,ifour4,tb(ayfnwr1v+1),dqlr5bse,work,g9fvdrbw,
     &three3)
      do 23051 iii = 1,4 
      yw2(iii) = g9fvdrbw(iii,3) - yw1(iii) 
23051 continue
      wpt = tb(ayfnwr1v+1) - tb(ayfnwr1v)
      if(.not.(dqlr5bse .ge. 4))goto 23053
      do 23055 iii = 1,4 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-4+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt * (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
      yq6lorbx = iii+1
      if(.not.(yq6lorbx .le. 4))goto 23057
      z4grbpiq(i2svdbx3tk) = z4grbpiq(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23057 continue
      yq6lorbx = iii+2
      if(.not.(yq6lorbx .le. 4))goto 23059
      d7glzhbj(i2svdbx3tk) = d7glzhbj(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23059 continue
      yq6lorbx = iii+3
      if(.not.(yq6lorbx .le. 4))goto 23061
      v2eydbxs(i2svdbx3tk) = v2eydbxs(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23061 continue
23055 continue
      goto 23054
23053 continue
      if(.not.(dqlr5bse .eq. 3))goto 23063
      do 23065 iii = 1,3 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-3+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
      yq6lorbx = iii+1
      if(.not.(yq6lorbx .le. 3))goto 23067
      z4grbpiq(i2svdbx3tk) = z4grbpiq(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23067 continue
      yq6lorbx = iii+2
      if(.not.(yq6lorbx .le. 3))goto 23069
      d7glzhbj(i2svdbx3tk) = d7glzhbj(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23069 continue
23065 continue
      goto 23064
23063 continue
      if(.not.(dqlr5bse .eq. 2))goto 23071
      do 23073 iii = 1,2 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-2+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
      yq6lorbx = iii+1
      if(.not.(yq6lorbx .le. 2))goto 23075
      z4grbpiq(i2svdbx3tk) = z4grbpiq(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23075 continue
23073 continue
      goto 23072
23071 continue
      if(.not.(dqlr5bse .eq. 1))goto 23077
      do 23079 iii = 1,1 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-1+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt* (yw1(iii)*yw1(
     &yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.
     &50 + yw2(iii)*yw2(yq6lorbx)*othird)
23079 continue
23077 continue
23072 continue
23064 continue
23054 continue
23047 continue
      return
      end
      subroutine vmnweiy2(buhyalv4,fulcp8wa,plj0trqx, e5knafcg,nk,
     &wep0oibc,iflag)
      implicit logical (a-z)
      integer e5knafcg,nk,wep0oibc,iflag
      double precision buhyalv4(e5knafcg,nk), fulcp8wa(e5knafcg,nk), 
     &plj0trqx(wep0oibc,nk)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh
      double precision wjm3(3),wjm2(2),wjm1(1),c0,c1,c2,c3
      double precision pcsuow9k, qdbgu6oi, upwkh5xz, rul5fnyd, ueydbrg6,
     & plce2srm, k3yvomnh, bfdjhu7l, ctfvwdu0
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      wjm3(1) = 0.0d0
      wjm3(2) = 0.0d0
      wjm3(3) = 0.0d0
      wjm2(1) = 0.0d0
      wjm2(2) = 0.0d0
      wjm1(1) = 0.0d0
      do 23081 ayfnwr1v = 1,nk 
      yq6lorbx = nk-ayfnwr1v+1
      c0 = 1.0d0 / buhyalv4(4,yq6lorbx)
      if(.not.(yq6lorbx .le. (nk-3)))goto 23083
      c1 = buhyalv4(1,yq6lorbx+3)*c0
      c2 = buhyalv4(2,yq6lorbx+2)*c0
      c3 = buhyalv4(3,yq6lorbx+1)*c0
      goto 23084
23083 continue
      if(.not.(yq6lorbx .eq. (nk-2)))goto 23085
      c1 = 0.0d0
      c2 = buhyalv4(2,yq6lorbx+2)*c0
      c3 = buhyalv4(3,yq6lorbx+1)*c0
      goto 23086
23085 continue
      if(.not.(yq6lorbx .eq. (nk-1)))goto 23087
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = buhyalv4(3,yq6lorbx+1)*c0
      goto 23088
23087 continue
      if(.not.(yq6lorbx .eq. nk))goto 23089
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
23089 continue
23088 continue
23086 continue
23084 continue
      pcsuow9k = c1*wjm3(1)
      qdbgu6oi = c2*wjm3(2)
      upwkh5xz = c3*wjm3(3)
      rul5fnyd = c1*wjm3(2)
      ueydbrg6 = c2*wjm2(1)
      plce2srm = c3*wjm2(2)
      k3yvomnh = c1*wjm3(3)
      bfdjhu7l = c2*wjm2(2)
      ctfvwdu0 = c3*wjm1(1)
      fulcp8wa(1,yq6lorbx) = 0.0d0 - (pcsuow9k+qdbgu6oi+upwkh5xz)
      fulcp8wa(2,yq6lorbx) = 0.0d0 - (rul5fnyd+ueydbrg6+plce2srm)
      fulcp8wa(3,yq6lorbx) = 0.0d0 - (k3yvomnh+bfdjhu7l+ctfvwdu0)
      fulcp8wa(4,yq6lorbx) = c0**2 + c1*(pcsuow9k + 2.0d0*(qdbgu6oi + 
     &upwkh5xz)) + c2*(ueydbrg6 + 2.0d0* plce2srm) + c3*ctfvwdu0
      wjm3(1) = wjm2(1)
      wjm3(2) = wjm2(2)
      wjm3(3) = fulcp8wa(2,yq6lorbx)
      wjm2(1) = wjm1(1)
      wjm2(2) = fulcp8wa(3,yq6lorbx)
      wjm1(1) = fulcp8wa(4,yq6lorbx)
23081 continue
      if(.not.(iflag .eq. 0))goto 23091
      return
23091 continue
      do 23093 ayfnwr1v = 1,nk 
      yq6lorbx = nk-ayfnwr1v+1
      gp1jxzuh = 1
23095 if(.not.(gp1jxzuh.le.4.and.yq6lorbx+gp1jxzuh-1.le.nk))goto 23097
      plj0trqx(yq6lorbx,yq6lorbx+gp1jxzuh-1) = fulcp8wa(5-gp1jxzuh,
     &yq6lorbx)
       gp1jxzuh = gp1jxzuh+1
      goto 23095
23097 continue
23093 continue
      do 23098 ayfnwr1v = 1,nk 
      yq6lorbx = nk-ayfnwr1v+1
      gp1jxzuh = yq6lorbx-4
23100 if(.not.(gp1jxzuh.ge.1))goto 23102
      c0 = 1.0 / buhyalv4(4,gp1jxzuh) 
      c1 = buhyalv4(1,gp1jxzuh+3)*c0
      c2 = buhyalv4(2,gp1jxzuh+2)*c0 
      c3 = buhyalv4(3,gp1jxzuh+1)*c0
      plj0trqx(gp1jxzuh,yq6lorbx) = 0.0d0- ( c1*plj0trqx(gp1jxzuh+3,
     &yq6lorbx) + c2*plj0trqx(gp1jxzuh+2,yq6lorbx) + c3*plj0trqx(
     &gp1jxzuh+1,yq6lorbx) )
       gp1jxzuh = gp1jxzuh-1
      goto 23100
23102 continue
23098 continue
      return
      end
      subroutine wmhctl9x(penalt,pjb6wfoq,x,y,w, kuzxj1lo,nk,icrit, 
     &ankcghz2,coef,sz,ifys6woa, qcpiaj7f, i9mwnvqt, xwy, zvau2lct,
     &f6lsuzax,fvh2rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, 
     &buhyalv4,fulcp8wa,plj0trqx, e5knafcg,wep0oibc,info)
      implicit logical (a-z)
      integer kuzxj1lo,nk,icrit, e5knafcg,wep0oibc,info
      double precision penalt,pjb6wfoq,x(kuzxj1lo),y(kuzxj1lo),w(
     &kuzxj1lo)
      double precision ankcghz2(nk+4), coef(nk),sz(kuzxj1lo),ifys6woa(
     &kuzxj1lo), qcpiaj7f, i9mwnvqt, xwy(nk)
      double precision zvau2lct(nk),f6lsuzax(nk),fvh2rwtc(nk),dcfir2no(
     &nk)
      double precision xecbg0pf(nk),z4grbpiq(nk),d7glzhbj(nk),v2eydbxs(
     &nk), buhyalv4(e5knafcg,nk),fulcp8wa(e5knafcg,nk),plj0trqx(
     &wep0oibc,nk)
      double precision resss, work(16), b0,b1,b2,b3,qaltf0nz, g9fvdrbw(
     &4,1), xv,eqdf
      double precision qtce8hzo
      double precision rxeqjn0y
      integer izero0, three3, ilo, pqzfxw4i, yq6lorbx, ayfnwr1v
      integer icoef, dqlr5bse, ifour4, hbsl0gto, nkp1
      ilo = 1
      qaltf0nz = 0.1d-10
      izero0 = 0
      three3 = 3
      ifour4 = 4
      hbsl0gto = 1
      nkp1 = nk + 1
      do 23103 ayfnwr1v = 1,nk 
      coef(ayfnwr1v) = xwy(ayfnwr1v) 
23103 continue
      do 23105 ayfnwr1v = 1,nk 
      buhyalv4(4,ayfnwr1v) = zvau2lct(ayfnwr1v)+i9mwnvqt*xecbg0pf(
     &ayfnwr1v) 
23105 continue
      do 23107 ayfnwr1v = 1,(nk-1) 
      buhyalv4(3,ayfnwr1v+1) = f6lsuzax(ayfnwr1v)+i9mwnvqt*z4grbpiq(
     &ayfnwr1v) 
23107 continue
      do 23109 ayfnwr1v = 1,(nk-2) 
      buhyalv4(2,ayfnwr1v+2) = fvh2rwtc(ayfnwr1v)+i9mwnvqt*d7glzhbj(
     &ayfnwr1v) 
23109 continue
      do 23111 ayfnwr1v = 1,(nk-3) 
      buhyalv4(1,ayfnwr1v+3) = dcfir2no(ayfnwr1v)+i9mwnvqt*v2eydbxs(
     &ayfnwr1v) 
23111 continue
      call dpbfa8(buhyalv4,e5knafcg,nk,three3,info)
      if(.not.(info .ne. 0))goto 23113
      return
23113 continue
      call dpbsl8(buhyalv4,e5knafcg,nk,three3,coef)
      icoef = 1
      do 23115 ayfnwr1v = 1,kuzxj1lo 
      xv = x(ayfnwr1v)
      call wbvalue(ankcghz2,coef, nk,ifour4,xv,izero0, sz(ayfnwr1v))
23115 continue
      if(.not.(icrit .eq. 0))goto 23117
      return
23117 continue
      call vmnweiy2(buhyalv4,fulcp8wa,plj0trqx, e5knafcg,nk,wep0oibc,
     &izero0)
      do 23119 ayfnwr1v = 1,kuzxj1lo 
      xv = x(ayfnwr1v)
      call vinterv(ankcghz2(1), nkp1 ,xv,dqlr5bse,pqzfxw4i)
      if(.not.(pqzfxw4i .eq. -1))goto 23121
      dqlr5bse = 4 
      xv = ankcghz2(4) + qaltf0nz 
23121 continue
      if(.not.(pqzfxw4i .eq. 1))goto 23123
      dqlr5bse = nk 
      xv = ankcghz2(nk+1) - qaltf0nz 
23123 continue
      yq6lorbx = dqlr5bse-3
      call vbsplvd(ankcghz2,ifour4,xv,dqlr5bse,work,g9fvdrbw,hbsl0gto)
      b0 = g9fvdrbw(1,1)
      b1 = g9fvdrbw(2,1)
      b2 = g9fvdrbw(3,1)
      b3 = g9fvdrbw(4,1)
      qtce8hzo = (b0 *(fulcp8wa(4,yq6lorbx)*b0 + 2.0d0*(fulcp8wa(3,
     &yq6lorbx)*b1 + fulcp8wa(2,yq6lorbx)*b2 + fulcp8wa(1,yq6lorbx)*b3))
     & + b1 *(fulcp8wa(4,yq6lorbx+1)*b1 + 2.0d0*(fulcp8wa(3,yq6lorbx+1)*
     &b2 + fulcp8wa(2,yq6lorbx+1)*b3)) + b2 *(fulcp8wa(4,yq6lorbx+2)*b2 
     &+ 2.0d0* fulcp8wa(3,yq6lorbx+2)*b3 )+ b3**2* fulcp8wa(4,yq6lorbx+
     &3)) * w(ayfnwr1v)**2
      ifys6woa(ayfnwr1v) = qtce8hzo
23119 continue
      if(.not.(icrit .eq. 1))goto 23125
      resss = 0.0d0 
      eqdf = 0.0d0 
      rxeqjn0y = 0.0d0
      do 23127 ayfnwr1v = 1,kuzxj1lo 
      resss = resss + ((y(ayfnwr1v)-sz(ayfnwr1v))*w(ayfnwr1v))**2
      eqdf = eqdf + ifys6woa(ayfnwr1v)
      rxeqjn0y = rxeqjn0y + w(ayfnwr1v)*w(ayfnwr1v)
23127 continue
      qcpiaj7f = (resss/rxeqjn0y)/((1.0d0-(pjb6wfoq+penalt*eqdf)/
     &rxeqjn0y)**2)
      goto 23126
23125 continue
      if(.not.(icrit .eq. 2))goto 23129
      qcpiaj7f = 0.0d0
      rxeqjn0y = 0.0d0
      do 23131 ayfnwr1v = 1,kuzxj1lo 
      qcpiaj7f = qcpiaj7f + (((y(ayfnwr1v)-sz(ayfnwr1v))*w(ayfnwr1v))/(
     &1.0d0-ifys6woa(ayfnwr1v)))**2
      rxeqjn0y = rxeqjn0y + w(ayfnwr1v)*w(ayfnwr1v)
23131 continue
      qcpiaj7f = qcpiaj7f / rxeqjn0y
      goto 23130
23129 continue
      qcpiaj7f = 0.0d0
      do 23133 ayfnwr1v = 1,kuzxj1lo 
      qcpiaj7f = qcpiaj7f+ifys6woa(ayfnwr1v)
23133 continue
      qcpiaj7f = 3.0d0 + (pjb6wfoq-qcpiaj7f)**2
23130 continue
23126 continue
      return
      end
      subroutine gt9iulbf(he7mqnvy,ghz9vuba,w,gkdx5jal, rvy1fpli,
     &kuzxj1lo, bhcji9glto,zvau2lct,f6lsuzax,fvh2rwtc,dcfir2no)
      implicit logical (a-z)
      integer rvy1fpli,kuzxj1lo
      double precision he7mqnvy(rvy1fpli),ghz9vuba(rvy1fpli),w(rvy1fpli)
     &,gkdx5jal(kuzxj1lo+4), bhcji9glto(kuzxj1lo), zvau2lct(kuzxj1lo),
     &f6lsuzax(kuzxj1lo),fvh2rwtc(kuzxj1lo),dcfir2no(kuzxj1lo)
      double precision qaltf0nz,g9fvdrbw(4,1),work(16)
      double precision w2svdbx3tk, wv2svdbx3tk
      integer yq6lorbx,ayfnwr1v,ilo,dqlr5bse,pqzfxw4i, nhnpt1zym1
      integer ifour4, hbsl0gto
      hbsl0gto = 1
      ifour4 = 4
      nhnpt1zym1 = kuzxj1lo + 1
      do 23135 ayfnwr1v = 1,kuzxj1lo 
      bhcji9glto(ayfnwr1v) = 0.0d0 
      zvau2lct(ayfnwr1v) = 0.0d0 
      f6lsuzax(ayfnwr1v) = 0.0d0
      fvh2rwtc(ayfnwr1v) = 0.0d0 
      dcfir2no(ayfnwr1v) = 0.0d0
23135 continue
      ilo = 1
      qaltf0nz = 0.1d-9
      do 23137 ayfnwr1v = 1,rvy1fpli 
      call vinterv(gkdx5jal(1), nhnpt1zym1 ,he7mqnvy(ayfnwr1v),dqlr5bse,
     &pqzfxw4i)
      if(.not.(pqzfxw4i .eq. 1))goto 23139
      if(.not.(he7mqnvy(ayfnwr1v) .le. (gkdx5jal(dqlr5bse)+qaltf0nz)))
     &goto 23141
      dqlr5bse = dqlr5bse-1
      goto 23142
23141 continue
      return
23142 continue
23139 continue
      call vbsplvd(gkdx5jal,ifour4,he7mqnvy(ayfnwr1v),dqlr5bse,work,
     &g9fvdrbw,hbsl0gto)
      yq6lorbx = dqlr5bse-4+1
      w2svdbx3tk = w(ayfnwr1v)**2
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(1,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*
     &ghz9vuba(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(1,
     &1)
      f6lsuzax(yq6lorbx) = f6lsuzax(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(2,
     &1)
      fvh2rwtc(yq6lorbx) = fvh2rwtc(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(3,
     &1)
      dcfir2no(yq6lorbx) = dcfir2no(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,
     &1)
      yq6lorbx = dqlr5bse-4+2
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(2,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*
     &ghz9vuba(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(2,
     &1)
      f6lsuzax(yq6lorbx) = f6lsuzax(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(3,
     &1)
      fvh2rwtc(yq6lorbx) = fvh2rwtc(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,
     &1)
      yq6lorbx = dqlr5bse-4+3
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(3,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*
     &ghz9vuba(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(3,
     &1)
      f6lsuzax(yq6lorbx) = f6lsuzax(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,
     &1)
      yq6lorbx = dqlr5bse
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(4,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*
     &ghz9vuba(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,
     &1)
23137 continue
      return
      end
