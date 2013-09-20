C Output from Public domain Ratfor, version 1.01
      subroutine dnaoqj0l(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk, ankcghz
     *2,coef,sz,ifys6woa, qcpiaj7f,wbkq9zyi,parms, scrtch, gp0xjetb,l3zp
     *bstu,e5knafcg,wep0oibc,fbd5yktj)
      implicit logical (a-z)
      integer kuzxj1lo, nk, gp0xjetb, l3zpbstu(3), e5knafcg, wep0oibc, f
     *bd5yktj
      double precision penalt, pjb6wfoq, xs(kuzxj1lo), ys(kuzxj1lo), ws(
     *kuzxj1lo), ankcghz2(nk+4), coef(nk), sz(kuzxj1lo), ifys6woa(kuzxj1
     *lo), qcpiaj7f, wbkq9zyi, parms(3), scrtch(*)
      call hbzuprs6(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk, ankcghz2,coef
     *,sz,ifys6woa, qcpiaj7f,l3zpbstu(1),wbkq9zyi,l3zpbstu(2), l3zpbstu(
     *3), parms(1),parms(2),parms(3), gp0xjetb, scrtch(1), scrtch(nk+1),
     *scrtch(2*nk+1),scrtch(3*nk+1),scrtch(4*nk+1), scrtch(5*nk+1),scrtc
     *h(6*nk+1),scrtch(7*nk+1),scrtch(8*nk+1), scrtch(9*nk+1),scrtch(9*n
     *k+e5knafcg*nk+1),scrtch(9*nk+2*e5knafcg*nk+1), e5knafcg,wep0oibc,f
     *bd5yktj)
      return
      end
      subroutine hbzuprs6(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk, ankcghz
     *2,coef,sz,ifys6woa, qcpiaj7f,icrit,i9mwnvqt,ispar, c5aesxku, mynl7
     *uaq,zustx4fw,tol, gp0xjetb, xwy, zvau2lct,f6lsuzax,fvh2rwtc,dcfir2
     *no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, buhyalv4,fulcp8wa,plj0trq
     *x, e5knafcg,wep0oibc,fbd5yktj)
      implicit logical (a-z)
      integer kuzxj1lo,nk, icrit,ispar, gp0xjetb, e5knafcg,wep0oibc,fbd5
     *yktj
      integer c5aesxku
      double precision penalt,pjb6wfoq,xs(kuzxj1lo),ys(kuzxj1lo),ws(kuzx
     *j1lo), ankcghz2(nk+4), coef(nk),sz(kuzxj1lo),ifys6woa(kuzxj1lo), q
     *cpiaj7f,i9mwnvqt,mynl7uaq,zustx4fw,tol, xwy(nk), zvau2lct(nk),f6ls
     *uzax(nk),fvh2rwtc(nk),dcfir2no(nk), xecbg0pf(nk),z4grbpiq(nk),d7gl
     *zhbj(nk),v2eydbxs(nk), buhyalv4(e5knafcg,nk),fulcp8wa(e5knafcg,nk)
     *,plj0trqx(wep0oibc,nk)
      double precision t1,t2,ratio, a,b,c,d,e,qaltf0nz,xm,p,q,r,tol1,tol
     *2,u,v,w, fu,fv,fw,fx,x, ax,bx
      integer ayfnwr1v, viter
      double precision yjpnro8d, hmayv1xt
      yjpnro8d = 8.0d88
      hmayv1xt = 0.0d0
      d = 0.5d0
      u = 0.5d0
      ratio = 0.5d0
      ayfnwr1v = 1
23000 if(.not.(ayfnwr1v .le. kuzxj1lo))goto 23002
      if(ws(ayfnwr1v).gt.0.0d0)then
      ws(ayfnwr1v) = dsqrt(ws(ayfnwr1v))
      endif
23001 ayfnwr1v = ayfnwr1v+1
      goto 23000
23002 continue
      if(gp0xjetb .eq. 0)then
      call zosq7hub(xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs,ankcghz2,nk)
      call gt9iulbf(xs,ys,ws,ankcghz2, kuzxj1lo,nk, xwy,zvau2lct,f6lsuza
     *x,fvh2rwtc,dcfir2no)
      t1 = 0.0d0 
      t2 = 0.0d0
      do23007 ayfnwr1v = 3,nk-3 
      t1 = t1 + zvau2lct(ayfnwr1v) 
23007 continue
23008 continue
      do23009 ayfnwr1v = 3,nk-3 
      t2 = t2 + xecbg0pf(ayfnwr1v) 
23009 continue
23010 continue
      ratio = t1/t2
      gp0xjetb = 1
      endif
      if(ispar .eq. 1)then
      call wmhctl9x(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk,icrit, ankcghz
     *2,coef,sz,ifys6woa,qcpiaj7f, i9mwnvqt, xwy, zvau2lct,f6lsuzax,fvh2
     *rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, buhyalv4,fulcp
     *8wa,plj0trqx,e5knafcg,wep0oibc,fbd5yktj)
      return
      endif
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
      call wmhctl9x(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk,icrit, ankcghz
     *2,coef,sz,ifys6woa,qcpiaj7f, i9mwnvqt, xwy, zvau2lct,f6lsuzax,fvh2
     *rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, buhyalv4,fulcp
     *8wa,plj0trqx,e5knafcg,wep0oibc,fbd5yktj)
      fx = qcpiaj7f
      fv = fx
      fw = fx
23013 if(fbd5yktj .eq. 0)then
      viter = viter + 1
      xm = 0.5d0*(a + b)
      tol1 = qaltf0nz*dabs(x) + tol/3.0d0
      tol2 = 2.0d0*tol1
      if((dabs(x - xm) .le. (tol2 - 0.5d0*(b - a))) .or. (viter .gt. c5a
     *esxku))then
      go to 90
      endif
      if((dabs(e) .le. tol1) .or. (fx .ge. yjpnro8d) .or. (fv .ge. yjpnr
     *o8d) .or. (fw .ge. yjpnro8d))then
      go to 40
      endif
      r = (x - w)*(fx - fv)
      q = (x - v)*(fx - fw)
      p = (x - v)*q - (x - w)*r
      q = 2.0d0 * (q - r)
      if(q .gt. 0.0d0)then
      p = -p
      endif
      q = dabs(q)
      r = e
      e = d
      if((dabs(p) .ge. dabs(0.5d0*q*r)) .or. (q .eq. 0.0d0))then
      go to 40
      endif
      if((p .le. q*(a - x)) .or. (p .ge. q*(b - x)))then
      go to 40
      endif
      d = p/q
      u = x + d
      if((u - a) .lt. tol2)then
      d = dsign(tol1, xm - x)
      endif
      if((b - u) .lt. tol2)then
      d = dsign(tol1, xm - x)
      endif
      go to 50
40    if(x .ge. xm)then
      e = a - x
      else
      e = b - x
      endif
      d = c*e
50    if(dabs(d) .ge. tol1)then
      u = x + d
      else
      u = x + dsign(tol1, d)
      endif
      i9mwnvqt = ratio * dexp((-2.0d0 + u*6.0) * dlog(16.0d0))
      call wmhctl9x(penalt,pjb6wfoq,xs,ys,ws, kuzxj1lo,nk,icrit, ankcghz
     *2,coef,sz,ifys6woa,qcpiaj7f, i9mwnvqt, xwy, zvau2lct,f6lsuzax,fvh2
     *rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, buhyalv4,fulcp
     *8wa,plj0trqx,e5knafcg,wep0oibc,fbd5yktj)
      fu = qcpiaj7f
      if(fu .gt. yjpnro8d)then
      fu = 2.0d0 * yjpnro8d
      endif
      if(fu .le. fx)then
      if(u .ge. x)then
      a = x
      else
      b = x
      endif
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
      else
      if(u .lt. x)then
      a = u
      else
      b = u
      endif
      if((fu .le. fw) .or. (w .eq. x))then
      v = w
      fv = fw
      w = u
      fw = fu
      else
      if((fu .le. fv) .or. (v .eq. x) .or. (v .eq. w))then
      v = u
      fv = fu
      endif
      endif
      endif
      goto 23013
      endif
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
      double precision xecbg0pf(nb),z4grbpiq(nb),d7glzhbj(nb),v2eydbxs(n
     *b),tb(nb+4)
      integer dqlr5bse,ilo,pqzfxw4i, three3, ifour4, nbp1
      integer ayfnwr1v,iii,yq6lorbx
      integer i2svdbx3tk
      double precision g9fvdrbw(4,3),work(16),yw1(4),yw2(4), wpt
      double precision othird
      othird = 1.0d0 / 3.0d0
      three3 = 3
      ifour4 = 4
      nbp1 = nb + 1
      do23045 ayfnwr1v = 1,nb 
      xecbg0pf(ayfnwr1v) = 0.0d0
      z4grbpiq(ayfnwr1v) = 0.0d0
      d7glzhbj(ayfnwr1v) = 0.0d0
      v2eydbxs(ayfnwr1v) = 0.0d0 
23045 continue
23046 continue
      ilo = 1
      do23047 ayfnwr1v = 1,nb 
      call vinterv(tb(1), nbp1 ,tb(ayfnwr1v),dqlr5bse,pqzfxw4i)
      call vbsplvd(tb,ifour4,tb(ayfnwr1v),dqlr5bse,work,g9fvdrbw,three3)
      do23049 iii = 1,4 
      yw1(iii) = g9fvdrbw(iii,3) 
23049 continue
23050 continue
      call vbsplvd(tb,ifour4,tb(ayfnwr1v+1),dqlr5bse,work,g9fvdrbw,three
     *3)
      do23051 iii = 1,4 
      yw2(iii) = g9fvdrbw(iii,3) - yw1(iii) 
23051 continue
23052 continue
      wpt = tb(ayfnwr1v+1) - tb(ayfnwr1v)
      if(dqlr5bse .ge. 4)then
      do23055 iii = 1,4 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-4+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt * (yw1(iii)*yw1(
     *yq6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50
     * + yw2(iii)*yw2(yq6lorbx)*othird)
      yq6lorbx = iii+1
      if(yq6lorbx .le. 4)then
      z4grbpiq(i2svdbx3tk) = z4grbpiq(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      endif
      yq6lorbx = iii+2
      if(yq6lorbx .le. 4)then
      d7glzhbj(i2svdbx3tk) = d7glzhbj(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      endif
      yq6lorbx = iii+3
      if(yq6lorbx .le. 4)then
      v2eydbxs(i2svdbx3tk) = v2eydbxs(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      endif
23055 continue
23056 continue
      else
      if(dqlr5bse .eq. 3)then
      do23065 iii = 1,3 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-3+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      yq6lorbx = iii+1
      if(yq6lorbx .le. 3)then
      z4grbpiq(i2svdbx3tk) = z4grbpiq(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      endif
      yq6lorbx = iii+2
      if(yq6lorbx .le. 3)then
      d7glzhbj(i2svdbx3tk) = d7glzhbj(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      endif
23065 continue
23066 continue
      else
      if(dqlr5bse .eq. 2)then
      do23073 iii = 1,2 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-2+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      yq6lorbx = iii+1
      if(yq6lorbx .le. 2)then
      z4grbpiq(i2svdbx3tk) = z4grbpiq(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
      endif
23073 continue
23074 continue
      else
      if(dqlr5bse .eq. 1)then
      do23079 iii = 1,1 
      yq6lorbx = iii
      i2svdbx3tk = dqlr5bse-1+iii
      xecbg0pf(i2svdbx3tk) = xecbg0pf(i2svdbx3tk) + wpt* (yw1(iii)*yw1(y
     *q6lorbx) + (yw2(iii)*yw1(yq6lorbx) + yw2(yq6lorbx)*yw1(iii))*0.50 
     *+ yw2(iii)*yw2(yq6lorbx)*othird)
23079 continue
23080 continue
      endif
      endif
      endif
      endif
23047 continue
23048 continue
      return
      end
      subroutine vmnweiy2(buhyalv4,fulcp8wa,plj0trqx, e5knafcg,nk,wep0oi
     *bc,iflag)
      implicit logical (a-z)
      integer e5knafcg,nk,wep0oibc,iflag
      double precision buhyalv4(e5knafcg,nk), fulcp8wa(e5knafcg,nk), plj
     *0trqx(wep0oibc,nk)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh
      double precision wjm3(3),wjm2(2),wjm1(1),c0,c1,c2,c3
      double precision pcsuow9k, qdbgu6oi, upwkh5xz, rul5fnyd, ueydbrg6,
     * plce2srm, k3yvomnh, bfdjhu7l, ctfvwdu0
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      wjm3(1) = 0.0d0
      wjm3(2) = 0.0d0
      wjm3(3) = 0.0d0
      wjm2(1) = 0.0d0
      wjm2(2) = 0.0d0
      wjm1(1) = 0.0d0
      do23081 ayfnwr1v = 1,nk 
      yq6lorbx = nk-ayfnwr1v+1
      c0 = 1.0d0 / buhyalv4(4,yq6lorbx)
      if(yq6lorbx .le. (nk-3))then
      c1 = buhyalv4(1,yq6lorbx+3)*c0
      c2 = buhyalv4(2,yq6lorbx+2)*c0
      c3 = buhyalv4(3,yq6lorbx+1)*c0
      else
      if(yq6lorbx .eq. (nk-2))then
      c1 = 0.0d0
      c2 = buhyalv4(2,yq6lorbx+2)*c0
      c3 = buhyalv4(3,yq6lorbx+1)*c0
      else
      if(yq6lorbx .eq. (nk-1))then
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = buhyalv4(3,yq6lorbx+1)*c0
      else
      if(yq6lorbx .eq. nk)then
      c1 = 0.0d0
      c2 = 0.0d0
      c3 = 0.0d0
      endif
      endif
      endif
      endif
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
      fulcp8wa(4,yq6lorbx) = c0**2 + c1*(pcsuow9k + 2.0d0*(qdbgu6oi + up
     *wkh5xz)) + c2*(ueydbrg6 + 2.0d0* plce2srm) + c3*ctfvwdu0
      wjm3(1) = wjm2(1)
      wjm3(2) = wjm2(2)
      wjm3(3) = fulcp8wa(2,yq6lorbx)
      wjm2(1) = wjm1(1)
      wjm2(2) = fulcp8wa(3,yq6lorbx)
      wjm1(1) = fulcp8wa(4,yq6lorbx)
23081 continue
23082 continue
      if(iflag .eq. 0)then
      return
      endif
      do23093 ayfnwr1v = 1,nk 
      yq6lorbx = nk-ayfnwr1v+1
      gp1jxzuh = 1
23095 if(.not.(gp1jxzuh .le. 4 .and. yq6lorbx+gp1jxzuh-1 .le. nk))goto 2
     *3097
      plj0trqx(yq6lorbx,yq6lorbx+gp1jxzuh-1) = fulcp8wa(5-gp1jxzuh,yq6lo
     *rbx)
23096 gp1jxzuh = gp1jxzuh+1
      goto 23095
23097 continue
23093 continue
23094 continue
      do23098 ayfnwr1v = 1,nk 
      yq6lorbx = nk-ayfnwr1v+1
      gp1jxzuh = yq6lorbx-4
23100 if(.not.(gp1jxzuh .ge. 1))goto 23102
      c0 = 1.0 / buhyalv4(4,gp1jxzuh) 
      c1 = buhyalv4(1,gp1jxzuh+3)*c0
      c2 = buhyalv4(2,gp1jxzuh+2)*c0 
      c3 = buhyalv4(3,gp1jxzuh+1)*c0
      plj0trqx(gp1jxzuh,yq6lorbx) = 0.0d0- ( c1*plj0trqx(gp1jxzuh+3,yq6l
     *orbx) + c2*plj0trqx(gp1jxzuh+2,yq6lorbx) + c3*plj0trqx(gp1jxzuh+1,
     *yq6lorbx) )
23101 gp1jxzuh = gp1jxzuh-1
      goto 23100
23102 continue
23098 continue
23099 continue
      return
      end
      subroutine wmhctl9x(penalt,pjb6wfoq,x,y,w, kuzxj1lo,nk,icrit, ankc
     *ghz2,coef,sz,ifys6woa, qcpiaj7f, i9mwnvqt, xwy, zvau2lct,f6lsuzax,
     *fvh2rwtc,dcfir2no, xecbg0pf,z4grbpiq,d7glzhbj,v2eydbxs, buhyalv4,f
     *ulcp8wa,plj0trqx, e5knafcg,wep0oibc,info)
      implicit logical (a-z)
      integer kuzxj1lo,nk,icrit, e5knafcg,wep0oibc,info
      double precision penalt,pjb6wfoq,x(kuzxj1lo),y(kuzxj1lo),w(kuzxj1l
     *o)
      double precision ankcghz2(nk+4), coef(nk),sz(kuzxj1lo),ifys6woa(ku
     *zxj1lo), qcpiaj7f, i9mwnvqt, xwy(nk)
      double precision zvau2lct(nk),f6lsuzax(nk),fvh2rwtc(nk),dcfir2no(n
     *k)
      double precision xecbg0pf(nk),z4grbpiq(nk),d7glzhbj(nk),v2eydbxs(n
     *k), buhyalv4(e5knafcg,nk),fulcp8wa(e5knafcg,nk),plj0trqx(wep0oibc,
     *nk)
      double precision resss, work(16), b0,b1,b2,b3,qaltf0nz, g9fvdrbw(4
     *,1), xv,eqdf
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
      do23103 ayfnwr1v = 1,nk 
      coef(ayfnwr1v) = xwy(ayfnwr1v) 
23103 continue
23104 continue
      do23105 ayfnwr1v = 1,nk 
      buhyalv4(4,ayfnwr1v) = zvau2lct(ayfnwr1v)+i9mwnvqt*xecbg0pf(ayfnwr
     *1v) 
23105 continue
23106 continue
      do23107 ayfnwr1v = 1,(nk-1) 
      buhyalv4(3,ayfnwr1v+1) = f6lsuzax(ayfnwr1v)+i9mwnvqt*z4grbpiq(ayfn
     *wr1v) 
23107 continue
23108 continue
      do23109 ayfnwr1v = 1,(nk-2) 
      buhyalv4(2,ayfnwr1v+2) = fvh2rwtc(ayfnwr1v)+i9mwnvqt*d7glzhbj(ayfn
     *wr1v) 
23109 continue
23110 continue
      do23111 ayfnwr1v = 1,(nk-3) 
      buhyalv4(1,ayfnwr1v+3) = dcfir2no(ayfnwr1v)+i9mwnvqt*v2eydbxs(ayfn
     *wr1v) 
23111 continue
23112 continue
      call dpbfa8(buhyalv4,e5knafcg,nk,three3,info)
      if(info .ne. 0)then
      return
      endif
      call dpbsl8(buhyalv4,e5knafcg,nk,three3,coef)
      icoef = 1
      do23115 ayfnwr1v = 1,kuzxj1lo 
      xv = x(ayfnwr1v)
      call wbvalue(ankcghz2,coef, nk,ifour4,xv,izero0, sz(ayfnwr1v))
23115 continue
23116 continue
      if(icrit .eq. 0)then
      return
      endif
      call vmnweiy2(buhyalv4,fulcp8wa,plj0trqx, e5knafcg,nk,wep0oibc,ize
     *ro0)
      do23119 ayfnwr1v = 1,kuzxj1lo 
      xv = x(ayfnwr1v)
      call vinterv(ankcghz2(1), nkp1 ,xv,dqlr5bse,pqzfxw4i)
      if(pqzfxw4i .eq. -1)then
      dqlr5bse = 4 
      xv = ankcghz2(4) + qaltf0nz 
      endif
      if(pqzfxw4i .eq. 1)then
      dqlr5bse = nk 
      xv = ankcghz2(nk+1) - qaltf0nz 
      endif
      yq6lorbx = dqlr5bse-3
      call vbsplvd(ankcghz2,ifour4,xv,dqlr5bse,work,g9fvdrbw,hbsl0gto)
      b0 = g9fvdrbw(1,1)
      b1 = g9fvdrbw(2,1)
      b2 = g9fvdrbw(3,1)
      b3 = g9fvdrbw(4,1)
      qtce8hzo = (b0 *(fulcp8wa(4,yq6lorbx)*b0 + 2.0d0*(fulcp8wa(3,yq6lo
     *rbx)*b1 + fulcp8wa(2,yq6lorbx)*b2 + fulcp8wa(1,yq6lorbx)*b3)) + b1
     * *(fulcp8wa(4,yq6lorbx+1)*b1 + 2.0d0*(fulcp8wa(3,yq6lorbx+1)*b2 + 
     *fulcp8wa(2,yq6lorbx+1)*b3)) + b2 *(fulcp8wa(4,yq6lorbx+2)*b2 + 2.0
     *d0* fulcp8wa(3,yq6lorbx+2)*b3 )+ b3**2* fulcp8wa(4,yq6lorbx+3)) * 
     *w(ayfnwr1v)**2
      ifys6woa(ayfnwr1v) = qtce8hzo
23119 continue
23120 continue
      if(icrit .eq. 1)then
      resss = 0.0d0 
      eqdf = 0.0d0 
      rxeqjn0y = 0.0d0
      do23127 ayfnwr1v = 1,kuzxj1lo 
      resss = resss + ((y(ayfnwr1v)-sz(ayfnwr1v))*w(ayfnwr1v))**2
      eqdf = eqdf + ifys6woa(ayfnwr1v)
      rxeqjn0y = rxeqjn0y + w(ayfnwr1v)*w(ayfnwr1v)
23127 continue
23128 continue
      qcpiaj7f = (resss/rxeqjn0y)/((1.0d0-(pjb6wfoq+penalt*eqdf)/rxeqjn0
     *y)**2)
      else
      if(icrit .eq. 2)then
      qcpiaj7f = 0.0d0
      rxeqjn0y = 0.0d0
      do23131 ayfnwr1v = 1,kuzxj1lo 
      qcpiaj7f = qcpiaj7f + (((y(ayfnwr1v)-sz(ayfnwr1v))*w(ayfnwr1v))/(1
     *.0d0-ifys6woa(ayfnwr1v)))**2
      rxeqjn0y = rxeqjn0y + w(ayfnwr1v)*w(ayfnwr1v)
23131 continue
23132 continue
      qcpiaj7f = qcpiaj7f / rxeqjn0y
      else
      qcpiaj7f = 0.0d0
      do23133 ayfnwr1v = 1,kuzxj1lo 
      qcpiaj7f = qcpiaj7f+ifys6woa(ayfnwr1v)
23133 continue
23134 continue
      qcpiaj7f = 3.0d0 + (pjb6wfoq-qcpiaj7f)**2
      endif
      endif
      return
      end
      subroutine gt9iulbf(he7mqnvy,ghz9vuba,w,gkdx5jal, rvy1fpli,kuzxj1l
     *o, bhcji9glto,zvau2lct,f6lsuzax,fvh2rwtc,dcfir2no)
      implicit logical (a-z)
      integer rvy1fpli,kuzxj1lo
      double precision he7mqnvy(rvy1fpli),ghz9vuba(rvy1fpli),w(rvy1fpli)
     *,gkdx5jal(kuzxj1lo+4), bhcji9glto(kuzxj1lo), zvau2lct(kuzxj1lo),f6
     *lsuzax(kuzxj1lo),fvh2rwtc(kuzxj1lo),dcfir2no(kuzxj1lo)
      double precision qaltf0nz,g9fvdrbw(4,1),work(16)
      double precision w2svdbx3tk, wv2svdbx3tk
      integer yq6lorbx,ayfnwr1v,ilo,dqlr5bse,pqzfxw4i, nhnpt1zym1
      integer ifour4, hbsl0gto
      hbsl0gto = 1
      ifour4 = 4
      nhnpt1zym1 = kuzxj1lo + 1
      do23135 ayfnwr1v = 1,kuzxj1lo 
      bhcji9glto(ayfnwr1v) = 0.0d0 
      zvau2lct(ayfnwr1v) = 0.0d0 
      f6lsuzax(ayfnwr1v) = 0.0d0
      fvh2rwtc(ayfnwr1v) = 0.0d0 
      dcfir2no(ayfnwr1v) = 0.0d0
23135 continue
23136 continue
      ilo = 1
      qaltf0nz = 0.1d-9
      do23137 ayfnwr1v = 1,rvy1fpli 
      call vinterv(gkdx5jal(1), nhnpt1zym1 ,he7mqnvy(ayfnwr1v),dqlr5bse,
     *pqzfxw4i)
      if(pqzfxw4i .eq. 1)then
      if(he7mqnvy(ayfnwr1v) .le. (gkdx5jal(dqlr5bse)+qaltf0nz))then
      dqlr5bse = dqlr5bse-1
      else
      return
      endif
      endif
      call vbsplvd(gkdx5jal,ifour4,he7mqnvy(ayfnwr1v),dqlr5bse,work,g9fv
     *drbw,hbsl0gto)
      yq6lorbx = dqlr5bse-4+1
      w2svdbx3tk = w(ayfnwr1v)**2
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(1,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*ghz9vuba
     *(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(1,1
     *)
      f6lsuzax(yq6lorbx) = f6lsuzax(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(2,1
     *)
      fvh2rwtc(yq6lorbx) = fvh2rwtc(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(3,1
     *)
      dcfir2no(yq6lorbx) = dcfir2no(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,1
     *)
      yq6lorbx = dqlr5bse-4+2
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(2,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*ghz9vuba
     *(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(2,1
     *)
      f6lsuzax(yq6lorbx) = f6lsuzax(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(3,1
     *)
      fvh2rwtc(yq6lorbx) = fvh2rwtc(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,1
     *)
      yq6lorbx = dqlr5bse-4+3
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(3,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*ghz9vuba
     *(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(3,1
     *)
      f6lsuzax(yq6lorbx) = f6lsuzax(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,1
     *)
      yq6lorbx = dqlr5bse
      wv2svdbx3tk = w2svdbx3tk * g9fvdrbw(4,1)
      bhcji9glto(yq6lorbx) = bhcji9glto(yq6lorbx) + wv2svdbx3tk*ghz9vuba
     *(ayfnwr1v)
      zvau2lct(yq6lorbx) = zvau2lct(yq6lorbx) + wv2svdbx3tk*g9fvdrbw(4,1
     *)
23137 continue
23138 continue
      return
      end
