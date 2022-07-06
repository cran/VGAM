C Output from Public domain Ratfor, version 1.01
      subroutine vbvs(kuzxj1lo,ankcghz2,rpyis2kc,nk,he7mqnvy,smat,order,
     *wy1vqfzu)
      integer kuzxj1lo, nk, order, wy1vqfzu
      double precision ankcghz2(nk+4), rpyis2kc(nk,wy1vqfzu), he7mqnvy(k
     *uzxj1lo), smat(kuzxj1lo,wy1vqfzu)
      double precision chw8lzty
      integer ayfnwr1v, yq6lorbx, ifour4
      ifour4 = 4
      do23000 yq6lorbx=1,wy1vqfzu 
      do23002 ayfnwr1v=1,kuzxj1lo 
      chw8lzty = he7mqnvy(ayfnwr1v)
      call wbvalue(ankcghz2, rpyis2kc(1,yq6lorbx), nk, ifour4, chw8lzty,
     * order, smat(ayfnwr1v,yq6lorbx))
23002 continue
C23003 continue
23000 continue
C23001 continue
      return
      end
      subroutine tfeswo7c(osiz4fxy, nk, wy1vqfzu, ldk, wbkq9zyi, sgmat)
      implicit logical (a-z)
      integer nk, wy1vqfzu, ldk
      double precision osiz4fxy(ldk,nk*wy1vqfzu), wbkq9zyi(wy1vqfzu), sg
     *mat(nk,4)
      integer ayfnwr1v, yq6lorbx
      do23004 ayfnwr1v=1,nk 
      do23006 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk,(ayfnwr1v-1)*wy1vqfzu+yq6lorbx) = osiz4fxy(ldk,(ayfnw
     *r1v-1)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorbx) * sgmat(ayfnwr1v,1)
23006 continue
C23007 continue
23004 continue
C23005 continue
      do23008 ayfnwr1v=1,(nk-1) 
      do23010 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk-wy1vqfzu,(ayfnwr1v-0)*wy1vqfzu+yq6lorbx) = osiz4fxy(l
     *dk-wy1vqfzu,(ayfnwr1v-0)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorbx) *
     * sgmat(ayfnwr1v,2)
23010 continue
C23011 continue
23008 continue
C23009 continue
      do23012 ayfnwr1v=1,(nk-2) 
      do23014 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk-2*wy1vqfzu,(ayfnwr1v+1)*wy1vqfzu+yq6lorbx) = osiz4fxy
     *(ldk-2*wy1vqfzu,(ayfnwr1v+1)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorb
     *x) * sgmat(ayfnwr1v,3)
23014 continue
C23015 continue
23012 continue
C23013 continue
      do23016 ayfnwr1v=1,(nk-3) 
      do23018 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk-3*wy1vqfzu,(ayfnwr1v+2)*wy1vqfzu+yq6lorbx) = osiz4fxy
     *(ldk-3*wy1vqfzu,(ayfnwr1v+2)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorb
     *x) * sgmat(ayfnwr1v,4)
23018 continue
C23019 continue
23016 continue
C23017 continue
      return
      end
      subroutine ybnagt8k(iii, cz8qdfyj, tesdm5kv, g9fvdrbw, osiz4fxy, w
     *mat, kxvq6sfw, nyfu9rod, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxd
     *w1, dufozmt7)
      implicit logical (a-z)
      integer iii, cz8qdfyj, tesdm5kv, kxvq6sfw, nyfu9rod, wy1vqfzu, ldk
     *, dimw, kuzxj1lo, nk, tgiyxdw1(*), dufozmt7(*)
      double precision g9fvdrbw(4,*), osiz4fxy(ldk, nk*wy1vqfzu), wmat(k
     *uzxj1lo,dimw)
      double precision obr6tcex
      integer urohxe6t, nead, bcol, brow, biuvowq2, nbj8tdsk
      bcol = cz8qdfyj + tesdm5kv
      brow = cz8qdfyj
      do23020 urohxe6t=1,dimw 
      obr6tcex = wmat(iii,urohxe6t) * g9fvdrbw(kxvq6sfw,1) * g9fvdrbw(ny
     *fu9rod,1)
      biuvowq2 = (brow-1)*wy1vqfzu + tgiyxdw1(urohxe6t)
      nbj8tdsk = (bcol-1)*wy1vqfzu + dufozmt7(urohxe6t)
      nead = nbj8tdsk - biuvowq2
      osiz4fxy(ldk-nead, nbj8tdsk) = osiz4fxy(ldk-nead, nbj8tdsk) + obr6
     *tcex
      if(tesdm5kv .gt. 0 .and. dufozmt7(urohxe6t) .ne. tgiyxdw1(urohxe6t
     *))then
      biuvowq2 = (brow-1)*wy1vqfzu + dufozmt7(urohxe6t)
      nbj8tdsk = (bcol-1)*wy1vqfzu + tgiyxdw1(urohxe6t)
      nead = nbj8tdsk - biuvowq2
      osiz4fxy(ldk-nead, nbj8tdsk) = osiz4fxy(ldk-nead, nbj8tdsk) + obr6
     *tcex
      endif
23020 continue
C23021 continue
      return
      end
      subroutine vsplin(he7mqnvy,rbne6ouj,wmat,kuzxj1lo,gkdx5jal, nk,ldk
     *,wy1vqfzu,dimw, tgiyxdw1,dufozmt7, wkmm, wbkq9zyi, info, t8hwvalr,
     * rpyis2kc, osiz4fxy, btwy, sgdub, ui8ysltq, yzoe1rsp, bmb, ifys6wo
     *a, dof, scrtch, fbd5yktj, truen)
      implicit logical (a-z)
      integer kuzxj1lo, nk, ldk, wy1vqfzu, dimw, tgiyxdw1(*), dufozmt7(*
     *), info, fbd5yktj, truen
      integer yzoe1rsp
      double precision he7mqnvy(kuzxj1lo), rbne6ouj(kuzxj1lo,wy1vqfzu), 
     *wmat(kuzxj1lo,dimw), gkdx5jal(nk+4), wkmm(wy1vqfzu,wy1vqfzu,16), w
     *bkq9zyi(wy1vqfzu), t8hwvalr(kuzxj1lo,wy1vqfzu), rpyis2kc(nk,wy1vqf
     *zu), osiz4fxy(ldk,nk*wy1vqfzu), btwy(wy1vqfzu,nk)
      double precision sgdub(nk,wy1vqfzu), ui8ysltq(truen,wy1vqfzu), bmb
     *(wy1vqfzu,wy1vqfzu), ifys6woa(kuzxj1lo,wy1vqfzu), dof(wy1vqfzu), s
     *crtch(*)
      integer yq6lorbx, ayfnwr1v, dqlr5bse, pqzfxw4i, urohxe6t, icrit
      integer gp0xjetb, e5knafcg, wep0oibc, l3zpbstu(3), ispar, i1loc
      double precision qaltf0nz, g9fvdrbw(4,1), ms0qypiw(16), penalt, qc
     *piaj7f, fp6nozvx, waiez6nt, toldf, parms(3)
      do23024 yq6lorbx=1,wy1vqfzu 
      if(wbkq9zyi(yq6lorbx) .eq. 0.0d0)then
      ispar=0
      icrit=3
      else
      ispar=1
      icrit=1
      endif
      if((wy1vqfzu .eq. 1) .or. (dimw.eq.wy1vqfzu) .or. (ispar .eq. 0))t
     *hen
      e5knafcg = 4
      fp6nozvx = 1.50d0
      waiez6nt = 0.00d0
      wep0oibc = 1
      toldf=0.001d0
      if(wy1vqfzu.eq.1)then
      toldf=0.005d0
      else
      if(wy1vqfzu.eq.2)then
      toldf=0.015d0
      else
      if(wy1vqfzu.eq.3)then
      toldf=0.025d0
      else
      toldf=0.045d0
      endif
      endif
      endif
      l3zpbstu(1) = icrit
      l3zpbstu(2) = ispar
      l3zpbstu(3) = 300
      parms(1) = waiez6nt
      parms(2) = fp6nozvx
      parms(3) = toldf
      gp0xjetb=0
      if((wy1vqfzu .eq. 1) .or. (dimw.eq.wy1vqfzu))then
      do23038 ayfnwr1v=1,kuzxj1lo 
      rbne6ouj(ayfnwr1v,yq6lorbx) = rbne6ouj(ayfnwr1v,yq6lorbx) / wmat(a
     *yfnwr1v,yq6lorbx)
23038 continue
C23039 continue
      call dnaoqj0l(penalt, dof(yq6lorbx), he7mqnvy, rbne6ouj(1,yq6lorbx
     *), wmat(1,yq6lorbx), kuzxj1lo,nk, gkdx5jal,rpyis2kc(1,yq6lorbx), t
     *8hwvalr(1,yq6lorbx), ifys6woa(1,yq6lorbx), qcpiaj7f,wbkq9zyi(yq6lo
     *rbx),parms, scrtch, gp0xjetb,l3zpbstu, e5knafcg,wep0oibc,fbd5yktj)
      if(fbd5yktj .ne. 0)then
      return
      endif
      do23042 ayfnwr1v=1,kuzxj1lo 
      wmat(ayfnwr1v,yq6lorbx) = wmat(ayfnwr1v,yq6lorbx) * wmat(ayfnwr1v,
     *yq6lorbx)
23042 continue
C23043 continue
      if(yzoe1rsp .ne. 0)then
      do23046 ayfnwr1v=1,kuzxj1lo 
      ui8ysltq(ayfnwr1v,yq6lorbx) = ifys6woa(ayfnwr1v,yq6lorbx) / wmat(a
     *yfnwr1v,yq6lorbx)
23046 continue
C23047 continue
      endif
      else
      call dnaoqj0l(penalt, dof(yq6lorbx), he7mqnvy, btwy(1,yq6lorbx), w
     *mat(1,yq6lorbx), kuzxj1lo,nk, gkdx5jal,rpyis2kc(1,yq6lorbx),t8hwva
     *lr(1,yq6lorbx), ifys6woa(1,yq6lorbx), qcpiaj7f,wbkq9zyi(yq6lorbx),
     *parms, scrtch, gp0xjetb,l3zpbstu, e5knafcg,wep0oibc,fbd5yktj)
      if(fbd5yktj .ne. 0)then
      return
      endif
      do23050 ayfnwr1v=1,kuzxj1lo 
      wmat(ayfnwr1v,yq6lorbx) = wmat(ayfnwr1v,yq6lorbx) * wmat(ayfnwr1v,
     *yq6lorbx)
23050 continue
C23051 continue
      endif
      if(fbd5yktj .ne. 0)then
      return
      endif
      endif
23024 continue
C23025 continue
      if((wy1vqfzu .eq. 1) .or. (dimw .eq. wy1vqfzu))then
      return
      endif
      do23056 ayfnwr1v=1,nk 
      do23058 yq6lorbx=1,wy1vqfzu 
      btwy(yq6lorbx,ayfnwr1v)=0.0d0
23058 continue
C23059 continue
23056 continue
C23057 continue
      do23060 ayfnwr1v=1,(nk*wy1vqfzu) 
      do23062 yq6lorbx=1,ldk 
      osiz4fxy(yq6lorbx,ayfnwr1v) = 0.0d0
23062 continue
C23063 continue
23060 continue
C23061 continue
      qaltf0nz = 0.1d-9
      do23064 ayfnwr1v=1,kuzxj1lo 
      call vinterv(gkdx5jal(1),(nk+1),he7mqnvy(ayfnwr1v),dqlr5bse,pqzfxw
     *4i)
      if(pqzfxw4i .eq. 1)then
      if(he7mqnvy(ayfnwr1v) .le. (gkdx5jal(dqlr5bse)+qaltf0nz))then
      dqlr5bse=dqlr5bse-1
      else
      return
      endif
      endif
      call vbsplvd(gkdx5jal,4,he7mqnvy(ayfnwr1v),dqlr5bse,ms0qypiw,g9fvd
     *rbw,1)
      yq6lorbx= dqlr5bse-4+1
      do23070 urohxe6t=1,wy1vqfzu 
      btwy(urohxe6t,yq6lorbx)=btwy(urohxe6t,yq6lorbx) + rbne6ouj(ayfnwr1
     *v,urohxe6t) * g9fvdrbw(1,1)
23070 continue
C23071 continue
      call ybnagt8k(ayfnwr1v, yq6lorbx, 0, g9fvdrbw, osiz4fxy, wmat, 1, 
     *1, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 1, g9fvdrbw, osiz4fxy, wmat, 1, 
     *2, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 2, g9fvdrbw, osiz4fxy, wmat, 1, 
     *3, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 3, g9fvdrbw, osiz4fxy, wmat, 1, 
     *4, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      yq6lorbx= dqlr5bse-4+2
      do23072 urohxe6t=1,wy1vqfzu 
      btwy(urohxe6t,yq6lorbx)=btwy(urohxe6t,yq6lorbx) + rbne6ouj(ayfnwr1
     *v,urohxe6t) * g9fvdrbw(2,1)
23072 continue
C23073 continue
      call ybnagt8k(ayfnwr1v, yq6lorbx, 0, g9fvdrbw, osiz4fxy, wmat, 2, 
     *2, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 1, g9fvdrbw, osiz4fxy, wmat, 2, 
     *3, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 2, g9fvdrbw, osiz4fxy, wmat, 2, 
     *4, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      yq6lorbx= dqlr5bse-4+3
      do23074 urohxe6t=1,wy1vqfzu 
      btwy(urohxe6t,yq6lorbx)=btwy(urohxe6t,yq6lorbx) + rbne6ouj(ayfnwr1
     *v,urohxe6t) * g9fvdrbw(3,1)
23074 continue
C23075 continue
      call ybnagt8k(ayfnwr1v, yq6lorbx, 0, g9fvdrbw, osiz4fxy, wmat, 3, 
     *3, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 1, g9fvdrbw, osiz4fxy, wmat, 3, 
     *4, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      yq6lorbx= dqlr5bse-4+4
      do23076 urohxe6t=1,wy1vqfzu 
      btwy(urohxe6t,yq6lorbx)=btwy(urohxe6t,yq6lorbx) + rbne6ouj(ayfnwr1
     *v,urohxe6t) * g9fvdrbw(4,1)
23076 continue
C23077 continue
      call ybnagt8k(ayfnwr1v, yq6lorbx, 0, g9fvdrbw, osiz4fxy, wmat, 4, 
     *4, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
23064 continue
C23065 continue
      call zosq7hub(sgdub(1,1), sgdub(1,2), sgdub(1,3), sgdub(1,4), gkdx
     *5jal, nk)
      call tfeswo7c(osiz4fxy, nk, wy1vqfzu, ldk, wbkq9zyi, sgdub)
      call vdpbfa7(osiz4fxy, ldk, nk*wy1vqfzu, ldk-1, info, sgdub)
      if(info .ne. 0)then
      return
      endif
      call vdpbsl7(osiz4fxy, ldk, nk*wy1vqfzu, ldk-1, btwy, sgdub)
      i1loc = 0
      do23080 ayfnwr1v=1,nk 
      do23082 yq6lorbx=1,wy1vqfzu 
      i1loc = i1loc + 1
      rpyis2kc(ayfnwr1v,yq6lorbx) = btwy(yq6lorbx,ayfnwr1v)
23082 continue
C23083 continue
23080 continue
C23081 continue
      call cn8kzpab(gkdx5jal, he7mqnvy, rpyis2kc, kuzxj1lo, nk, wy1vqfzu
     *, t8hwvalr)
      call vicb2(osiz4fxy, osiz4fxy, sgdub, wkmm, ldk-1, nk*wy1vqfzu)
      call icpd0omv(osiz4fxy, he7mqnvy, gkdx5jal, ui8ysltq, ldk, kuzxj1l
     *o, nk, wy1vqfzu, yzoe1rsp, bmb, wkmm, wmat, ifys6woa, dimw, tgiyxd
     *w1, dufozmt7, truen)
      return
      end
      subroutine cn8kzpab(ankcghz2, he7mqnvy, rpyis2kc, kuzxj1lo, nk, wy
     *1vqfzu, t8hwvalr)
      implicit logical (a-z)
      integer kuzxj1lo, nk, wy1vqfzu
      double precision ankcghz2(nk+4), he7mqnvy(kuzxj1lo), rpyis2kc(nk,w
     *y1vqfzu), t8hwvalr(kuzxj1lo,wy1vqfzu)
      double precision chw8lzty
      integer ayfnwr1v, yq6lorbx, izero0, ifour4
      izero0 = 0
      ifour4 = 4
      do23084 ayfnwr1v=1,kuzxj1lo 
      chw8lzty = he7mqnvy(ayfnwr1v)
      do23086 yq6lorbx=1,wy1vqfzu 
      call wbvalue(ankcghz2, rpyis2kc(1,yq6lorbx), nk, ifour4, chw8lzty,
     * izero0, t8hwvalr(ayfnwr1v,yq6lorbx))
23086 continue
C23087 continue
23084 continue
C23085 continue
      return
      end
      subroutine vsuff9(kuzxj1lo,nef,ezlgm2up, he7mqnvy,tlgduey8,wmat, p
     *ygsw6ko,pasjmo8g,wbar,uwbar,wpasjmo8g, wy1vqfzu, dimw, dimu, tgiyx
     *dw1, dufozmt7, work, work2, hjm2ktyr, kgwmz4ip, iz2nbfjc, wuwbar, 
     *dvhw1ulq)
      implicit logical (a-z)
      integer kuzxj1lo, nef, ezlgm2up(kuzxj1lo), wy1vqfzu, dimw, dimu, k
     *gwmz4ip, iz2nbfjc, wuwbar, dvhw1ulq, tgiyxdw1(*),dufozmt7(*)
      double precision he7mqnvy(kuzxj1lo), tlgduey8(kuzxj1lo,wy1vqfzu), 
     *wmat(kuzxj1lo,dimw), pygsw6ko(nef), pasjmo8g(nef,wy1vqfzu), wbar(n
     *ef,*), uwbar(dimu,nef), wpasjmo8g(nef,wy1vqfzu), work(wy1vqfzu,wy1
     *vqfzu+1), work2(kgwmz4ip,kgwmz4ip+1), hjm2ktyr(wy1vqfzu,kgwmz4ip)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, urohxe6t, bpvaqm5z, imk5wjxg
      integer oneint
      oneint = 1
      if(iz2nbfjc .eq. 1)then
      if((dimu .ne. dimw) .or. (kgwmz4ip .ne. wy1vqfzu))then
      dvhw1ulq = 0
      return
      endif
      endif
      imk5wjxg = wy1vqfzu * (wy1vqfzu+1) / 2
      if(dimw .gt. imk5wjxg)then
      endif
      call qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      do23094 ayfnwr1v=1,kuzxj1lo 
      pygsw6ko(ezlgm2up(ayfnwr1v))=he7mqnvy(ayfnwr1v)
23094 continue
C23095 continue
      do23096 yq6lorbx=1,wy1vqfzu 
      do23098 ayfnwr1v=1,nef 
      wpasjmo8g(ayfnwr1v,yq6lorbx) = 0.0d0
23098 continue
C23099 continue
23096 continue
C23097 continue
      do23100 yq6lorbx=1,dimw 
      do23102 ayfnwr1v=1,nef 
      wbar(ayfnwr1v,yq6lorbx) = 0.0d0
23102 continue
C23103 continue
23100 continue
C23101 continue
      if(dimw .ne. imk5wjxg)then
      do23106 gp1jxzuh=1,wy1vqfzu 
      do23108 yq6lorbx=1,wy1vqfzu 
      work(yq6lorbx,gp1jxzuh) = 0.0d0
23108 continue
C23109 continue
23106 continue
C23107 continue
      endif
      do23110 ayfnwr1v=1,kuzxj1lo 
      do23112 yq6lorbx=1,dimw 
      work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx)) = wmat(ayfnwr1v,yq6lor
     *bx)
      work(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx)) = work(tgiyxdw1(yq6lor
     *bx),dufozmt7(yq6lorbx))
23112 continue
C23113 continue
      do23114 yq6lorbx=1,wy1vqfzu 
      do23116 gp1jxzuh=1,wy1vqfzu 
      wpasjmo8g(ezlgm2up(ayfnwr1v),yq6lorbx) = wpasjmo8g(ezlgm2up(ayfnwr
     *1v),yq6lorbx) + work(yq6lorbx,gp1jxzuh)*tlgduey8(ayfnwr1v,gp1jxzuh
     *)
23116 continue
C23117 continue
23114 continue
C23115 continue
      do23118 yq6lorbx=1,dimw 
      wbar(ezlgm2up(ayfnwr1v),yq6lorbx) = wbar(ezlgm2up(ayfnwr1v),yq6lor
     *bx) + wmat(ayfnwr1v,yq6lorbx)
23118 continue
C23119 continue
23110 continue
C23111 continue
      dvhw1ulq = 1
      if(iz2nbfjc .eq. 1)then
      do23122 ayfnwr1v=1,nef 
      do23124 yq6lorbx=1,dimw 
      work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx)) = wbar(ayfnwr1v,yq6lor
     *bx)
      work(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx)) = work(tgiyxdw1(yq6lor
     *bx),dufozmt7(yq6lorbx))
23124 continue
C23125 continue
      do23126 yq6lorbx=1,wy1vqfzu 
      work(yq6lorbx,wy1vqfzu+1)=wpasjmo8g(ayfnwr1v,yq6lorbx)
23126 continue
C23127 continue
      call vcholf(work, work(1,wy1vqfzu+1), wy1vqfzu, dvhw1ulq, oneint)
      if(dvhw1ulq .ne. 1)then
      return
      endif
      if(wuwbar .ne. 0)then
      do23132 yq6lorbx=1,dimw 
      uwbar(yq6lorbx,ayfnwr1v) = work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lor
     *bx))
23132 continue
C23133 continue
      endif
      do23134 yq6lorbx=1,wy1vqfzu 
      pasjmo8g(ayfnwr1v,yq6lorbx)=work(yq6lorbx,wy1vqfzu+1)
23134 continue
C23135 continue
23122 continue
C23123 continue
      else
      if(dimw .ne. imk5wjxg)then
      do23138 yq6lorbx=1,wy1vqfzu 
      do23140 gp1jxzuh=1,wy1vqfzu 
      work(yq6lorbx,gp1jxzuh) = 0.0d0
23140 continue
C23141 continue
23138 continue
C23139 continue
      endif
      do23142 ayfnwr1v=1,nef 
      call qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      do23144 yq6lorbx=1,dimw 
      work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx)) = wbar(ayfnwr1v,yq6lor
     *bx)
      work(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx)) = work(tgiyxdw1(yq6lor
     *bx),dufozmt7(yq6lorbx))
23144 continue
C23145 continue
      do23146 yq6lorbx=1,wy1vqfzu 
      work(yq6lorbx,wy1vqfzu+1)=wpasjmo8g(ayfnwr1v,yq6lorbx)
23146 continue
C23147 continue
      do23148 yq6lorbx=1,kgwmz4ip 
      do23150 gp1jxzuh=yq6lorbx,kgwmz4ip 
      work2(yq6lorbx,gp1jxzuh) = 0.0d0
      do23152 urohxe6t=1,wy1vqfzu 
      do23154 bpvaqm5z=1,wy1vqfzu 
      work2(yq6lorbx,gp1jxzuh) = work2(yq6lorbx,gp1jxzuh) + hjm2ktyr(uro
     *hxe6t,yq6lorbx) * work(urohxe6t,bpvaqm5z) * hjm2ktyr(bpvaqm5z,gp1j
     *xzuh)
23154 continue
C23155 continue
23152 continue
C23153 continue
23150 continue
C23151 continue
23148 continue
C23149 continue
      call qpsedg8xf(tgiyxdw1, dufozmt7, kgwmz4ip)
      do23156 yq6lorbx=1,dimu 
      wbar(ayfnwr1v,yq6lorbx) = work2(tgiyxdw1(yq6lorbx),dufozmt7(yq6lor
     *bx))
23156 continue
C23157 continue
      do23158 yq6lorbx=1,kgwmz4ip 
      work2(yq6lorbx,kgwmz4ip+1) = 0.0d0
      do23160 urohxe6t=1,wy1vqfzu 
      work2(yq6lorbx,kgwmz4ip+1) = work2(yq6lorbx,kgwmz4ip+1) + hjm2ktyr
     *(urohxe6t,yq6lorbx) * work(urohxe6t,wy1vqfzu+1)
23160 continue
C23161 continue
23158 continue
C23159 continue
      do23162 yq6lorbx=1,kgwmz4ip 
      wpasjmo8g(ayfnwr1v,yq6lorbx) = work2(yq6lorbx,kgwmz4ip+1)
23162 continue
C23163 continue
      call vcholf(work2, work2(1,kgwmz4ip+1), kgwmz4ip, dvhw1ulq, oneint
     *)
      if(dvhw1ulq .ne. 1)then
      return
      endif
      if(wuwbar .ne. 0)then
      do23168 yq6lorbx=1,dimu 
      uwbar(yq6lorbx,ayfnwr1v) = work2(tgiyxdw1(yq6lorbx),dufozmt7(yq6lo
     *rbx))
23168 continue
C23169 continue
      endif
      do23170 yq6lorbx=1,kgwmz4ip 
      pasjmo8g(ayfnwr1v,yq6lorbx) = work2(yq6lorbx,kgwmz4ip+1)
23170 continue
C23171 continue
23142 continue
C23143 continue
      endif
      return
      end
      subroutine icpd0omv(enaqpzk9, he7mqnvy, gkdx5jal, grmuyvx9, ldk, k
     *uzxj1lo, nk, wy1vqfzu, jzwsy6tp, bmb, work, wmat, ifys6woa, dimw, 
     *tgiyxdw1, dufozmt7, truen)
      implicit logical (a-z)
      integer ldk, kuzxj1lo, nk, wy1vqfzu, jzwsy6tp, dimw, tgiyxdw1(*), 
     *dufozmt7(*), truen
      double precision enaqpzk9(ldk,nk*wy1vqfzu), he7mqnvy(kuzxj1lo), gk
     *dx5jal(nk+4), grmuyvx9(truen,wy1vqfzu), bmb(wy1vqfzu,wy1vqfzu), wo
     *rk(wy1vqfzu,wy1vqfzu), wmat(kuzxj1lo,dimw), ifys6woa(kuzxj1lo,wy1v
     *qfzu)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, dqlr5bse, pqzfxw4i, urohxe6t
     *, bpvaqm5z
      double precision qaltf0nz, ms0qypiw(16), g9fvdrbw(4,1)
      if(jzwsy6tp .ne. 0)then
      do23174 gp1jxzuh=1,wy1vqfzu 
      do23176 ayfnwr1v=1,kuzxj1lo 
      grmuyvx9(ayfnwr1v,gp1jxzuh) = 0.0d0
23176 continue
C23177 continue
23174 continue
C23175 continue
      endif
      qaltf0nz = 0.10d-9
      call qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      do23178 ayfnwr1v=1,kuzxj1lo 
      do23180 yq6lorbx=1,wy1vqfzu 
      do23182 gp1jxzuh=1,wy1vqfzu 
      bmb(yq6lorbx,gp1jxzuh)=0.0d0
23182 continue
C23183 continue
23180 continue
C23181 continue
      call vinterv(gkdx5jal(1), (nk+1), he7mqnvy(ayfnwr1v), dqlr5bse, pq
     *zfxw4i)
      if(pqzfxw4i.eq. 1)then
      if(he7mqnvy(ayfnwr1v) .le. (gkdx5jal(dqlr5bse)+qaltf0nz))then
      dqlr5bse=dqlr5bse-1
      else
      return
      endif
      endif
      call vbsplvd(gkdx5jal, 4, he7mqnvy(ayfnwr1v), dqlr5bse, ms0qypiw, 
     *g9fvdrbw, 1)
      yq6lorbx= dqlr5bse-4+1
      do23188 urohxe6t=yq6lorbx,yq6lorbx+3 
      call vsel(urohxe6t, urohxe6t, wy1vqfzu, nk, ldk, enaqpzk9, work)
      call o0xlszqr(wy1vqfzu, g9fvdrbw(urohxe6t-yq6lorbx+1,1) * g9fvdrbw
     *(urohxe6t-yq6lorbx+1,1), work, bmb)
23188 continue
C23189 continue
      do23190 urohxe6t=yq6lorbx,yq6lorbx+3 
      do23192 bpvaqm5z=urohxe6t+1,yq6lorbx+3 
      call vsel(urohxe6t, bpvaqm5z, wy1vqfzu, nk, ldk, enaqpzk9, work)
      call o0xlszqr(wy1vqfzu, 2.0d0 * g9fvdrbw(urohxe6t-yq6lorbx+1,1) * 
     *g9fvdrbw(bpvaqm5z-yq6lorbx+1,1), work, bmb)
23192 continue
C23193 continue
23190 continue
C23191 continue
      if(jzwsy6tp .ne. 0)then
      do23196 yq6lorbx=1,wy1vqfzu 
      grmuyvx9(ayfnwr1v,yq6lorbx) = bmb(yq6lorbx,yq6lorbx)
23196 continue
C23197 continue
      endif
      call ovjnsmt2(bmb, wmat, work, ifys6woa, wy1vqfzu, kuzxj1lo, dimw,
     * tgiyxdw1, dufozmt7, ayfnwr1v)
23178 continue
C23179 continue
      return
      end
      subroutine o0xlszqr(wy1vqfzu, g9fvdrbw, work, bmb)
      implicit logical (a-z)
      integer wy1vqfzu
      double precision g9fvdrbw, work(wy1vqfzu,wy1vqfzu), bmb(wy1vqfzu,w
     *y1vqfzu)
      integer yq6lorbx, gp1jxzuh
      do23198 yq6lorbx=1,wy1vqfzu 
      do23200 gp1jxzuh=1,wy1vqfzu 
      work(yq6lorbx,gp1jxzuh) = work(yq6lorbx,gp1jxzuh) * g9fvdrbw
23200 continue
C23201 continue
23198 continue
C23199 continue
      do23202 yq6lorbx=1,wy1vqfzu 
      do23204 gp1jxzuh=1,wy1vqfzu 
      bmb(gp1jxzuh,yq6lorbx) = bmb(gp1jxzuh,yq6lorbx) + work(gp1jxzuh,yq
     *6lorbx)
23204 continue
C23205 continue
23202 continue
C23203 continue
      return
      end
      subroutine vsel(s, t, wy1vqfzu, nk, ldk, minv, work)
      implicit logical (a-z)
      integer s, t, wy1vqfzu, nk, ldk
      double precision minv(ldk,nk*wy1vqfzu), work(wy1vqfzu,wy1vqfzu)
      integer ayfnwr1v, yq6lorbx, biuvowq2, nbj8tdsk
      do23206 ayfnwr1v=1,wy1vqfzu 
      do23208 yq6lorbx=1,wy1vqfzu 
      work(ayfnwr1v,yq6lorbx) = 0.0d0
23208 continue
C23209 continue
23206 continue
C23207 continue
      if(s .ne. t)then
      do23212 ayfnwr1v=1,wy1vqfzu 
      biuvowq2 = (s-1)*wy1vqfzu + ayfnwr1v
      do23214 yq6lorbx=1,wy1vqfzu 
      nbj8tdsk = (t-1)*wy1vqfzu + yq6lorbx
      work(ayfnwr1v,yq6lorbx) = minv(ldk-(nbj8tdsk-biuvowq2), nbj8tdsk)
23214 continue
C23215 continue
23212 continue
C23213 continue
      else
      do23216 ayfnwr1v=1,wy1vqfzu 
      biuvowq2 = (s-1)*wy1vqfzu + ayfnwr1v
      do23218 yq6lorbx=ayfnwr1v,wy1vqfzu 
      nbj8tdsk = (t-1)*wy1vqfzu + yq6lorbx
      work(ayfnwr1v,yq6lorbx) = minv(ldk-(nbj8tdsk-biuvowq2), nbj8tdsk)
23218 continue
C23219 continue
23216 continue
C23217 continue
      do23220 ayfnwr1v=1,wy1vqfzu 
      do23222 yq6lorbx=ayfnwr1v+1,wy1vqfzu 
      work(yq6lorbx,ayfnwr1v) = work(ayfnwr1v,yq6lorbx)
23222 continue
C23223 continue
23220 continue
C23221 continue
      endif
      return
      end
      subroutine ovjnsmt2(bmb, wmat, work, ifys6woa, wy1vqfzu, kuzxj1lo,
     * dimw, tgiyxdw1, dufozmt7, iii)
      implicit logical (a-z)
      integer wy1vqfzu, kuzxj1lo, dimw, tgiyxdw1(*), dufozmt7(*), iii
      double precision bmb(wy1vqfzu,wy1vqfzu), wmat(kuzxj1lo,dimw), work
     *(wy1vqfzu,wy1vqfzu), ifys6woa(kuzxj1lo,wy1vqfzu)
      double precision q6zdcwxk, obr6tcex
      integer yq6lorbx, gp1jxzuh, urohxe6t, bpvaqm5z
      do23224 bpvaqm5z=1,wy1vqfzu 
      do23226 yq6lorbx=1,wy1vqfzu 
      do23228 gp1jxzuh=1,wy1vqfzu 
      work(gp1jxzuh,yq6lorbx) = 0.0d0
23228 continue
C23229 continue
23226 continue
C23227 continue
      do23230 urohxe6t=1,dimw 
      obr6tcex = wmat(iii,urohxe6t)
      work(tgiyxdw1(urohxe6t),dufozmt7(urohxe6t)) = obr6tcex
      work(dufozmt7(urohxe6t),tgiyxdw1(urohxe6t)) = obr6tcex
23230 continue
C23231 continue
      q6zdcwxk = 0.0d0
      do23232 yq6lorbx=1,wy1vqfzu 
      q6zdcwxk = q6zdcwxk + bmb(bpvaqm5z,yq6lorbx) * work(yq6lorbx,bpvaq
     *m5z)
23232 continue
C23233 continue
      ifys6woa(iii,bpvaqm5z) = q6zdcwxk
23224 continue
C23225 continue
      return
      end
      subroutine vicb2(enaqpzk9, wpuarq2m, d, uu, wy1vqfzu, kuzxj1lo)
      implicit logical (a-z)
      integer wy1vqfzu, kuzxj1lo
      double precision enaqpzk9(wy1vqfzu+1,kuzxj1lo), wpuarq2m(wy1vqfzu+
     *1,kuzxj1lo), d(kuzxj1lo), uu(wy1vqfzu+1,wy1vqfzu+1)
      integer ayfnwr1v, gp1jxzuh, lsvdbx3tk, uplim, sedf7mxb, hofjnx2e, 
     *kij0gwer
      enaqpzk9(wy1vqfzu+1,kuzxj1lo) = 1.0d0 / d(kuzxj1lo)
      hofjnx2e = wy1vqfzu+1
      sedf7mxb = kuzxj1lo+1 - hofjnx2e
      do23234 kij0gwer=sedf7mxb,kuzxj1lo 
      do23236 ayfnwr1v=1,hofjnx2e 
      uu(ayfnwr1v, kij0gwer-sedf7mxb+1) = wpuarq2m(ayfnwr1v, kij0gwer)
23236 continue
C23237 continue
23234 continue
C23235 continue
      ayfnwr1v = kuzxj1lo-1 
23238 if(.not.(ayfnwr1v .ge. 1))goto 23240
      if(wy1vqfzu .lt. kuzxj1lo-ayfnwr1v)then
      uplim = wy1vqfzu
      else
      uplim = kuzxj1lo-ayfnwr1v
      endif
      lsvdbx3tk=1
23243 if(.not.(lsvdbx3tk .le. uplim))goto 23245
      enaqpzk9(-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) = 0.0d0
      gp1jxzuh=1
23246 if(.not.(gp1jxzuh .le. lsvdbx3tk))goto 23248
      enaqpzk9(-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) = enaqpzk9(-lsv
     *dbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) - uu(-gp1jxzuh+wy1vqfzu+1,ay
     *fnwr1v+gp1jxzuh -sedf7mxb+1) * enaqpzk9(gp1jxzuh-lsvdbx3tk+wy1vqfz
     *u+1,ayfnwr1v+lsvdbx3tk)
      gp1jxzuh=gp1jxzuh+1
      goto 23246
23248 continue
23249 if(.not.(gp1jxzuh .le. uplim))goto 23251
      enaqpzk9(-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) = enaqpzk9(-lsv
     *dbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) - uu(-gp1jxzuh+wy1vqfzu+1,ay
     *fnwr1v+gp1jxzuh -sedf7mxb+1) * enaqpzk9(lsvdbx3tk-gp1jxzuh+wy1vqfz
     *u+1,ayfnwr1v+gp1jxzuh)
      gp1jxzuh=gp1jxzuh+1
      goto 23249
23251 continue
      lsvdbx3tk=lsvdbx3tk+1
      goto 23243
23245 continue
      enaqpzk9(wy1vqfzu+1,ayfnwr1v) = 1.0d0 / d(ayfnwr1v)
      lsvdbx3tk = 1
23252 if(.not.(lsvdbx3tk .le. uplim))goto 23254
      enaqpzk9(wy1vqfzu+1,ayfnwr1v) = enaqpzk9(wy1vqfzu+1,ayfnwr1v) - uu
     *(-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk -sedf7mxb+1) * enaqpzk9(
     *-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk)
      lsvdbx3tk=lsvdbx3tk+1
      goto 23252
23254 continue
      if(ayfnwr1v .eq. sedf7mxb)then
      sedf7mxb = sedf7mxb-1
      if(sedf7mxb .lt. 1)then
      sedf7mxb = 1
      else
      kij0gwer=hofjnx2e-1
23259 if(.not.(kij0gwer .ge. 1))goto 23261
      gp1jxzuh=1
23262 if(.not.(gp1jxzuh .le. hofjnx2e))goto 23264
      uu(gp1jxzuh,kij0gwer+1) = uu(gp1jxzuh,kij0gwer)
      gp1jxzuh=gp1jxzuh+1
      goto 23262
23264 continue
      kij0gwer=kij0gwer-1
      goto 23259
23261 continue
      gp1jxzuh=1
23265 if(.not.(gp1jxzuh .le. hofjnx2e))goto 23267
      uu(gp1jxzuh,1) = wpuarq2m(gp1jxzuh,sedf7mxb)
      gp1jxzuh=gp1jxzuh+1
      goto 23265
23267 continue
      endif
      endif
      ayfnwr1v = ayfnwr1v-1
      goto 23238
23240 continue
      return
      end
      subroutine ewg7qruh(sjwyig9tto,tlgduey8,wmat, kuzxj1lo,wy1vqfzu,ez
     *lgm2up,nef, wbkq9zyi,dof,smo,cov, s0, xin,yin,rbne6ouj,win, work1,
     *work3, dimw, fbd5yktj, ldk, info, yzoe1rsp, sgdub, rpyis2kc, zv2xf
     *hei, acpios9q,tgiyxdw1,dufozmt7, bmb, ifys6woa, wkmm, iz2nbfjc,kgw
     *mz4ip,ges1xpkr, hjm2ktyr, beta, fasrkub3, sout, r0oydcxb, ub4xioar
     *, effect, uwin)
      implicit logical (a-z)
      integer kuzxj1lo,wy1vqfzu,ezlgm2up(kuzxj1lo),nef, dimw, fbd5yktj, 
     *ldk, info, yzoe1rsp, acpios9q,tgiyxdw1(*),dufozmt7(*), iz2nbfjc, k
     *gwmz4ip, ges1xpkr(kgwmz4ip*2)
      double precision sjwyig9tto(kuzxj1lo), tlgduey8(kuzxj1lo,wy1vqfzu)
     *, wmat(kuzxj1lo,dimw), wbkq9zyi(kgwmz4ip), dof(kgwmz4ip), smo(kuzx
     *j1lo,kgwmz4ip), cov(kuzxj1lo,kgwmz4ip)
      double precision s0(2*kgwmz4ip, 2*kgwmz4ip,2)
      double precision work1(*), work3(*), sgdub(*), rpyis2kc(*), zv2xfh
     *ei(acpios9q+4)
      double precision xin(nef), yin(nef,wy1vqfzu), rbne6ouj(nef,wy1vqfz
     *u), win(nef,*), bmb(*), ifys6woa(nef,kgwmz4ip), wkmm(wy1vqfzu,wy1v
     *qfzu,16), hjm2ktyr(wy1vqfzu,kgwmz4ip)
      double precision beta(2*kgwmz4ip), fasrkub3(2*kgwmz4ip), sout(nef,
     *kgwmz4ip), r0oydcxb(kgwmz4ip,nef), ub4xioar(kgwmz4ip,nef), effect(
     *nef*kgwmz4ip), uwin(*)
      integer dimwin
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, rutyk8mg, xjc4ywlh, job, qem
     *j9asg, dvhw1ulq
      integer oneint
      double precision xmin, xrange, pvofyg8z
      oneint = 1
      if(iz2nbfjc .eq. 1)then
      dimwin = dimw
      else
      dimwin = kgwmz4ip*(kgwmz4ip+1)/2
      endif
      call qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      call vsuff9(kuzxj1lo,nef,ezlgm2up, sjwyig9tto,tlgduey8,wmat, xin,y
     *in,win,uwin,rbne6ouj, wy1vqfzu, dimw, dimwin, tgiyxdw1, dufozmt7, 
     *wkmm, wkmm(1,1,3), hjm2ktyr, kgwmz4ip, iz2nbfjc, oneint, dvhw1ulq)
      if(dvhw1ulq .ne. 1)then
      return
      endif
      xmin = xin(1)
      xrange = xin(nef)-xin(1)
      do23272 ayfnwr1v=1,nef 
      xin(ayfnwr1v) = (xin(ayfnwr1v)-xmin)/xrange
23272 continue
C23273 continue
      ldk = 4*kgwmz4ip
      fbd5yktj = 0
      do23274 yq6lorbx=1,kgwmz4ip 
      if(wbkq9zyi(yq6lorbx) .eq. 0.0d0)then
      dof(yq6lorbx) = dof(yq6lorbx) + 1.0d0
      endif
23274 continue
C23275 continue
      call qpsedg8xf(tgiyxdw1, dufozmt7, kgwmz4ip)
      call vsplin(xin,rbne6ouj,win,nef,zv2xfhei, acpios9q,ldk,kgwmz4ip,d
     *imwin, tgiyxdw1,dufozmt7, wkmm, wbkq9zyi, info, sout, rpyis2kc, wo
     *rk3(1), work3(1+acpios9q*kgwmz4ip*ldk), sgdub, cov, yzoe1rsp, bmb,
     * ifys6woa, dof, work1, fbd5yktj, kuzxj1lo)
      do23278 yq6lorbx=1,kgwmz4ip 
      dof(yq6lorbx) = -1.0d0
      do23280 ayfnwr1v=1,nef 
      dof(yq6lorbx)=dof(yq6lorbx)+ifys6woa(ayfnwr1v,yq6lorbx)
23280 continue
C23281 continue
23278 continue
C23279 continue
      if(kgwmz4ip .ge. 1)then
      pvofyg8z = 1.0d-7
      rutyk8mg = nef*kgwmz4ip
      xjc4ywlh = 2*kgwmz4ip
      job = 101
      info = 1
      call x6kanjdh(xin, work3, nef, kgwmz4ip)
      call qpsedg8xf(tgiyxdw1, dufozmt7, kgwmz4ip)
      call mux17f(uwin, work3, kgwmz4ip, xjc4ywlh, nef, wkmm(1,1,1), wkm
     *m(1,1,2), tgiyxdw1, dufozmt7, dimwin, rutyk8mg)
      do23284 gp1jxzuh=1,xjc4ywlh 
      ges1xpkr(gp1jxzuh) = gp1jxzuh
23284 continue
C23285 continue
      call vqrdca(work3,rutyk8mg,rutyk8mg,xjc4ywlh,fasrkub3,ges1xpkr,wor
     *k1,qemj9asg,pvofyg8z)
      call qpsedg8xf(tgiyxdw1, dufozmt7, kgwmz4ip)
      call mux22f(uwin,sout,r0oydcxb,dimwin,tgiyxdw1,dufozmt7,nef,kgwmz4
     *ip,wkmm)
      call vdqrsl(work3,rutyk8mg,rutyk8mg,qemj9asg,fasrkub3,r0oydcxb,wor
     *k1(1),effect,beta, work1(1),ub4xioar,job,info)
      call vbksf(uwin,ub4xioar,kgwmz4ip,nef,wkmm,tgiyxdw1,dufozmt7,dimwi
     *n)
      if(yzoe1rsp .ne. 0)then
      call vrinvf9(work3, rutyk8mg, xjc4ywlh, dvhw1ulq, s0(1,1,1), s0(1,
     *1,2))
      if(dvhw1ulq .ne. 1)then
      return
      endif
      do23290 yq6lorbx=1,kgwmz4ip 
      do23292 ayfnwr1v=1,nef 
      cov(ayfnwr1v,yq6lorbx) = cov(ayfnwr1v,yq6lorbx) - s0(yq6lorbx,yq6l
     *orbx,1) - xin(ayfnwr1v) * (2.0d0 * s0(yq6lorbx,yq6lorbx+kgwmz4ip,1
     *) + xin(ayfnwr1v) * s0(yq6lorbx+kgwmz4ip,yq6lorbx+kgwmz4ip,1))
23292 continue
C23293 continue
23290 continue
C23291 continue
      endif
      else
      call dsrt0gem(nef, xin, win, sout, ub4xioar, cov, yzoe1rsp)
      endif
      do23294 ayfnwr1v=1,nef 
      do23296 yq6lorbx=1,kgwmz4ip 
      sout(ayfnwr1v,yq6lorbx) = sout(ayfnwr1v,yq6lorbx) - ub4xioar(yq6lo
     *rbx,ayfnwr1v)
23296 continue
C23297 continue
23294 continue
C23295 continue
      do23298 yq6lorbx=1,kgwmz4ip 
      call shm8ynte(kuzxj1lo, nef, ezlgm2up, sout(1,yq6lorbx), smo(1,yq6
     *lorbx))
23298 continue
C23299 continue
      return
      end
      subroutine x6kanjdh(he7mqnvy, xout, kuzxj1lo, wy1vqfzu)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu
      double precision he7mqnvy(kuzxj1lo), xout(*)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, iptr
      iptr=1
      do23300 yq6lorbx=1,wy1vqfzu 
      do23302 ayfnwr1v=1,kuzxj1lo 
      do23304 gp1jxzuh=1,wy1vqfzu 
      if(yq6lorbx .eq. gp1jxzuh)then
      xout(iptr) = 1.0d0
      else
      xout(iptr) = 0.0d0
      endif
      iptr=iptr+1
23304 continue
C23305 continue
23302 continue
C23303 continue
23300 continue
C23301 continue
      do23308 yq6lorbx=1,wy1vqfzu 
      do23310 ayfnwr1v=1,kuzxj1lo 
      do23312 gp1jxzuh=1,wy1vqfzu 
      if(yq6lorbx .eq. gp1jxzuh)then
      xout(iptr) = he7mqnvy(ayfnwr1v)
      else
      xout(iptr) = 0.0d0
      endif
      iptr=iptr+1
23312 continue
C23313 continue
23310 continue
C23311 continue
23308 continue
C23309 continue
      return
      end
      double precision function rd9beyfk(kuzxj1lo, bhcji9gl, m0ibglfx, p
     *o8rwsmy)
      integer kuzxj1lo
      double precision bhcji9gl(kuzxj1lo), m0ibglfx(kuzxj1lo), po8rwsmy(
     *kuzxj1lo)
      integer ayfnwr1v
      double precision lm9vcjob, rxeqjn0y, work
      rxeqjn0y = 0.0d0
      lm9vcjob = 0.0d0
      do23316 ayfnwr1v=1,kuzxj1lo 
      work = bhcji9gl(ayfnwr1v) - m0ibglfx(ayfnwr1v)
      rxeqjn0y = rxeqjn0y + po8rwsmy(ayfnwr1v)*work*work
      lm9vcjob = lm9vcjob + po8rwsmy(ayfnwr1v)
23316 continue
C23317 continue
      if(lm9vcjob .gt. 0.0d0)then
      rd9beyfk=rxeqjn0y/lm9vcjob
      else
      rd9beyfk=0.0d0
      endif
      return
      end
      subroutine pitmeh0q(kuzxj1lo, bhcji9gl, po8rwsmy, lfu2qhid, lm9vcj
     *ob)
      implicit logical (a-z)
      integer kuzxj1lo
      double precision bhcji9gl(kuzxj1lo), po8rwsmy(kuzxj1lo), lfu2qhid,
     * lm9vcjob
      double precision rxeqjn0y
      integer ayfnwr1v
      lm9vcjob = 0.0d0
      rxeqjn0y = 0.0d0
      do23320 ayfnwr1v=1,kuzxj1lo 
      rxeqjn0y = rxeqjn0y + bhcji9gl(ayfnwr1v) * po8rwsmy(ayfnwr1v)
      lm9vcjob = lm9vcjob + po8rwsmy(ayfnwr1v)
23320 continue
C23321 continue
      if(lm9vcjob .gt. 0.0d0)then
      lfu2qhid = rxeqjn0y / lm9vcjob
      else
      lfu2qhid = 0.0d0
      endif
      return
      end
      subroutine dsrt0gem(kuzxj1lo, x, w, bhcji9gl, ub4xioar, cov, yzoe1
     *rsp)
      implicit logical (a-z)
      integer kuzxj1lo
      integer yzoe1rsp
      double precision x(kuzxj1lo), w(kuzxj1lo), bhcji9gl(kuzxj1lo), ub4
     *xioar(kuzxj1lo)
      double precision cov(kuzxj1lo,*)
      integer ayfnwr1v
      double precision pasjmo8g, pygsw6ko, q6zdcwxk, nsum, eck8vubt, int
     *erc, bzmd6ftv, hofjnx2e, lm9vcjob
      call pitmeh0q(kuzxj1lo,bhcji9gl,w,pasjmo8g, lm9vcjob)
      call pitmeh0q(kuzxj1lo,x,w,pygsw6ko, lm9vcjob)
      nsum = 0.0d0
      q6zdcwxk = 0.0d0
      do23324 ayfnwr1v=1,kuzxj1lo 
      hofjnx2e = x(ayfnwr1v)-pygsw6ko
      nsum = nsum + hofjnx2e * (bhcji9gl(ayfnwr1v)-pasjmo8g) * w(ayfnwr1
     *v)
      hofjnx2e = hofjnx2e * hofjnx2e
      q6zdcwxk = q6zdcwxk + hofjnx2e * w(ayfnwr1v)
23324 continue
C23325 continue
      eck8vubt = nsum/q6zdcwxk
      interc = pasjmo8g - eck8vubt * pygsw6ko
      do23326 ayfnwr1v=1,kuzxj1lo 
      ub4xioar(ayfnwr1v) = interc + eck8vubt * x(ayfnwr1v)
23326 continue
C23327 continue
      bzmd6ftv = interc + eck8vubt * x(1)
      if(yzoe1rsp .ne. 0)then
      do23330 ayfnwr1v=1,kuzxj1lo 
      hofjnx2e = x(ayfnwr1v)-pygsw6ko
      if(w(ayfnwr1v) .gt. 0.0d0)then
      cov(ayfnwr1v,1) = cov(ayfnwr1v,1) - 1.0d0/lm9vcjob - hofjnx2e * ho
     *fjnx2e / q6zdcwxk
      else
      cov(ayfnwr1v,1) = 0.0d0
      endif
23330 continue
C23331 continue
      endif
      return
      end
      subroutine shm8ynte(kuzxj1lo, p, ezlgm2up, pygsw6ko, x)
      implicit logical (a-z)
      integer kuzxj1lo, p, ezlgm2up(kuzxj1lo)
      double precision pygsw6ko(p), x(kuzxj1lo)
      integer ayfnwr1v
      do23334 ayfnwr1v=1,kuzxj1lo 
      x(ayfnwr1v) = pygsw6ko(ezlgm2up(ayfnwr1v))
23334 continue
C23335 continue
      return
      end
      subroutine vankcghz2l2(x, kuzxj1lo, ankcghz2, rvy1fpli, ukgwt7na)
      implicit logical (a-z)
      integer kuzxj1lo, rvy1fpli, ukgwt7na
      double precision x(kuzxj1lo), ankcghz2(kuzxj1lo)
      integer ndk, yq6lorbx
      if(ukgwt7na .eq. 0)then
      if(kuzxj1lo .le. 40)then
      ndk = kuzxj1lo
      else
      ndk = 40 + nint(dexp(0.25d0 * dlog(kuzxj1lo-40.0d0)))
      endif
      else
      ndk = rvy1fpli - 6
      endif
      rvy1fpli = ndk + 6
      do23340 yq6lorbx = 1,3 
      ankcghz2(yq6lorbx) = x(1) 
23340 continue
C23341 continue
      do23342 yq6lorbx = 1,ndk 
      ankcghz2(yq6lorbx+3) = x( 1 + (yq6lorbx-1)*(kuzxj1lo-1)/(ndk-1) ) 
23342 continue
C23343 continue
      do23344 yq6lorbx = 1,3 
      ankcghz2(ndk+3+yq6lorbx) = x(kuzxj1lo) 
23344 continue
C23345 continue
      return
      end
      subroutine pankcghz2l2(ankcghz2, kuzxj1lo, zo8wpibx, tol)
      implicit logical (a-z)
      integer kuzxj1lo, zo8wpibx(kuzxj1lo)
      double precision ankcghz2(kuzxj1lo), tol
      integer ayfnwr1v, cjop5bwm
      do23346 ayfnwr1v=1,4 
      zo8wpibx(ayfnwr1v) = 1
23346 continue
C23347 continue
      cjop5bwm = 4
      do23348 ayfnwr1v=5,(kuzxj1lo-4) 
      if((ankcghz2(ayfnwr1v) - ankcghz2(cjop5bwm) .ge. tol) .and. (ankcg
     *hz2(kuzxj1lo) - ankcghz2(ayfnwr1v) .ge. tol))then
      zo8wpibx(ayfnwr1v) = 1
      cjop5bwm = ayfnwr1v
      else
      zo8wpibx(ayfnwr1v) = 0
      endif
23348 continue
C23349 continue
      do23352 ayfnwr1v=(kuzxj1lo-3),kuzxj1lo 
      zo8wpibx(ayfnwr1v) = 1
23352 continue
C23353 continue
      return
      end
