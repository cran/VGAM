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
23003 continue
23000 continue
23001 continue
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
23007 continue
23004 continue
23005 continue
      do23008 ayfnwr1v=1,(nk-1) 
      do23010 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk-wy1vqfzu,(ayfnwr1v-0)*wy1vqfzu+yq6lorbx) = osiz4fxy(l
     *dk-wy1vqfzu,(ayfnwr1v-0)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorbx) *
     * sgmat(ayfnwr1v,2)
23010 continue
23011 continue
23008 continue
23009 continue
      do23012 ayfnwr1v=1,(nk-2) 
      do23014 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk-2*wy1vqfzu,(ayfnwr1v+1)*wy1vqfzu+yq6lorbx) = osiz4fxy
     *(ldk-2*wy1vqfzu,(ayfnwr1v+1)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorb
     *x) * sgmat(ayfnwr1v,3)
23014 continue
23015 continue
23012 continue
23013 continue
      do23016 ayfnwr1v=1,(nk-3) 
      do23018 yq6lorbx=1,wy1vqfzu 
      osiz4fxy(ldk-3*wy1vqfzu,(ayfnwr1v+2)*wy1vqfzu+yq6lorbx) = osiz4fxy
     *(ldk-3*wy1vqfzu,(ayfnwr1v+2)*wy1vqfzu+yq6lorbx) + wbkq9zyi(yq6lorb
     *x) * sgmat(ayfnwr1v,4)
23018 continue
23019 continue
23016 continue
23017 continue
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
23021 continue
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
23039 continue
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
23043 continue
      if(yzoe1rsp .ne. 0)then
      do23046 ayfnwr1v=1,kuzxj1lo 
      ui8ysltq(ayfnwr1v,yq6lorbx) = ifys6woa(ayfnwr1v,yq6lorbx) / wmat(a
     *yfnwr1v,yq6lorbx)
23046 continue
23047 continue
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
23051 continue
      endif
      if(fbd5yktj .ne. 0)then
      return
      endif
      endif
23024 continue
23025 continue
      if((wy1vqfzu .eq. 1) .or. (dimw .eq. wy1vqfzu))then
      return
      endif
      do23056 ayfnwr1v=1,nk 
      do23058 yq6lorbx=1,wy1vqfzu 
      btwy(yq6lorbx,ayfnwr1v)=0.0d0
23058 continue
23059 continue
23056 continue
23057 continue
      do23060 ayfnwr1v=1,(nk*wy1vqfzu) 
      do23062 yq6lorbx=1,ldk 
      osiz4fxy(yq6lorbx,ayfnwr1v) = 0.0d0
23062 continue
23063 continue
23060 continue
23061 continue
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
23071 continue
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
23073 continue
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
23075 continue
      call ybnagt8k(ayfnwr1v, yq6lorbx, 0, g9fvdrbw, osiz4fxy, wmat, 3, 
     *3, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      call ybnagt8k(ayfnwr1v, yq6lorbx, 1, g9fvdrbw, osiz4fxy, wmat, 3, 
     *4, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
      yq6lorbx= dqlr5bse-4+4
      do23076 urohxe6t=1,wy1vqfzu 
      btwy(urohxe6t,yq6lorbx)=btwy(urohxe6t,yq6lorbx) + rbne6ouj(ayfnwr1
     *v,urohxe6t) * g9fvdrbw(4,1)
23076 continue
23077 continue
      call ybnagt8k(ayfnwr1v, yq6lorbx, 0, g9fvdrbw, osiz4fxy, wmat, 4, 
     *4, wy1vqfzu, ldk, dimw, kuzxj1lo, nk, tgiyxdw1, dufozmt7)
23064 continue
23065 continue
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
23083 continue
23080 continue
23081 continue
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
23087 continue
23084 continue
23085 continue
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
23095 continue
      do23096 yq6lorbx=1,wy1vqfzu 
      do23098 ayfnwr1v=1,nef 
      wpasjmo8g(ayfnwr1v,yq6lorbx) = 0.0d0
23098 continue
23099 continue
23096 continue
23097 continue
      do23100 yq6lorbx=1,dimw 
      do23102 ayfnwr1v=1,nef 
      wbar(ayfnwr1v,yq6lorbx) = 0.0d0
23102 continue
23103 continue
23100 continue
23101 continue
      if(dimw .ne. imk5wjxg)then
      do23106 gp1jxzuh=1,wy1vqfzu 
      do23108 yq6lorbx=1,wy1vqfzu 
      work(yq6lorbx,gp1jxzuh) = 0.0d0
23108 continue
23109 continue
23106 continue
23107 continue
      endif
      do23110 ayfnwr1v=1,kuzxj1lo 
      do23112 yq6lorbx=1,dimw 
      work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx)) = wmat(ayfnwr1v,yq6lor
     *bx)
      work(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx)) = work(tgiyxdw1(yq6lor
     *bx),dufozmt7(yq6lorbx))
23112 continue
23113 continue
      do23114 yq6lorbx=1,wy1vqfzu 
      do23116 gp1jxzuh=1,wy1vqfzu 
      wpasjmo8g(ezlgm2up(ayfnwr1v),yq6lorbx) = wpasjmo8g(ezlgm2up(ayfnwr
     *1v),yq6lorbx) + work(yq6lorbx,gp1jxzuh)*tlgduey8(ayfnwr1v,gp1jxzuh
     *)
23116 continue
23117 continue
23114 continue
23115 continue
      do23118 yq6lorbx=1,dimw 
      wbar(ezlgm2up(ayfnwr1v),yq6lorbx) = wbar(ezlgm2up(ayfnwr1v),yq6lor
     *bx) + wmat(ayfnwr1v,yq6lorbx)
23118 continue
23119 continue
23110 continue
23111 continue
      dvhw1ulq = 1
      if(iz2nbfjc .eq. 1)then
      do23122 ayfnwr1v=1,nef 
      do23124 yq6lorbx=1,dimw 
      work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx)) = wbar(ayfnwr1v,yq6lor
     *bx)
      work(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx)) = work(tgiyxdw1(yq6lor
     *bx),dufozmt7(yq6lorbx))
23124 continue
23125 continue
      do23126 yq6lorbx=1,wy1vqfzu 
      work(yq6lorbx,wy1vqfzu+1)=wpasjmo8g(ayfnwr1v,yq6lorbx)
23126 continue
23127 continue
      call vcholf(work, work(1,wy1vqfzu+1), wy1vqfzu, dvhw1ulq, oneint)
      if(dvhw1ulq .ne. 1)then
      return
      endif
      if(wuwbar .ne. 0)then
      do23132 yq6lorbx=1,dimw 
      uwbar(yq6lorbx,ayfnwr1v) = work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lor
     *bx))
23132 continue
23133 continue
      endif
      do23134 yq6lorbx=1,wy1vqfzu 
      pasjmo8g(ayfnwr1v,yq6lorbx)=work(yq6lorbx,wy1vqfzu+1)
23134 continue
23135 continue
23122 continue
23123 continue
      else
      if(dimw .ne. imk5wjxg)then
      do23138 yq6lorbx=1,wy1vqfzu 
      do23140 gp1jxzuh=1,wy1vqfzu 
      work(yq6lorbx,gp1jxzuh) = 0.0d0
23140 continue
23141 continue
23138 continue
23139 continue
      endif
      do23142 ayfnwr1v=1,nef 
      call qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      do23144 yq6lorbx=1,dimw 
      work(tgiyxdw1(yq6lorbx),dufozmt7(yq6lorbx)) = wbar(ayfnwr1v,yq6lor
     *bx)
      work(dufozmt7(yq6lorbx),tgiyxdw1(yq6lorbx)) = work(tgiyxdw1(yq6lor
     *bx),dufozmt7(yq6lorbx))
23144 continue
23145 continue
      do23146 yq6lorbx=1,wy1vqfzu 
      work(yq6lorbx,wy1vqfzu+1)=wpasjmo8g(ayfnwr1v,yq6lorbx)
23146 continue
23147 continue
      do23148 yq6lorbx=1,kgwmz4ip 
      do23150 gp1jxzuh=yq6lorbx,kgwmz4ip 
      work2(yq6lorbx,gp1jxzuh) = 0.0d0
      do23152 urohxe6t=1,wy1vqfzu 
      do23154 bpvaqm5z=1,wy1vqfzu 
      work2(yq6lorbx,gp1jxzuh) = work2(yq6lorbx,gp1jxzuh) + hjm2ktyr(uro
     *hxe6t,yq6lorbx) * work(urohxe6t,bpvaqm5z) * hjm2ktyr(bpvaqm5z,gp1j
     *xzuh)
23154 continue
23155 continue
23152 continue
23153 continue
23150 continue
23151 continue
23148 continue
23149 continue
      call qpsedg8xf(tgiyxdw1, dufozmt7, kgwmz4ip)
      do23156 yq6lorbx=1,dimu 
      wbar(ayfnwr1v,yq6lorbx) = work2(tgiyxdw1(yq6lorbx),dufozmt7(yq6lor
     *bx))
23156 continue
23157 continue
      do23158 yq6lorbx=1,kgwmz4ip 
      work2(yq6lorbx,kgwmz4ip+1) = 0.0d0
      do23160 urohxe6t=1,wy1vqfzu 
      work2(yq6lorbx,kgwmz4ip+1) = work2(yq6lorbx,kgwmz4ip+1) + hjm2ktyr
     *(urohxe6t,yq6lorbx) * work(urohxe6t,wy1vqfzu+1)
23160 continue
23161 continue
23158 continue
23159 continue
      do23162 yq6lorbx=1,kgwmz4ip 
      wpasjmo8g(ayfnwr1v,yq6lorbx) = work2(yq6lorbx,kgwmz4ip+1)
23162 continue
23163 continue
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
23169 continue
      endif
      do23170 yq6lorbx=1,kgwmz4ip 
      pasjmo8g(ayfnwr1v,yq6lorbx) = work2(yq6lorbx,kgwmz4ip+1)
23170 continue
23171 continue
23142 continue
23143 continue
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
23177 continue
23174 continue
23175 continue
      endif
      qaltf0nz = 0.10d-9
      call qpsedg8xf(tgiyxdw1, dufozmt7, wy1vqfzu)
      do23178 ayfnwr1v=1,kuzxj1lo 
      do23180 yq6lorbx=1,wy1vqfzu 
      do23182 gp1jxzuh=1,wy1vqfzu 
      bmb(yq6lorbx,gp1jxzuh)=0.0d0
23182 continue
23183 continue
23180 continue
23181 continue
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
23189 continue
      do23190 urohxe6t=yq6lorbx,yq6lorbx+3 
      do23192 bpvaqm5z=urohxe6t+1,yq6lorbx+3 
      call vsel(urohxe6t, bpvaqm5z, wy1vqfzu, nk, ldk, enaqpzk9, work)
      call o0xlszqr(wy1vqfzu, 2.0d0 * g9fvdrbw(urohxe6t-yq6lorbx+1,1) * 
     *g9fvdrbw(bpvaqm5z-yq6lorbx+1,1), work, bmb)
23192 continue
23193 continue
23190 continue
23191 continue
      if(jzwsy6tp .ne. 0)then
      do23196 yq6lorbx=1,wy1vqfzu 
      grmuyvx9(ayfnwr1v,yq6lorbx) = bmb(yq6lorbx,yq6lorbx)
23196 continue
23197 continue
      endif
      call ovjnsmt2(bmb, wmat, work, ifys6woa, wy1vqfzu, kuzxj1lo, dimw,
     * tgiyxdw1, dufozmt7, ayfnwr1v)
23178 continue
23179 continue
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
23201 continue
23198 continue
23199 continue
      do23202 yq6lorbx=1,wy1vqfzu 
      do23204 gp1jxzuh=1,wy1vqfzu 
      bmb(gp1jxzuh,yq6lorbx) = bmb(gp1jxzuh,yq6lorbx) + work(gp1jxzuh,yq
     *6lorbx)
23204 continue
23205 continue
23202 continue
23203 continue
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
23209 continue
23206 continue
23207 continue
      if(s .ne. t)then
      do23212 ayfnwr1v=1,wy1vqfzu 
      biuvowq2 = (s-1)*wy1vqfzu + ayfnwr1v
      do23214 yq6lorbx=1,wy1vqfzu 
      nbj8tdsk = (t-1)*wy1vqfzu + yq6lorbx
      work(ayfnwr1v,yq6lorbx) = minv(ldk-(nbj8tdsk-biuvowq2), nbj8tdsk)
23214 continue
23215 continue
23212 continue
23213 continue
      else
      do23216 ayfnwr1v=1,wy1vqfzu 
      biuvowq2 = (s-1)*wy1vqfzu + ayfnwr1v
      do23218 yq6lorbx=ayfnwr1v,wy1vqfzu 
      nbj8tdsk = (t-1)*wy1vqfzu + yq6lorbx
      work(ayfnwr1v,yq6lorbx) = minv(ldk-(nbj8tdsk-biuvowq2), nbj8tdsk)
23218 continue
23219 continue
23216 continue
23217 continue
      do23220 ayfnwr1v=1,wy1vqfzu 
      do23222 yq6lorbx=ayfnwr1v+1,wy1vqfzu 
      work(yq6lorbx,ayfnwr1v) = work(ayfnwr1v,yq6lorbx)
23222 continue
23223 continue
23220 continue
23221 continue
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
23229 continue
23226 continue
23227 continue
      do23230 urohxe6t=1,dimw 
      obr6tcex = wmat(iii,urohxe6t)
      work(tgiyxdw1(urohxe6t),dufozmt7(urohxe6t)) = obr6tcex
      work(dufozmt7(urohxe6t),tgiyxdw1(urohxe6t)) = obr6tcex
23230 continue
23231 continue
      q6zdcwxk = 0.0d0
      do23232 yq6lorbx=1,wy1vqfzu 
      q6zdcwxk = q6zdcwxk + bmb(bpvaqm5z,yq6lorbx) * work(yq6lorbx,bpvaq
     *m5z)
23232 continue
23233 continue
      ifys6woa(iii,bpvaqm5z) = q6zdcwxk
23224 continue
23225 continue
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
23237 continue
23234 continue
23235 continue
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
23247 gp1jxzuh=gp1jxzuh+1
      goto 23246
23248 continue
23249 if(.not.(gp1jxzuh .le. uplim))goto 23251
      enaqpzk9(-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) = enaqpzk9(-lsv
     *dbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk) - uu(-gp1jxzuh+wy1vqfzu+1,ay
     *fnwr1v+gp1jxzuh -sedf7mxb+1) * enaqpzk9(lsvdbx3tk-gp1jxzuh+wy1vqfz
     *u+1,ayfnwr1v+gp1jxzuh)
23250 gp1jxzuh=gp1jxzuh+1
      goto 23249
23251 continue
23244 lsvdbx3tk=lsvdbx3tk+1
      goto 23243
23245 continue
      enaqpzk9(wy1vqfzu+1,ayfnwr1v) = 1.0d0 / d(ayfnwr1v)
      lsvdbx3tk = 1
23252 if(.not.(lsvdbx3tk .le. uplim))goto 23254
      enaqpzk9(wy1vqfzu+1,ayfnwr1v) = enaqpzk9(wy1vqfzu+1,ayfnwr1v) - uu
     *(-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk -sedf7mxb+1) * enaqpzk9(
     *-lsvdbx3tk+wy1vqfzu+1,ayfnwr1v+lsvdbx3tk)
23253 lsvdbx3tk=lsvdbx3tk+1
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
23263 gp1jxzuh=gp1jxzuh+1
      goto 23262
23264 continue
23260 kij0gwer=kij0gwer-1
      goto 23259
23261 continue
      gp1jxzuh=1
23265 if(.not.(gp1jxzuh .le. hofjnx2e))goto 23267
      uu(gp1jxzuh,1) = wpuarq2m(gp1jxzuh,sedf7mxb)
23266 gp1jxzuh=gp1jxzuh+1
      goto 23265
23267 continue
      endif
      endif
23239 ayfnwr1v = ayfnwr1v-1
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
23273 continue
      ldk = 4*kgwmz4ip
      fbd5yktj = 0
      do23274 yq6lorbx=1,kgwmz4ip 
      if(wbkq9zyi(yq6lorbx) .eq. 0.0d0)then
      dof(yq6lorbx) = dof(yq6lorbx) + 1.0d0
      endif
23274 continue
23275 continue
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
23281 continue
23278 continue
23279 continue
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
23285 continue
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
23293 continue
23290 continue
23291 continue
      endif
      else
      call dsrt0gem(nef, xin, win, sout, ub4xioar, cov, yzoe1rsp)
      endif
      do23294 ayfnwr1v=1,nef 
      do23296 yq6lorbx=1,kgwmz4ip 
      sout(ayfnwr1v,yq6lorbx) = sout(ayfnwr1v,yq6lorbx) - ub4xioar(yq6lo
     *rbx,ayfnwr1v)
23296 continue
23297 continue
23294 continue
23295 continue
      do23298 yq6lorbx=1,kgwmz4ip 
      call shm8ynte(kuzxj1lo, nef, ezlgm2up, sout(1,yq6lorbx), smo(1,yq6
     *lorbx))
23298 continue
23299 continue
      return
      end
      subroutine vbfa( n,wy1vqfzu,psdvgce3, he7mqnvy,tlgduey8,wmat,wbkq9
     *zyi,dof, ezlgm2up,nef,which, ub4xioar,kispwgx3,m0ibglfx,s0, beta,c
     *ov,zpcqv3uj, vc6hatuj,fasrkub3, ges1xpkr, xbig, wpuarq2m, hjm2ktyr
     *, jnxpuym2, hnpt1zym, fzm1ihwj, iz2nbfjc, work1, wk2, wkmm, work3,
     * sgdub, bmb, ifys6woa, mwk, twk, rpyis2kc, zv2xfhei, resss, nbzjkp
     *i3, acpios9q, itwk, jwbkl9fp)
      implicit logical (a-z)
      integer irhm4cfa, n, wy1vqfzu, psdvgce3(15), ezlgm2up(*),nef(*),wh
     *ich(*), ges1xpkr(*)
      integer jnxpuym2(*), hnpt1zym(*), fzm1ihwj(*), iz2nbfjc(*), nbzjkp
     *i3(*), acpios9q(*), itwk(*), jwbkl9fp(*)
      double precision he7mqnvy(*),tlgduey8(*),wmat(*),wbkq9zyi(*),dof(*
     *), ub4xioar(*),kispwgx3(*), m0ibglfx(*), s0(wy1vqfzu), beta(*),cov
     *(*),zpcqv3uj, vc6hatuj(*),fasrkub3(*)
      double precision xbig(*), wpuarq2m(*), hjm2ktyr(*), work1(*), wk2(
     *n,wy1vqfzu,3), wkmm(wy1vqfzu,wy1vqfzu,16), work3(*), sgdub(*), bmb
     *(*), ifys6woa(*), mwk(*), twk(*), rpyis2kc(*), zv2xfhei(*), resss
      integer p,q,yzoe1rsp,niter,gtrlbz3e, rutyk8mg, xjc4ywlh, lyma1kwc,
     * dimw, dimu, fbd5yktj,ldk
      integer iter
      integer xs4wtvlg
      integer ayfnwr1v, imk5wjxg, qemj9asg
      irhm4cfa = 0
      imk5wjxg = wy1vqfzu*(wy1vqfzu+1)/2
      p=psdvgce3(2)
      q=psdvgce3(3)
      yzoe1rsp= 0
      if(psdvgce3(4) .eq. 1)then
      yzoe1rsp = 1
      endif
      gtrlbz3e=psdvgce3(6)
      qemj9asg=psdvgce3(7)
      rutyk8mg=psdvgce3(9)
      xjc4ywlh=psdvgce3(10)
      lyma1kwc=psdvgce3(11)
      dimw=psdvgce3(12)
      dimu=psdvgce3(13)
      fbd5yktj = 0
      ldk=psdvgce3(15)
      xs4wtvlg = 1
      if(lyma1kwc .gt. 0)then
      do23304 ayfnwr1v=1,lyma1kwc 
      work1(ayfnwr1v) = dof(ayfnwr1v)
      work1(ayfnwr1v+lyma1kwc) = wbkq9zyi(ayfnwr1v)
      work1(ayfnwr1v+2*lyma1kwc) = dof(ayfnwr1v)
23304 continue
23305 continue
      endif
      iter = 0
23306 if(xs4wtvlg .ne. 0)then
      iter = iter+1
      if(iter .gt. 1)then
      if(lyma1kwc .gt. 0)then
      do23312 ayfnwr1v=1,lyma1kwc 
      if(work1(ayfnwr1v+lyma1kwc).eq.0.0d0 .and. (dabs(work1(ayfnwr1v+2*
     *lyma1kwc)-dof(ayfnwr1v))/dof(ayfnwr1v).gt.0.05d0))then
      work1(ayfnwr1v+2*lyma1kwc) = dof(ayfnwr1v)
      dof(ayfnwr1v)=work1(ayfnwr1v)
      wbkq9zyi(ayfnwr1v)=0.0d0
      else
      work1(ayfnwr1v+2*lyma1kwc) = dof(ayfnwr1v)
      endif
23312 continue
23313 continue
      endif
      endif
      call vbfa1(irhm4cfa,n,wy1vqfzu, he7mqnvy,tlgduey8,wmat,wbkq9zyi,do
     *f, ezlgm2up,nef,which, ub4xioar,kispwgx3,m0ibglfx,s0, beta,cov,zpc
     *qv3uj, vc6hatuj,fasrkub3, qemj9asg,ges1xpkr, xbig, wpuarq2m, hjm2k
     *tyr, jnxpuym2, hnpt1zym, fzm1ihwj(1), fzm1ihwj(1 + imk5wjxg), iz2n
     *bfjc, work1(1+3*lyma1kwc), wkmm, work3, sgdub, bmb, ifys6woa, mwk,
     * twk, rpyis2kc, zv2xfhei, resss, nbzjkpi3, acpios9q, itwk, jwbkl9f
     *p, p,q,yzoe1rsp,niter,gtrlbz3e, wk2(1,1,1), wk2(1,1,2), wk2(1,1,3)
     *, rutyk8mg, xjc4ywlh, lyma1kwc, dimw, dimu, fbd5yktj, ldk)
      if(irhm4cfa .ne. 0)then
      call vcall2(xs4wtvlg,w,y,m0ibglfx,beta,wpuarq2m)
      else
      xs4wtvlg = 0
      endif
      if(xs4wtvlg .ne. 0)then
      qemj9asg=0
      endif
      goto 23306
      endif
23307 continue
      psdvgce3(7) = qemj9asg
      psdvgce3(5) = niter
      psdvgce3(14) = fbd5yktj
      return
      end
      subroutine vbfa1(irhm4cfa,kuzxj1lo,wy1vqfzu, he7mqnvy,tlgduey8,wma
     *t,wbkq9zyi,dof, ezlgm2up,nef,which, ub4xioar,kispwgx3,m0ibglfx,s0,
     * beta,cov,zpcqv3uj, vc6hatuj,fasrkub3, qemj9asg,ges1xpkr, xbig, wp
     *uarq2m, hjm2ktyr, jnxpuym2, hnpt1zym, tgiyxdw1, dufozmt7, iz2nbfjc
     *, work1, wkmm, work3, sgdub, bmb, ifys6woa, mwk, twk, rpyis2kc, zv
     *2xfhei, resss, nbzjkpi3, acpios9q, itwk, jwbkl9fp, p, q, yzoe1rsp,
     * niter, gtrlbz3e, ghz9vuba, oldmat, wk2, rutyk8mg, xjc4ywlh, lyma1
     *kwc, dimw, dimu, fbd5yktj, ldk)
      implicit logical (a-z)
      integer qemj9asg
      integer dufozmt7(*), tgiyxdw1(*)
      integer p, q, yzoe1rsp, niter, gtrlbz3e, rutyk8mg, xjc4ywlh, lyma1
     *kwc, dimw, dimu, fbd5yktj, ldk
      integer irhm4cfa, kuzxj1lo, wy1vqfzu, ezlgm2up(kuzxj1lo,q),nef(q),
     *which(q), ges1xpkr(xjc4ywlh)
      integer jnxpuym2(q), hnpt1zym(q), iz2nbfjc(q), nbzjkpi3(q+1), acpi
     *os9q(q), itwk(*), jwbkl9fp(q+1)
      double precision he7mqnvy(kuzxj1lo,p), tlgduey8(kuzxj1lo,wy1vqfzu)
     *, wmat(kuzxj1lo,dimw), wbkq9zyi(lyma1kwc), dof(lyma1kwc)
      double precision ub4xioar(wy1vqfzu,kuzxj1lo), kispwgx3(kuzxj1lo,ly
     *ma1kwc), m0ibglfx(wy1vqfzu,kuzxj1lo), s0(wy1vqfzu), beta(xjc4ywlh)
     *, cov(kuzxj1lo,lyma1kwc), zpcqv3uj, vc6hatuj(rutyk8mg,xjc4ywlh), f
     *asrkub3(xjc4ywlh)
      double precision xbig(rutyk8mg,xjc4ywlh), wpuarq2m(dimu,kuzxj1lo),
     * hjm2ktyr(wy1vqfzu,lyma1kwc), work1(*), wk2(kuzxj1lo,wy1vqfzu), wk
     *mm(wy1vqfzu,wy1vqfzu,16), work3(*), sgdub(*), bmb(*), ifys6woa(*),
     * mwk(*), twk(*), rpyis2kc(*), zv2xfhei(*), resss
      double precision ghz9vuba(kuzxj1lo,wy1vqfzu), oldmat(kuzxj1lo,wy1v
     *qfzu)
      integer job,info,nefk
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, wg1xifdy
      double precision vo4mtexk, rd9beyfk,ratio, deltaf, z4vrscot,pvofyg
     *8z
      pvofyg8z = 1.0d-7
      job = 101
      info = 1
      if(q .eq. 0)then
      gtrlbz3e = 1
      endif
      if(irhm4cfa .ne. 0)then
      do23324 yq6lorbx=1,xjc4ywlh 
      do23326 ayfnwr1v=1,rutyk8mg 
      vc6hatuj(ayfnwr1v,yq6lorbx)=xbig(ayfnwr1v,yq6lorbx)
23326 continue
23327 continue
23324 continue
23325 continue
      endif
      if(qemj9asg.eq.0)then
      call qpsedg8xf(tgiyxdw1,dufozmt7,wy1vqfzu)
      call mux17f(wpuarq2m, vc6hatuj, wy1vqfzu, xjc4ywlh, kuzxj1lo, wkmm
     *(1,1,1), wkmm(1,1,2), tgiyxdw1, dufozmt7, dimu, rutyk8mg)
      do23330 gp1jxzuh=1,xjc4ywlh 
      ges1xpkr(gp1jxzuh) = gp1jxzuh
23330 continue
23331 continue
      call vqrdca(vc6hatuj,rutyk8mg,rutyk8mg,xjc4ywlh,fasrkub3,ges1xpkr,
     *twk,qemj9asg,pvofyg8z)
      endif
      do23332 yq6lorbx=1,wy1vqfzu 
      do23334 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(yq6lorbx,ayfnwr1v)=0.0d0
23334 continue
23335 continue
      if(q .gt. 0)then
      do23338 gp1jxzuh=1,q 
      if(iz2nbfjc(gp1jxzuh).eq.1)then
      do23342 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) + kispwg
     *x3(ayfnwr1v,hnpt1zym(gp1jxzuh)+yq6lorbx-1)
23342 continue
23343 continue
      else
      do23344 wg1xifdy=1,jnxpuym2(gp1jxzuh) 
      do23346 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) + hjm2kt
     *yr(yq6lorbx,hnpt1zym(gp1jxzuh)+wg1xifdy-1) * kispwgx3(ayfnwr1v,hnp
     *t1zym(gp1jxzuh)+wg1xifdy-1)
23346 continue
23347 continue
23344 continue
23345 continue
      endif
23338 continue
23339 continue
      endif
23332 continue
23333 continue
      niter = 0
      ratio = 1.0d0
23348 if((ratio .gt. zpcqv3uj ) .and. (niter .lt. gtrlbz3e))then
      niter = niter + 1
      deltaf = 0.0d0
      do23350 yq6lorbx=1,wy1vqfzu 
      do23352 ayfnwr1v=1,kuzxj1lo 
      ghz9vuba(ayfnwr1v,yq6lorbx)=tlgduey8(ayfnwr1v,yq6lorbx)-m0ibglfx(y
     *q6lorbx,ayfnwr1v)
23352 continue
23353 continue
23350 continue
23351 continue
      call qpsedg8xf(tgiyxdw1,dufozmt7,wy1vqfzu)
      call mux22f(wpuarq2m,ghz9vuba, twk, dimu,tgiyxdw1,dufozmt7,kuzxj1l
     *o,wy1vqfzu,wkmm)
      call vdqrsl(vc6hatuj,rutyk8mg,rutyk8mg,qemj9asg,fasrkub3, twk, wk2
     *,wk2, beta, wk2,ub4xioar,job,info)
      resss=0.0d0
      do23354 ayfnwr1v=1,kuzxj1lo 
      do23356 yq6lorbx=1,wy1vqfzu 
      vo4mtexk = twk((ayfnwr1v-1)*wy1vqfzu+yq6lorbx) - ub4xioar(yq6lorbx
     *,ayfnwr1v)
      resss = resss + vo4mtexk * vo4mtexk
23356 continue
23357 continue
23354 continue
23355 continue
      call vbksf(wpuarq2m,ub4xioar,wy1vqfzu,kuzxj1lo,wkmm,tgiyxdw1,dufoz
     *mt7,dimu)
      if(q .gt. 0)then
      do23360 gp1jxzuh=1,q 
      do23362 yq6lorbx=1,wy1vqfzu 
      if(iz2nbfjc(gp1jxzuh).eq.1)then
      do23366 ayfnwr1v=1,kuzxj1lo 
      oldmat(ayfnwr1v,yq6lorbx)=kispwgx3(ayfnwr1v,hnpt1zym(gp1jxzuh)+yq6
     *lorbx-1)
      ghz9vuba(ayfnwr1v,yq6lorbx) = tlgduey8(ayfnwr1v,yq6lorbx) - ub4xio
     *ar(yq6lorbx,ayfnwr1v) - m0ibglfx(yq6lorbx,ayfnwr1v) + oldmat(ayfnw
     *r1v,yq6lorbx)
23366 continue
23367 continue
      else
      do23368 ayfnwr1v=1,kuzxj1lo 
      oldmat(ayfnwr1v,yq6lorbx)=0.0d0
      do23370 wg1xifdy=1,jnxpuym2(gp1jxzuh) 
      oldmat(ayfnwr1v,yq6lorbx)=oldmat(ayfnwr1v,yq6lorbx) + hjm2ktyr(yq6
     *lorbx,hnpt1zym(gp1jxzuh)+wg1xifdy-1) * kispwgx3(ayfnwr1v,hnpt1zym(
     *gp1jxzuh)+wg1xifdy-1)
23370 continue
23371 continue
      ghz9vuba(ayfnwr1v,yq6lorbx) = tlgduey8(ayfnwr1v,yq6lorbx) - ub4xio
     *ar(yq6lorbx,ayfnwr1v) - m0ibglfx(yq6lorbx,ayfnwr1v) + oldmat(ayfnw
     *r1v,yq6lorbx)
23368 continue
23369 continue
      endif
23362 continue
23363 continue
      nefk = nef(gp1jxzuh)
      call ewg7qruh(he7mqnvy(1,which(gp1jxzuh)),ghz9vuba,wmat, kuzxj1lo,
     *wy1vqfzu,ezlgm2up(1,gp1jxzuh),nefk, wbkq9zyi(hnpt1zym(gp1jxzuh)), 
     *dof(hnpt1zym(gp1jxzuh)), kispwgx3(1,hnpt1zym(gp1jxzuh)), cov(1,hnp
     *t1zym(gp1jxzuh)), s0, mwk(1), mwk(1+nefk), mwk(1+nefk*(wy1vqfzu+1)
     *), mwk(1+nefk*(2*wy1vqfzu+1)), work1, work3, dimw, fbd5yktj, ldk, 
     *info, yzoe1rsp, sgdub, rpyis2kc(nbzjkpi3(gp1jxzuh)), zv2xfhei(jwbk
     *l9fp(gp1jxzuh)), acpios9q(gp1jxzuh),tgiyxdw1, dufozmt7, bmb, ifys6
     *woa, wkmm, iz2nbfjc(gp1jxzuh),jnxpuym2(gp1jxzuh),itwk, hjm2ktyr(1,
     *hnpt1zym(gp1jxzuh)), twk(1), twk(1+2*jnxpuym2(gp1jxzuh)), twk(1+4*
     *jnxpuym2(gp1jxzuh)), twk(1+(4+nefk)*jnxpuym2(gp1jxzuh)), twk(1+(4+
     *2*nefk)*jnxpuym2(gp1jxzuh)), twk(1+(4+3*nefk)*jnxpuym2(gp1jxzuh)),
     * twk(1+(4+4*nefk)*jnxpuym2(gp1jxzuh)))
      do23372 yq6lorbx=1,wy1vqfzu 
      if(iz2nbfjc(gp1jxzuh).eq.1)then
      do23376 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) + kispwg
     *x3(ayfnwr1v,hnpt1zym(gp1jxzuh)+yq6lorbx-1)
23376 continue
23377 continue
      else
      do23378 wg1xifdy=1,jnxpuym2(gp1jxzuh) 
      do23380 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(yq6lorbx,ayfnwr1v)=m0ibglfx(yq6lorbx,ayfnwr1v) + hjm2ktyr
     *(yq6lorbx,hnpt1zym(gp1jxzuh)+wg1xifdy-1) * kispwgx3(ayfnwr1v,hnpt1
     *zym(gp1jxzuh)+wg1xifdy-1)
23380 continue
23381 continue
23378 continue
23379 continue
      endif
      do23382 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) - oldmat
     *(ayfnwr1v,yq6lorbx)
23382 continue
23383 continue
23372 continue
23373 continue
      do23384 yq6lorbx=1,wy1vqfzu 
      if(iz2nbfjc(gp1jxzuh) .eq. 1)then
      deltaf = deltaf + rd9beyfk(kuzxj1lo,oldmat(1,yq6lorbx),kispwgx3(1,
     *hnpt1zym(gp1jxzuh)+yq6lorbx-1), wmat(1,yq6lorbx))
      else
      do23388 ayfnwr1v=1,kuzxj1lo 
      twk(ayfnwr1v) = 0.0d0
      do23390 wg1xifdy=1,jnxpuym2(gp1jxzuh) 
      twk(ayfnwr1v) = twk(ayfnwr1v) + hjm2ktyr(yq6lorbx,hnpt1zym(gp1jxzu
     *h)+wg1xifdy-1) * kispwgx3(ayfnwr1v,hnpt1zym(gp1jxzuh)+wg1xifdy-1)
23390 continue
23391 continue
23388 continue
23389 continue
      deltaf = deltaf + rd9beyfk(kuzxj1lo, oldmat(1,yq6lorbx), twk, wmat
     *(1,yq6lorbx))
      endif
23384 continue
23385 continue
      do23392 yq6lorbx=1,wy1vqfzu 
      do23394 ayfnwr1v=1,kuzxj1lo 
      ghz9vuba(ayfnwr1v,yq6lorbx)=tlgduey8(ayfnwr1v,yq6lorbx)-m0ibglfx(y
     *q6lorbx,ayfnwr1v)
23394 continue
23395 continue
23392 continue
23393 continue
      call qpsedg8xf(tgiyxdw1,dufozmt7,wy1vqfzu)
      call mux22f(wpuarq2m,ghz9vuba, twk, dimu,tgiyxdw1,dufozmt7,kuzxj1l
     *o,wy1vqfzu,wkmm)
      call vdqrsl(vc6hatuj,rutyk8mg,rutyk8mg,qemj9asg,fasrkub3, twk, wk2
     *,wk2, beta, wk2,ub4xioar,job,info)
      call vbksf(wpuarq2m,ub4xioar,wy1vqfzu,kuzxj1lo,wkmm,tgiyxdw1,dufoz
     *mt7,dimu)
23360 continue
23361 continue
      endif
      if(q .gt. 0)then
      z4vrscot=0.0d0
      do23398 yq6lorbx=1,wy1vqfzu 
      do23400 ayfnwr1v=1,kuzxj1lo 
      z4vrscot = z4vrscot + wmat(ayfnwr1v,yq6lorbx) * m0ibglfx(yq6lorbx,
     *ayfnwr1v)**2
23400 continue
23401 continue
23398 continue
23399 continue
      if(z4vrscot .gt. 0.0d0)then
      ratio = dsqrt(deltaf/z4vrscot)
      else
      ratio = 0.0d0
      endif
      endif
      if(niter .eq. 1)then
      ratio = 1.0d0
      endif
      goto 23348
      endif
23349 continue
      do23406 yq6lorbx=1,xjc4ywlh 
      twk(yq6lorbx)=beta(yq6lorbx)
23406 continue
23407 continue
      do23408 yq6lorbx=1,xjc4ywlh 
      beta(ges1xpkr(yq6lorbx))=twk(yq6lorbx)
23408 continue
23409 continue
      do23410 ayfnwr1v=1,kuzxj1lo 
      do23412 yq6lorbx=1,wy1vqfzu 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) + ub4xio
     *ar(yq6lorbx,ayfnwr1v)
23412 continue
23413 continue
23410 continue
23411 continue
      if((yzoe1rsp .ne. 0) .and. (q .gt. 0))then
      do23416 gp1jxzuh=1,q 
      do23418 wg1xifdy=1,jnxpuym2(gp1jxzuh) 
      call shm8ynte(kuzxj1lo,nef(gp1jxzuh),ezlgm2up(1,gp1jxzuh), cov(1,h
     *npt1zym(gp1jxzuh)+wg1xifdy-1),oldmat)
      do23420 ayfnwr1v=1,kuzxj1lo 
      cov(ayfnwr1v,hnpt1zym(gp1jxzuh)+wg1xifdy-1) = oldmat(ayfnwr1v,1)
23420 continue
23421 continue
23418 continue
23419 continue
23416 continue
23417 continue
      endif
      return
      end
      subroutine x6kanjdh(he7mqnvy, xout, kuzxj1lo, wy1vqfzu)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu
      double precision he7mqnvy(kuzxj1lo), xout(*)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, iptr
      iptr=1
      do23422 yq6lorbx=1,wy1vqfzu 
      do23424 ayfnwr1v=1,kuzxj1lo 
      do23426 gp1jxzuh=1,wy1vqfzu 
      if(yq6lorbx .eq. gp1jxzuh)then
      xout(iptr) = 1.0d0
      else
      xout(iptr) = 0.0d0
      endif
      iptr=iptr+1
23426 continue
23427 continue
23424 continue
23425 continue
23422 continue
23423 continue
      do23430 yq6lorbx=1,wy1vqfzu 
      do23432 ayfnwr1v=1,kuzxj1lo 
      do23434 gp1jxzuh=1,wy1vqfzu 
      if(yq6lorbx .eq. gp1jxzuh)then
      xout(iptr) = he7mqnvy(ayfnwr1v)
      else
      xout(iptr) = 0.0d0
      endif
      iptr=iptr+1
23434 continue
23435 continue
23432 continue
23433 continue
23430 continue
23431 continue
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
      do23438 ayfnwr1v=1,kuzxj1lo 
      work = bhcji9gl(ayfnwr1v) - m0ibglfx(ayfnwr1v)
      rxeqjn0y = rxeqjn0y + po8rwsmy(ayfnwr1v)*work*work
      lm9vcjob = lm9vcjob + po8rwsmy(ayfnwr1v)
23438 continue
23439 continue
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
      do23442 ayfnwr1v=1,kuzxj1lo 
      rxeqjn0y = rxeqjn0y + bhcji9gl(ayfnwr1v) * po8rwsmy(ayfnwr1v)
      lm9vcjob = lm9vcjob + po8rwsmy(ayfnwr1v)
23442 continue
23443 continue
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
      do23446 ayfnwr1v=1,kuzxj1lo 
      hofjnx2e = x(ayfnwr1v)-pygsw6ko
      nsum = nsum + hofjnx2e * (bhcji9gl(ayfnwr1v)-pasjmo8g) * w(ayfnwr1
     *v)
      hofjnx2e = hofjnx2e * hofjnx2e
      q6zdcwxk = q6zdcwxk + hofjnx2e * w(ayfnwr1v)
23446 continue
23447 continue
      eck8vubt = nsum/q6zdcwxk
      interc = pasjmo8g - eck8vubt * pygsw6ko
      do23448 ayfnwr1v=1,kuzxj1lo 
      ub4xioar(ayfnwr1v) = interc + eck8vubt * x(ayfnwr1v)
23448 continue
23449 continue
      bzmd6ftv = interc + eck8vubt * x(1)
      if(yzoe1rsp .ne. 0)then
      do23452 ayfnwr1v=1,kuzxj1lo 
      hofjnx2e = x(ayfnwr1v)-pygsw6ko
      if(w(ayfnwr1v) .gt. 0.0d0)then
      cov(ayfnwr1v,1) = cov(ayfnwr1v,1) - 1.0d0/lm9vcjob - hofjnx2e * ho
     *fjnx2e / q6zdcwxk
      else
      cov(ayfnwr1v,1) = 0.0d0
      endif
23452 continue
23453 continue
      endif
      return
      end
      subroutine shm8ynte(kuzxj1lo, p, ezlgm2up, pygsw6ko, x)
      implicit logical (a-z)
      integer kuzxj1lo, p, ezlgm2up(kuzxj1lo)
      double precision pygsw6ko(p), x(kuzxj1lo)
      integer ayfnwr1v
      do23456 ayfnwr1v=1,kuzxj1lo 
      x(ayfnwr1v) = pygsw6ko(ezlgm2up(ayfnwr1v))
23456 continue
23457 continue
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
      ndk = 40 + dexp(0.25d0 * dlog(kuzxj1lo-40.0d0))
      endif
      else
      ndk = rvy1fpli - 6
      endif
      rvy1fpli = ndk + 6
      do23462 yq6lorbx = 1,3 
      ankcghz2(yq6lorbx) = x(1) 
23462 continue
23463 continue
      do23464 yq6lorbx = 1,ndk 
      ankcghz2(yq6lorbx+3) = x( 1 + (yq6lorbx-1)*(kuzxj1lo-1)/(ndk-1) ) 
23464 continue
23465 continue
      do23466 yq6lorbx = 1,3 
      ankcghz2(ndk+3+yq6lorbx) = x(kuzxj1lo) 
23466 continue
23467 continue
      return
      end
      subroutine pankcghz2l2(ankcghz2, kuzxj1lo, zo8wpibx, tol)
      implicit logical (a-z)
      integer kuzxj1lo, zo8wpibx(kuzxj1lo)
      double precision ankcghz2(kuzxj1lo), tol
      integer ayfnwr1v, cjop5bwm
      do23468 ayfnwr1v=1,4 
      zo8wpibx(ayfnwr1v) = 1
23468 continue
23469 continue
      cjop5bwm = 4
      do23470 ayfnwr1v=5,(kuzxj1lo-4) 
      if((ankcghz2(ayfnwr1v) - ankcghz2(cjop5bwm) .ge. tol) .and. (ankcg
     *hz2(kuzxj1lo) - ankcghz2(ayfnwr1v) .ge. tol))then
      zo8wpibx(ayfnwr1v) = 1
      cjop5bwm = ayfnwr1v
      else
      zo8wpibx(ayfnwr1v) = 0
      endif
23470 continue
23471 continue
      do23474 ayfnwr1v=(kuzxj1lo-3),kuzxj1lo 
      zo8wpibx(ayfnwr1v) = 1
23474 continue
23475 continue
      return
      end
