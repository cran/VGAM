C Output from Public domain Ratfor, version 1.01
      subroutine pnm1or(objzgdk0, lfu2qhid)
      implicit logical (a-z)
      double precision objzgdk0, lfu2qhid
      integer sn
      double precision r1, r2, y, y2, y3, y4, y5, y6, y7
      double precision erf, erfc, z, z2, z3, z4
      double precision sqrt2, sqrtpi, ulimit, p10,p11,p12,p13, q10,q11,q
     *12,q13
      double precision p20,p21,p22,p23,p24,p25,p26,p27
      double precision q20,q21,q22,q23,q24,q25,q26,q27
      double precision p30,p31,p32,p33,p34
      double precision q30,q31,q32,q33,q34
      sqrt2 = 1.414213562373095049d0
      sqrtpi = 1.772453850905516027d0
      ulimit = 20.0d0
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
      if(objzgdk0 .lt. -ulimit)then
      lfu2qhid = 2.753624d-89
      return
      endif
      if(objzgdk0 .gt. ulimit)then
      lfu2qhid = 1.0d0
      return
      endif
      y = objzgdk0 / sqrt2
      if(y .lt. 0.0d0)then
      y = -y
      sn = -1
      else
      sn = 1
      endif
      y2 = y * y
      y4 = y2 * y2
      y6 = y4 * y2
      if(y .lt. 0.46875d0)then
      r1 = p10 + p11 * y2 + p12 * y4 + p13 * y6
      r2 = q10 + q11 * y2 + q12 * y4 + q13 * y6
      erf = y * r1 / r2
      if(sn .eq. 1)then
      lfu2qhid = 0.5d0 + 0.5*erf
      else
      lfu2qhid = 0.5d0 - 0.5*erf
      endif
      else
      if(y .lt. 4.0d0)then
      y3 = y2 * y
      y5 = y4 * y
      y7 = y6 * y
      r1 = p20 + p21 * y + p22 * y2 + p23 * y3 + p24 * y4 + p25 * y5 + p
     *26 * y6 + p27 * y7
      r2 = q20 + q21 * y + q22 * y2 + q23 * y3 + q24 * y4 + q25 * y5 + q
     *26 * y6 + q27 * y7
      erfc = dexp(-y2) * r1 / r2
      if(sn .eq. 1)then
      lfu2qhid = 1.0 - 0.5*erfc
      else
      lfu2qhid = 0.5*erfc
      endif
      else
      z = y4
      z2 = z * z
      z3 = z2 * z
      z4 = z2 * z2
      r1 = p30 + p31 * z + p32 * z2 + p33 * z3 + p34 * z4
      r2 = q30 + q31 * z + q32 * z2 + q33 * z3 + q34 * z4
      erfc = (dexp(-y2)/y) * (1.0 / sqrtpi + r1 / (r2 * y2))
      if(sn .eq. 1)then
      lfu2qhid = 1.0d0 - 0.5*erfc
      else
      lfu2qhid = 0.5*erfc
      endif
      endif
      endif
      return
      end
      subroutine pnm1ow(objzgdk0, lfu2qhid, kuzxj1lo)
      implicit logical (a-z)
      integer kuzxj1lo, ayfnwr1v
      double precision objzgdk0(kuzxj1lo), lfu2qhid(kuzxj1lo)
      do23016 ayfnwr1v=1,kuzxj1lo 
      call pnm1or(objzgdk0(ayfnwr1v), lfu2qhid(ayfnwr1v))
23016 continue
23017 continue
      return
      end
      subroutine n2howibc2a(objzgdk0, i9mwnvqt, lfu2qhid)
      implicit logical (a-z)
      double precision objzgdk0, i9mwnvqt, lfu2qhid
      double precision xd4mybgja
      if(1.0d0 - objzgdk0 .ge. 1.0d0)then
      lfu2qhid = -8.12589d0 / (3.0*dsqrt(i9mwnvqt))
      else
      if(1.0d0 - objzgdk0 .le. 0.0d0)then
      lfu2qhid = 8.12589d0 / (3.0*dsqrt(i9mwnvqt))
      else
      call pnm1or(1.0d0-objzgdk0, xd4mybgja)
      xd4mybgja = xd4mybgja / (3.0*dsqrt(i9mwnvqt))
      lfu2qhid = -3.0d0 * dlog(1.0d0 + xd4mybgja)
      endif
      endif
      return
      end
      subroutine zi8qrpsb(objzgdk0, lfu2qhid)
      implicit logical (a-z)
      double precision objzgdk0, lfu2qhid
      if(1.0d0 - objzgdk0 .ge. 1.0d0)then
      lfu2qhid = -35.0d0
      else
      if(1.0d0 - objzgdk0 .le. 0.0d0)then
      lfu2qhid = 3.542106d0
      else
      lfu2qhid = dlog(-dlog(1.0d0 - objzgdk0))
      endif
      endif
      return
      end
      subroutine g2vwexyk9(objzgdk0, lfu2qhid)
      implicit logical (a-z)
      double precision objzgdk0, lfu2qhid
      if(1.0d0 - objzgdk0 .ge. 1.0d0)then
      lfu2qhid = -34.53958d0
      else
      if(1.0d0 - objzgdk0 .le. 0.0d0)then
      lfu2qhid = 34.53958d0
      else
      lfu2qhid = dlog(objzgdk0 / (1.0d0 - objzgdk0))
      endif
      endif
      return
      end
      subroutine pkc4ejib(w8znmyce, beta, m0ibglfx, kuzxj1lo, wy1vqfzu, 
     *br5ovgcj, xlpjcg3s, vtsou9pz, hj3ftvzu, qfx3vhct, unhycz0e, vm4xjo
     *sb)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu, br5ovgcj, xlpjcg3s, vtsou9pz, hj3ftvzu
     *, qfx3vhct, unhycz0e
      double precision w8znmyce(br5ovgcj,xlpjcg3s), beta(xlpjcg3s), m0ib
     *glfx(wy1vqfzu,kuzxj1lo), vm4xjosb(kuzxj1lo)
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, i1loc, sedf7mxb
      double precision vogkfwt8
      if(vtsou9pz .eq. 1)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      sedf7mxb = 2*hj3ftvzu-1
      do23034 ayfnwr1v=1,kuzxj1lo 
      vogkfwt8 = 0.0d0
      do23036 gp1jxzuh=1,xlpjcg3s 
      vogkfwt8 = vogkfwt8 + w8znmyce(2*ayfnwr1v-1,gp1jxzuh) * beta(gp1jx
     *zuh)
23036 continue
23037 continue
      m0ibglfx(sedf7mxb,ayfnwr1v) = vogkfwt8
23034 continue
23035 continue
      sedf7mxb = 2*hj3ftvzu
      do23038 ayfnwr1v=1,kuzxj1lo 
      vogkfwt8 = 0.0d0
      do23040 gp1jxzuh=1,xlpjcg3s 
      vogkfwt8 = vogkfwt8 + w8znmyce(2*ayfnwr1v ,gp1jxzuh) * beta(gp1jxz
     *uh)
23040 continue
23041 continue
      m0ibglfx(sedf7mxb,ayfnwr1v) = vogkfwt8
23038 continue
23039 continue
      else
      do23042 ayfnwr1v=1,br5ovgcj 
      vogkfwt8 = 0.0d0
      do23044 gp1jxzuh=1,xlpjcg3s 
      vogkfwt8 = vogkfwt8 + w8znmyce(ayfnwr1v,gp1jxzuh) * beta(gp1jxzuh)
23044 continue
23045 continue
      m0ibglfx(hj3ftvzu,ayfnwr1v) = vogkfwt8
23042 continue
23043 continue
      endif
      else
      i1loc = 1
      do23046 ayfnwr1v=1,kuzxj1lo 
      do23048 yq6lorbx=1,wy1vqfzu 
      vogkfwt8 = 0.0d0
      do23050 gp1jxzuh=1,xlpjcg3s 
      vogkfwt8 = vogkfwt8 + w8znmyce(i1loc,gp1jxzuh) * beta(gp1jxzuh)
23050 continue
23051 continue
      i1loc = i1loc + 1
      m0ibglfx(yq6lorbx,ayfnwr1v) = vogkfwt8
23048 continue
23049 continue
23046 continue
23047 continue
      endif
      if(unhycz0e .eq. 1)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23056 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) 
     *+ vm4xjosb(ayfnwr1v)
23056 continue
23057 continue
      else
      do23058 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = m0ibglfx(hj3ftvzu,ayfnwr1v) + vm4xjo
     *sb(ayfnwr1v)
23058 continue
23059 continue
      endif
      endif
      return
      end
      subroutine nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0k
     *ns, qfx3vhct, hj3ftvzu)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, qfx3vhct, hj3ftvzu
      double precision m0ibglfx(wy1vqfzu,kuzxj1lo), t8hwvalr(afpc0kns,ku
     *zxj1lo)
      integer ayfnwr1v, yq6lorbx
      double precision o3jyipdf0
      if(hj3ftvzu .eq. 0)then
      if(qfx3vhct .eq. 1)then
      do23064 ayfnwr1v=1,kuzxj1lo 
      do23066 yq6lorbx=1,wy1vqfzu 
      o3jyipdf0 = dexp(m0ibglfx(yq6lorbx,ayfnwr1v))
      t8hwvalr(yq6lorbx,ayfnwr1v) = o3jyipdf0 / (1.0d0 + o3jyipdf0)
23066 continue
23067 continue
23064 continue
23065 continue
      endif
      if(qfx3vhct .eq. 2)then
      do23070 ayfnwr1v=1,kuzxj1lo 
      do23072 yq6lorbx=1,wy1vqfzu 
      t8hwvalr(yq6lorbx,ayfnwr1v) = dexp(m0ibglfx(yq6lorbx,ayfnwr1v))
23072 continue
23073 continue
23070 continue
23071 continue
      endif
      if(qfx3vhct .eq. 4)then
      do23076 ayfnwr1v=1,kuzxj1lo 
      do23078 yq6lorbx=1,wy1vqfzu 
      t8hwvalr(yq6lorbx,ayfnwr1v) = 1.0d0-dexp(-dexp(m0ibglfx(yq6lorbx,a
     *yfnwr1v)))
23078 continue
23079 continue
23076 continue
23077 continue
      endif
      if(qfx3vhct .eq. 5)then
      do23082 ayfnwr1v=1,kuzxj1lo 
      do23084 yq6lorbx=1,afpc0kns 
      t8hwvalr(yq6lorbx,ayfnwr1v) = dexp(m0ibglfx(2*yq6lorbx-1,ayfnwr1v)
     *)
23084 continue
23085 continue
23082 continue
23083 continue
      endif
      if(qfx3vhct .eq. 3)then
      do23088 ayfnwr1v=1,kuzxj1lo 
      do23090 yq6lorbx=1,afpc0kns 
      t8hwvalr(yq6lorbx,ayfnwr1v) = dexp(m0ibglfx(2*yq6lorbx-1,ayfnwr1v)
     *)
23090 continue
23091 continue
23088 continue
23089 continue
      endif
      if(qfx3vhct .eq. 8)then
      do23094 ayfnwr1v=1,kuzxj1lo 
      do23096 yq6lorbx=1,wy1vqfzu 
      t8hwvalr(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v)
23096 continue
23097 continue
23094 continue
23095 continue
      endif
      else
      if(qfx3vhct .eq. 1)then
      do23100 ayfnwr1v=1,kuzxj1lo 
      o3jyipdf0 = dexp(m0ibglfx(hj3ftvzu,ayfnwr1v))
      t8hwvalr(hj3ftvzu,ayfnwr1v) = o3jyipdf0 / (1.0d0 + o3jyipdf0)
23100 continue
23101 continue
      endif
      if(qfx3vhct .eq. 2)then
      do23104 ayfnwr1v=1,kuzxj1lo 
      t8hwvalr(hj3ftvzu,ayfnwr1v) = dexp(m0ibglfx(hj3ftvzu,ayfnwr1v))
23104 continue
23105 continue
      endif
      if(qfx3vhct .eq. 4)then
      do23108 ayfnwr1v=1,kuzxj1lo 
      t8hwvalr(hj3ftvzu,ayfnwr1v) = 1.0d0 - dexp(-dexp(m0ibglfx(hj3ftvzu
     *,ayfnwr1v)))
23108 continue
23109 continue
      endif
      if(qfx3vhct .eq. 5)then
      do23112 ayfnwr1v=1,kuzxj1lo 
      t8hwvalr(hj3ftvzu,ayfnwr1v) = dexp(m0ibglfx(2*hj3ftvzu-1,ayfnwr1v)
     *)
23112 continue
23113 continue
      endif
      if(qfx3vhct .eq. 3)then
      do23116 ayfnwr1v=1,kuzxj1lo 
      t8hwvalr(hj3ftvzu,ayfnwr1v) = dexp(m0ibglfx(2*hj3ftvzu-1,ayfnwr1v)
     *)
23116 continue
23117 continue
      endif
      if(qfx3vhct .eq. 8)then
      do23120 ayfnwr1v=1,kuzxj1lo 
      t8hwvalr(hj3ftvzu,ayfnwr1v) = m0ibglfx(hj3ftvzu,ayfnwr1v)
23120 continue
23121 continue
      endif
      endif
      return
      end
      subroutine shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, 
     *wy1vqfzu, afpc0kns, dimw, m0ibglfx, dev, hj3ftvzu, n3iasxug, vsoih
     *n1r, cll)
      implicit logical (a-z)
      integer qfx3vhct, kuzxj1lo, wy1vqfzu, afpc0kns, dimw, hj3ftvzu, cl
     *l
      double precision tlgduey8(kuzxj1lo, afpc0kns), wmat(kuzxj1lo, dimw
     *), t8hwvalr(afpc0kns, kuzxj1lo), m0ibglfx(wy1vqfzu,kuzxj1lo), dev,
     * n3iasxug, vsoihn1r
      integer ayfnwr1v, yq6lorbx
      double precision bzmd6ftv, txlvcey5, xd4mybgj, uqnkc6zg, hofjnx2e,
     * smu, afwp5imx, ivqk2ywz, qvd7yktm
      double precision hdqsx7bk, anopu9vi, jtnbu2hz
      logical lbgwvp3q
      bzmd6ftv = 0.0d0
      if(hj3ftvzu .eq. 0)then
      if((qfx3vhct .eq. 1) .or. (qfx3vhct .eq. 4))then
      do23126 yq6lorbx=1,wy1vqfzu 
      do23128 ayfnwr1v=1,kuzxj1lo 
      if(tlgduey8(ayfnwr1v,yq6lorbx) .gt. 0.0d0)then
      ivqk2ywz = tlgduey8(ayfnwr1v,yq6lorbx) * dlog(tlgduey8(ayfnwr1v,yq
     *6lorbx))
      else
      ivqk2ywz = 0.0d0
      endif
      if(tlgduey8(ayfnwr1v,yq6lorbx) .lt. 1.0d0)then
      ivqk2ywz = ivqk2ywz + (1.0d0 - tlgduey8(ayfnwr1v,yq6lorbx)) * dlog
     *(1.0d0 - tlgduey8(ayfnwr1v,yq6lorbx))
      endif
      xd4mybgj = t8hwvalr(yq6lorbx,ayfnwr1v) * (1.0d0 - t8hwvalr(yq6lorb
     *x,ayfnwr1v))
      if(xd4mybgj .lt. n3iasxug)then
      smu = t8hwvalr(yq6lorbx,ayfnwr1v)
      if(smu .lt. n3iasxug)then
      qvd7yktm = tlgduey8(ayfnwr1v,yq6lorbx) * vsoihn1r
      else
      qvd7yktm = tlgduey8(ayfnwr1v,yq6lorbx) * dlog(smu)
      endif
      afwp5imx = 1.0d0 - smu
      if(afwp5imx .lt. n3iasxug)then
      qvd7yktm = qvd7yktm + (1.0d0 - tlgduey8(ayfnwr1v,yq6lorbx)) * vsoi
     *hn1r
      else
      qvd7yktm = qvd7yktm + (1.0d0 - tlgduey8(ayfnwr1v,yq6lorbx)) * dlog
     *(afwp5imx)
      endif
      else
      qvd7yktm = (tlgduey8(ayfnwr1v,yq6lorbx) * dlog(t8hwvalr(yq6lorbx,a
     *yfnwr1v)) + (1.0d0 - tlgduey8(ayfnwr1v,yq6lorbx)) * dlog(1.0d0 - t
     *8hwvalr(yq6lorbx,ayfnwr1v)))
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * (ivqk2ywz - qvd7yktm)
23128 continue
23129 continue
23126 continue
23127 continue
      endif
      if(qfx3vhct .eq. 2)then
      do23142 yq6lorbx=1,wy1vqfzu 
      do23144 ayfnwr1v=1,kuzxj1lo 
      if(tlgduey8(ayfnwr1v,yq6lorbx) .gt. 0.0d0)then
      xd4mybgj = t8hwvalr(yq6lorbx,ayfnwr1v) - tlgduey8(ayfnwr1v,yq6lorb
     *x) + tlgduey8(ayfnwr1v,yq6lorbx) * dlog(tlgduey8(ayfnwr1v,yq6lorbx
     *) / t8hwvalr(yq6lorbx,ayfnwr1v))
      else
      xd4mybgj = t8hwvalr(yq6lorbx,ayfnwr1v) - tlgduey8(ayfnwr1v,yq6lorb
     *x)
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj
23144 continue
23145 continue
23142 continue
23143 continue
      endif
      if(qfx3vhct .eq. 5)then
      do23150 yq6lorbx=1,afpc0kns 
      do23152 ayfnwr1v=1,kuzxj1lo 
      jtnbu2hz = dexp(m0ibglfx(2*yq6lorbx,ayfnwr1v))
      call tldz5ion(jtnbu2hz, uqnkc6zg)
      if(tlgduey8(ayfnwr1v,yq6lorbx) .gt. 0.0d0)then
      xd4mybgj = (jtnbu2hz - 1.0d0) * dlog(tlgduey8(ayfnwr1v,yq6lorbx)) 
     *+ (dlog(jtnbu2hz)-tlgduey8(ayfnwr1v,yq6lorbx) / t8hwvalr(yq6lorbx,
     *ayfnwr1v) - dlog(t8hwvalr(yq6lorbx,ayfnwr1v)) ) * jtnbu2hz - uqnkc
     *6zg
      else
      xd4mybgj = -1000.0d0
      endif
      xd4mybgj = -xd4mybgj
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj
23152 continue
23153 continue
23150 continue
23151 continue
      endif
      if(qfx3vhct .eq. 3)then
      if(cll .eq. 0)then
      anopu9vi = 34.0d0
      do23160 yq6lorbx=1,afpc0kns 
      do23162 ayfnwr1v=1,kuzxj1lo 
      if(m0ibglfx(2*yq6lorbx,ayfnwr1v) .gt. anopu9vi)then
      hdqsx7bk = dexp(anopu9vi)
      lbgwvp3q = .true.
      else
      if(m0ibglfx(2*yq6lorbx,ayfnwr1v) .lt. -anopu9vi)then
      hdqsx7bk = dexp(-anopu9vi)
      lbgwvp3q = .true.
      else
      hdqsx7bk = dexp(m0ibglfx(2*yq6lorbx,ayfnwr1v))
      lbgwvp3q = .false.
      endif
      endif
      if(tlgduey8(ayfnwr1v,yq6lorbx) .lt. 1.0d0)then
      xd4mybgj = 1.0d0
      else
      xd4mybgj = tlgduey8(ayfnwr1v,yq6lorbx)
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * (tlgduey8(ayfnwr1v,yq6lor
     *bx) * dlog(xd4mybgj/t8hwvalr(yq6lorbx,ayfnwr1v)) + (tlgduey8(ayfnw
     *r1v,yq6lorbx) + hdqsx7bk) * dlog((t8hwvalr(yq6lorbx,ayfnwr1v)+hdqs
     *x7bk) / (hdqsx7bk+ tlgduey8(ayfnwr1v,yq6lorbx))))
23162 continue
23163 continue
23160 continue
23161 continue
      else
      anopu9vi = 34.0d0
      do23170 yq6lorbx=1,afpc0kns 
      do23172 ayfnwr1v=1,kuzxj1lo 
      if(m0ibglfx(2*yq6lorbx,ayfnwr1v) .gt. anopu9vi)then
      hdqsx7bk = dexp(anopu9vi)
      lbgwvp3q = .true.
      else
      if(m0ibglfx(2*yq6lorbx,ayfnwr1v) .lt. -anopu9vi)then
      hdqsx7bk = dexp(-anopu9vi)
      lbgwvp3q = .true.
      else
      hdqsx7bk = dexp(m0ibglfx(2*yq6lorbx,ayfnwr1v))
      lbgwvp3q = .false.
      endif
      endif
      if( lbgwvp3q )then
      uqnkc6zg = 0.0d0
      hofjnx2e = 0.0d0
      else
      call tldz5ion(hdqsx7bk + tlgduey8(ayfnwr1v,yq6lorbx), uqnkc6zg)
      call tldz5ion(hdqsx7bk, hofjnx2e)
      endif
      call tldz5ion(1.0d0 + tlgduey8(ayfnwr1v,yq6lorbx), txlvcey5)
      xd4mybgj = hdqsx7bk * dlog(hdqsx7bk / (hdqsx7bk + t8hwvalr(yq6lorb
     *x,ayfnwr1v))) + uqnkc6zg - hofjnx2e - txlvcey5
      if(tlgduey8(ayfnwr1v,yq6lorbx) .gt. 0.0d0)then
      xd4mybgj = xd4mybgj + tlgduey8(ayfnwr1v,yq6lorbx) * dlog(t8hwvalr(
     *yq6lorbx,ayfnwr1v) / (hdqsx7bk + t8hwvalr(yq6lorbx,ayfnwr1v)))
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj
23172 continue
23173 continue
23170 continue
23171 continue
      bzmd6ftv = -bzmd6ftv / 2.0d0
      endif
      endif
      if(qfx3vhct .eq. 8)then
      do23184 yq6lorbx=1,wy1vqfzu 
      do23186 ayfnwr1v=1,kuzxj1lo 
      xd4mybgj = tlgduey8(ayfnwr1v,yq6lorbx) - t8hwvalr(yq6lorbx,ayfnwr1
     *v)
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj**2
23186 continue
23187 continue
23184 continue
23185 continue
      endif
      else
      if((qfx3vhct .eq. 1) .or. (qfx3vhct .eq. 4))then
      do23190 ayfnwr1v=1,kuzxj1lo 
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .gt. 0.0d0)then
      ivqk2ywz = tlgduey8(ayfnwr1v,hj3ftvzu) * dlog(tlgduey8(ayfnwr1v,hj
     *3ftvzu))
      else
      ivqk2ywz = 0.0d0
      endif
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .lt. 1.0d0)then
      ivqk2ywz = ivqk2ywz + (1.0d0 - tlgduey8(ayfnwr1v,hj3ftvzu)) * dlog
     *(1.0d0 - tlgduey8(ayfnwr1v,hj3ftvzu))
      endif
      xd4mybgj = t8hwvalr(hj3ftvzu,ayfnwr1v) * (1.0d0 - t8hwvalr(hj3ftvz
     *u,ayfnwr1v))
      if(xd4mybgj .lt. n3iasxug)then
      smu = t8hwvalr(hj3ftvzu,ayfnwr1v)
      if(smu .lt. n3iasxug)then
      qvd7yktm = tlgduey8(ayfnwr1v,hj3ftvzu) * vsoihn1r
      else
      qvd7yktm = tlgduey8(ayfnwr1v,hj3ftvzu) * dlog(smu)
      endif
      afwp5imx = 1.0d0 - smu
      if(afwp5imx .lt. n3iasxug)then
      qvd7yktm = qvd7yktm + (1.0d0-tlgduey8(ayfnwr1v,hj3ftvzu))*vsoihn1r
      else
      qvd7yktm = qvd7yktm + (1.0d0-tlgduey8(ayfnwr1v,hj3ftvzu))*dlog(afw
     *p5imx)
      endif
      else
      qvd7yktm = (tlgduey8(ayfnwr1v,hj3ftvzu) * dlog(t8hwvalr(hj3ftvzu,a
     *yfnwr1v)) + (1.0d0 - tlgduey8(ayfnwr1v,hj3ftvzu)) * dlog(1.0d0 - t
     *8hwvalr(hj3ftvzu,ayfnwr1v)))
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * (ivqk2ywz - qvd7yktm)
23190 continue
23191 continue
      endif
      if(qfx3vhct .eq. 2)then
      do23204 ayfnwr1v=1,kuzxj1lo 
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .gt. 0.0d0)then
      xd4mybgj = t8hwvalr(hj3ftvzu,ayfnwr1v) - tlgduey8(ayfnwr1v,hj3ftvz
     *u) + tlgduey8(ayfnwr1v,hj3ftvzu) * dlog(tlgduey8(ayfnwr1v,hj3ftvzu
     *) / t8hwvalr(hj3ftvzu,ayfnwr1v))
      else
      xd4mybgj = t8hwvalr(hj3ftvzu,ayfnwr1v) - tlgduey8(ayfnwr1v,hj3ftvz
     *u)
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj
23204 continue
23205 continue
      endif
      if(qfx3vhct .eq. 5)then
      do23210 ayfnwr1v=1,kuzxj1lo 
      jtnbu2hz = dexp(m0ibglfx(2*hj3ftvzu,ayfnwr1v))
      call tldz5ion(jtnbu2hz, uqnkc6zg)
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .gt. 0.0d0)then
      xd4mybgj = (jtnbu2hz - 1.0d0) * dlog(tlgduey8(ayfnwr1v,hj3ftvzu)) 
     *+ jtnbu2hz * (dlog(jtnbu2hz) - tlgduey8(ayfnwr1v,hj3ftvzu) / t8hwv
     *alr(hj3ftvzu,ayfnwr1v) - dlog(t8hwvalr(hj3ftvzu,ayfnwr1v))) - uqnk
     *c6zg
      else
      xd4mybgj = -1000.0d0
      endif
      xd4mybgj = -xd4mybgj
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj
23210 continue
23211 continue
      endif
      if(qfx3vhct .eq. 3)then
      if(cll .eq. 0)then
      anopu9vi = 34.0d0
      do23218 ayfnwr1v=1,kuzxj1lo 
      if(m0ibglfx(2*hj3ftvzu,ayfnwr1v) .gt. anopu9vi)then
      hdqsx7bk = dexp(anopu9vi)
      lbgwvp3q = .true.
      else
      if(m0ibglfx(2*hj3ftvzu,ayfnwr1v) .lt. -anopu9vi)then
      hdqsx7bk = dexp(-anopu9vi)
      lbgwvp3q = .true.
      else
      hdqsx7bk = dexp(m0ibglfx(2*hj3ftvzu,ayfnwr1v))
      lbgwvp3q = .false.
      endif
      endif
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .lt. 1.0d0)then
      xd4mybgj = 1.0d0
      else
      xd4mybgj = tlgduey8(ayfnwr1v,hj3ftvzu)
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * (tlgduey8(ayfnwr1v,hj3ftv
     *zu) * dlog(xd4mybgj/t8hwvalr(hj3ftvzu,ayfnwr1v)) + (tlgduey8(ayfnw
     *r1v,hj3ftvzu)+hdqsx7bk) * dlog((t8hwvalr(hj3ftvzu,ayfnwr1v) + hdqs
     *x7bk) / ( hdqsx7bk+tlgduey8(ayfnwr1v,hj3ftvzu))))
23218 continue
23219 continue
      else
      do23226 ayfnwr1v=1,kuzxj1lo 
      hdqsx7bk = dexp(m0ibglfx(2*hj3ftvzu,ayfnwr1v))
      call tldz5ion(hdqsx7bk + tlgduey8(ayfnwr1v,hj3ftvzu), uqnkc6zg)
      call tldz5ion(hdqsx7bk, hofjnx2e)
      call tldz5ion(1.0d0 + tlgduey8(ayfnwr1v,hj3ftvzu), txlvcey5)
      xd4mybgj = hdqsx7bk * dlog(hdqsx7bk / (hdqsx7bk + t8hwvalr(hj3ftvz
     *u,ayfnwr1v))) + uqnkc6zg - hofjnx2e - txlvcey5
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .gt. 0.0d0)then
      xd4mybgj = xd4mybgj + tlgduey8(ayfnwr1v,hj3ftvzu) * dlog(t8hwvalr(
     *hj3ftvzu,ayfnwr1v) / (hdqsx7bk + t8hwvalr(hj3ftvzu,ayfnwr1v)))
      endif
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj
23226 continue
23227 continue
      bzmd6ftv = -bzmd6ftv / 2.0d0
      endif
      endif
      if(qfx3vhct .eq. 8)then
      do23232 ayfnwr1v=1,kuzxj1lo 
      xd4mybgj = tlgduey8(ayfnwr1v,hj3ftvzu) - t8hwvalr(hj3ftvzu,ayfnwr1
     *v)
      bzmd6ftv = bzmd6ftv + wmat(ayfnwr1v,1) * xd4mybgj**2
23232 continue
23233 continue
      endif
      endif
      dev = 2.0d0 * bzmd6ftv
      return
      end
      subroutine flncwkfq76(lncwkfq7, w8znmyce, kuzxj1lo, br5ovgcj, xwdf
     *5ltg, qfx3vhct)
      implicit logical (a-z)
      integer kuzxj1lo, br5ovgcj, xwdf5ltg, qfx3vhct
      double precision lncwkfq7(kuzxj1lo,xwdf5ltg), w8znmyce(br5ovgcj,*)
      integer ayfnwr1v, sedf7mxb, hpmwnav2
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq.5 ))then
      sedf7mxb = 1
      do23236 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = 1.0d0
      w8znmyce(2*ayfnwr1v, sedf7mxb) = 0.0d0
23236 continue
23237 continue
      sedf7mxb = sedf7mxb + 1
      do23238 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = 0.0d0
      w8znmyce(2*ayfnwr1v, sedf7mxb) = 1.0d0
23238 continue
23239 continue
      sedf7mxb = sedf7mxb + 1
      do23240 hpmwnav2=1,xwdf5ltg 
      do23242 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = lncwkfq7(ayfnwr1v,hpmwnav2)
      w8znmyce(2*ayfnwr1v, sedf7mxb) = 0.0d0
23242 continue
23243 continue
      sedf7mxb = sedf7mxb + 1
23240 continue
23241 continue
      else
      sedf7mxb = 1
      do23244 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(ayfnwr1v,sedf7mxb) = 1.0d0
23244 continue
23245 continue
      sedf7mxb = sedf7mxb + 1
      do23246 hpmwnav2=1,xwdf5ltg 
      do23248 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(ayfnwr1v,sedf7mxb)=lncwkfq7(ayfnwr1v,hpmwnav2)
23248 continue
23249 continue
      sedf7mxb = sedf7mxb + 1
23246 continue
23247 continue
      endif
      return
      end
      subroutine flncwkfq71(lncwkfq7, w8znmyce, kuzxj1lo, xwdf5ltg, qfx3
     *vhct, vm4xjosb, br5ovgcj, xlpjcg3s, hyqwtp6i, tgiyxdw1, dufozmt7, 
     *kifxa0he, p1, unhycz0e)
      implicit logical (a-z)
      integer kuzxj1lo, xwdf5ltg, qfx3vhct, br5ovgcj, xlpjcg3s, hyqwtp6i
     *, tgiyxdw1(hyqwtp6i), dufozmt7(hyqwtp6i), p1, unhycz0e
      double precision lncwkfq7(kuzxj1lo,xwdf5ltg), w8znmyce(br5ovgcj,xl
     *pjcg3s), kifxa0he(kuzxj1lo,p1)
      double precision vm4xjosb(kuzxj1lo)
      integer i0spbklx, ayfnwr1v, sedf7mxb, hpmwnav2
      double precision tad5vhsu, uqnkc6zg
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23252 hpmwnav2=1,xwdf5ltg 
      do23254 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,hpmwnav2) = lncwkfq7(ayfnwr1v,hpmwnav2)
      w8znmyce(2*ayfnwr1v ,hpmwnav2) = 0.0d0
23254 continue
23255 continue
23252 continue
23253 continue
      sedf7mxb = xwdf5ltg + 1
      if(unhycz0e .eq. 0)then
      do23258 i0spbklx=1,hyqwtp6i 
      do23260 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = lncwkfq7(ayfnwr1v,tgiyxdw1(i0spb
     *klx)) * lncwkfq7(ayfnwr1v,dufozmt7(i0spbklx))
      w8znmyce(2*ayfnwr1v ,sedf7mxb) = 0.0d0
23260 continue
23261 continue
      sedf7mxb = sedf7mxb + 1
23258 continue
23259 continue
      else
      do23262 ayfnwr1v=1,kuzxj1lo 
      tad5vhsu = 0.0d0
      do23264 hpmwnav2=1,xwdf5ltg 
      uqnkc6zg = lncwkfq7(ayfnwr1v,hpmwnav2)
      tad5vhsu = tad5vhsu + uqnkc6zg * uqnkc6zg
23264 continue
23265 continue
      vm4xjosb(ayfnwr1v) = -0.50d0 * tad5vhsu
23262 continue
23263 continue
      endif
      else
      do23266 hpmwnav2=1,xwdf5ltg 
      do23268 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(ayfnwr1v,hpmwnav2) = lncwkfq7(ayfnwr1v,hpmwnav2)
23268 continue
23269 continue
23266 continue
23267 continue
      sedf7mxb = xwdf5ltg + 1
      if(unhycz0e .eq. 0)then
      do23272 i0spbklx=1,hyqwtp6i 
      do23274 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(ayfnwr1v,sedf7mxb) = lncwkfq7(ayfnwr1v,tgiyxdw1(i0spbklx)
     *) * lncwkfq7(ayfnwr1v,dufozmt7(i0spbklx))
23274 continue
23275 continue
      sedf7mxb = sedf7mxb + 1
23272 continue
23273 continue
      else
      do23276 ayfnwr1v=1,kuzxj1lo 
      tad5vhsu = 0.0d0
      do23278 hpmwnav2=1,xwdf5ltg 
      uqnkc6zg = lncwkfq7(ayfnwr1v,hpmwnav2)
      tad5vhsu = tad5vhsu + uqnkc6zg * uqnkc6zg
23278 continue
23279 continue
      vm4xjosb(ayfnwr1v) = -0.50d0 * tad5vhsu
23276 continue
23277 continue
      endif
      endif
      if(p1 .gt. 0)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23284 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = 1.0d0
      w8znmyce(2*ayfnwr1v, sedf7mxb) = 0.0d0
23284 continue
23285 continue
      sedf7mxb = sedf7mxb + 1
      do23286 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = 0.0d0
      w8znmyce(2*ayfnwr1v, sedf7mxb) = 1.0d0
23286 continue
23287 continue
      sedf7mxb = sedf7mxb + 1
      if(p1 .gt. 1)then
      do23290 i0spbklx=2,p1 
      do23292 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(2*ayfnwr1v-1,sedf7mxb) = kifxa0he(ayfnwr1v,i0spbklx)
      w8znmyce(2*ayfnwr1v, sedf7mxb) = 0.0d0
23292 continue
23293 continue
      sedf7mxb = sedf7mxb + 1
23290 continue
23291 continue
      endif
      else
      do23294 i0spbklx=1,p1 
      do23296 ayfnwr1v=1,kuzxj1lo 
      w8znmyce(ayfnwr1v,sedf7mxb) = kifxa0he(ayfnwr1v,i0spbklx)
23296 continue
23297 continue
      sedf7mxb = sedf7mxb + 1
23294 continue
23295 continue
      endif
      endif
      return
      end
      subroutine flncwkfq72(lncwkfq7, w8znmyce, kuzxj1lo, wy1vqfzu, br5o
     *vgcj, xwdf5ltg, qfx3vhct, afpc0kns, fmzq7aob, eu3oxvyb, hyqwtp6i, 
     *tgiyxdw1, dufozmt7, unhycz0e, vm4xjosb)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu, br5ovgcj, xwdf5ltg, qfx3vhct, afpc0kns
     *, fmzq7aob, eu3oxvyb, hyqwtp6i, tgiyxdw1(hyqwtp6i), dufozmt7(hyqwt
     *p6i), unhycz0e
      double precision lncwkfq7(kuzxj1lo,xwdf5ltg), w8znmyce(br5ovgcj,*)
     *, vm4xjosb(kuzxj1lo)
      integer i0spbklx, ayfnwr1v, yq6lorbx, gp1jxzuh, ptr, sedf7mxb, hpm
     *wnav2
      double precision uqnkc6zg, tad5vhsu
      do23298 gp1jxzuh=1,eu3oxvyb 
      do23300 ayfnwr1v=1,br5ovgcj 
      w8znmyce(ayfnwr1v,gp1jxzuh) = 0.0d0
23300 continue
23301 continue
23298 continue
23299 continue
      sedf7mxb = 0
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23304 hpmwnav2=1,xwdf5ltg 
      ptr = 1
      do23306 ayfnwr1v=1,kuzxj1lo 
      do23308 yq6lorbx=1,afpc0kns 
      w8znmyce(ptr,sedf7mxb+yq6lorbx) = lncwkfq7(ayfnwr1v,hpmwnav2)
      ptr = ptr + 2
23308 continue
23309 continue
23306 continue
23307 continue
      sedf7mxb = sedf7mxb + afpc0kns
23304 continue
23305 continue
      else
      do23310 hpmwnav2=1,xwdf5ltg 
      ptr = 0
      do23312 ayfnwr1v=1,kuzxj1lo 
      do23314 yq6lorbx=1,wy1vqfzu 
      ptr = ptr + 1
      w8znmyce(ptr,sedf7mxb+yq6lorbx) = lncwkfq7(ayfnwr1v,hpmwnav2)
23314 continue
23315 continue
23312 continue
23313 continue
      sedf7mxb = sedf7mxb + wy1vqfzu
23310 continue
23311 continue
      endif
      if(fmzq7aob .eq. 0)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23320 i0spbklx=1,hyqwtp6i 
      ptr = 1
      do23322 ayfnwr1v=1,kuzxj1lo 
      uqnkc6zg = lncwkfq7(ayfnwr1v,tgiyxdw1(i0spbklx)) * lncwkfq7(ayfnwr
     *1v,dufozmt7(i0spbklx))
      do23324 yq6lorbx=1,afpc0kns 
      w8znmyce(ptr,sedf7mxb+yq6lorbx) = uqnkc6zg
      ptr = ptr + 2
23324 continue
23325 continue
23322 continue
23323 continue
      sedf7mxb = sedf7mxb + afpc0kns
23320 continue
23321 continue
      else
      do23326 i0spbklx=1,hyqwtp6i 
      ptr = 0
      do23328 ayfnwr1v=1,kuzxj1lo 
      uqnkc6zg = lncwkfq7(ayfnwr1v,tgiyxdw1(i0spbklx)) * lncwkfq7(ayfnwr
     *1v,dufozmt7(i0spbklx))
      do23330 yq6lorbx=1,wy1vqfzu 
      ptr = ptr + 1
      w8znmyce(ptr,sedf7mxb+yq6lorbx) = uqnkc6zg
23330 continue
23331 continue
23328 continue
23329 continue
      sedf7mxb = sedf7mxb + wy1vqfzu
23326 continue
23327 continue
      endif
      else
      if(unhycz0e .eq. 1)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23336 ayfnwr1v=1,kuzxj1lo 
      tad5vhsu = 0.0d0
      do23338 hpmwnav2=1,xwdf5ltg 
      uqnkc6zg = lncwkfq7(ayfnwr1v,hpmwnav2)
      tad5vhsu = tad5vhsu + uqnkc6zg * uqnkc6zg
23338 continue
23339 continue
      vm4xjosb(ayfnwr1v) = -0.50d0 * tad5vhsu
23336 continue
23337 continue
      else
      do23340 ayfnwr1v=1,kuzxj1lo 
      tad5vhsu = 0.0d0
      do23342 hpmwnav2=1,xwdf5ltg 
      uqnkc6zg = lncwkfq7(ayfnwr1v,hpmwnav2)
      tad5vhsu = tad5vhsu + uqnkc6zg * uqnkc6zg
23342 continue
23343 continue
      vm4xjosb(ayfnwr1v) = -0.50d0 * tad5vhsu
23340 continue
23341 continue
      endif
      else
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23346 i0spbklx=1,hyqwtp6i 
      ptr = 1
      do23348 ayfnwr1v=1,kuzxj1lo 
      uqnkc6zg = lncwkfq7(ayfnwr1v,tgiyxdw1(i0spbklx)) * lncwkfq7(ayfnwr
     *1v,dufozmt7(i0spbklx))
      do23350 yq6lorbx=1,afpc0kns 
      w8znmyce(ptr,sedf7mxb+i0spbklx) = uqnkc6zg
      ptr = ptr + 2
23350 continue
23351 continue
23348 continue
23349 continue
23346 continue
23347 continue
      sedf7mxb = sedf7mxb + hyqwtp6i
      else
      do23352 i0spbklx=1,hyqwtp6i 
      ptr = 0
      do23354 ayfnwr1v=1,kuzxj1lo 
      uqnkc6zg = lncwkfq7(ayfnwr1v,tgiyxdw1(i0spbklx)) * lncwkfq7(ayfnwr
     *1v,dufozmt7(i0spbklx))
      do23356 yq6lorbx=1,wy1vqfzu 
      ptr = ptr + 1
      w8znmyce(ptr,sedf7mxb+i0spbklx) = uqnkc6zg
23356 continue
23357 continue
23354 continue
23355 continue
23352 continue
23353 continue
      sedf7mxb = sedf7mxb + hyqwtp6i
      endif
      endif
      endif
      return
      end
      subroutine ietam6(tlgduey8, m0ibglfx, y7sdgtqi, kuzxj1lo, wy1vqfzu
     *, afpc0kns, qfx3vhct, hj3ftvzu, wmat, wr0lbopv)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, qfx3vhct, hj3ftvzu, wr0lbopv
      double precision tlgduey8(kuzxj1lo,afpc0kns), m0ibglfx(wy1vqfzu,ku
     *zxj1lo), y7sdgtqi(15)
      double precision wmat(kuzxj1lo,*)
      double precision vogkfwt8, cumw, gyuq8dex, g2vwexykp, qa8ltuhj, kw
     *vo4ury, cpz4fgkx, fguvm9tyi, kinit
      integer ayfnwr1v
      if((qfx3vhct .eq. 1) .or. (qfx3vhct .eq. 4) .or. (qfx3vhct .eq. 3)
     * .or. (qfx3vhct .eq. 5))then
      vogkfwt8 = 0.0d0
      cumw = 0.0d0
      do23360 ayfnwr1v=1,kuzxj1lo 
      vogkfwt8 = vogkfwt8 + tlgduey8(ayfnwr1v,hj3ftvzu) * wmat(ayfnwr1v,
     *1)
      cumw = cumw + wmat(ayfnwr1v,1)
23360 continue
23361 continue
      gyuq8dex = vogkfwt8 / cumw
      endif
      if(qfx3vhct .eq. 1)then
      call g2vwexyk9(gyuq8dex, g2vwexykp)
      do23364 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = g2vwexykp
23364 continue
23365 continue
      endif
      if(qfx3vhct .eq. 2)then
      do23368 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = dlog(tlgduey8(ayfnwr1v,hj3ftvzu) + 0
     *.125d0)
23368 continue
23369 continue
      endif
      if(qfx3vhct .eq. 4)then
      call zi8qrpsb(gyuq8dex, qa8ltuhj)
      do23372 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = qa8ltuhj
23372 continue
23373 continue
      endif
      if(qfx3vhct .eq. 5)then
      if(wr0lbopv .eq. 1)then
      kwvo4ury = dlog(gyuq8dex + 0.03125d0)
      cpz4fgkx = dlog(y7sdgtqi(3+afpc0kns+hj3ftvzu)+0.01d0)
      do23378 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = kwvo4ury
      m0ibglfx(2*hj3ftvzu, ayfnwr1v) = cpz4fgkx
23378 continue
23379 continue
      else
      if(wr0lbopv .eq. 2)then
      kwvo4ury = dlog((6.0/8.0)*gyuq8dex+0.000d0)
      cpz4fgkx = dlog(y7sdgtqi(3+afpc0kns+hj3ftvzu)+0.01d0)
      do23382 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = kwvo4ury
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = cpz4fgkx
23382 continue
23383 continue
      else
      cpz4fgkx = dlog(y7sdgtqi(3+afpc0kns+hj3ftvzu)+0.01d0)
      do23384 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = dlog(tlgduey8(ayfnwr1v,hj3ftvzu)
     * + 0.03125d0)
      m0ibglfx(2*hj3ftvzu, ayfnwr1v) = cpz4fgkx
23384 continue
23385 continue
      endif
      endif
      endif
      if(qfx3vhct .eq. 3)then
      if(wr0lbopv .eq. 1)then
      kwvo4ury = dlog(gyuq8dex + 0.03125d0)
      cpz4fgkx = dlog(y7sdgtqi(3+hj3ftvzu)+0.03125d0)
      do23390 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = kwvo4ury
      m0ibglfx(2*hj3ftvzu,ayfnwr1v) = cpz4fgkx
23390 continue
23391 continue
      else
      if(wr0lbopv .eq. 2)then
      kwvo4ury = dlog(gyuq8dex + 0.03125d0)
      kinit = y7sdgtqi(3+hj3ftvzu)
      cpz4fgkx = dlog(kinit)
      do23394 ayfnwr1v=1,kuzxj1lo 
      fguvm9tyi = tlgduey8(ayfnwr1v,hj3ftvzu) - gyuq8dex
      if(fguvm9tyi .gt. 3.0 * gyuq8dex)then
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = dlog(dsqrt(tlgduey8(ayfnwr1v,hj3
     *ftvzu)))
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = cpz4fgkx
      else
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = kwvo4ury
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = cpz4fgkx
      endif
23394 continue
23395 continue
      else
      if(wr0lbopv .eq. 3)then
      kwvo4ury = dlog(gyuq8dex + 0.03125d0)
      kinit = y7sdgtqi(3+hj3ftvzu)
      cpz4fgkx = dlog(kinit)
      do23400 ayfnwr1v=1,kuzxj1lo 
      fguvm9tyi = tlgduey8(ayfnwr1v,hj3ftvzu) - gyuq8dex
      if(fguvm9tyi .gt. gyuq8dex)then
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = dlog(0.5*(tlgduey8(ayfnwr1v,hj3f
     *tvzu)+gyuq8dex))
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = dlog(kinit / (fguvm9tyi / gyuq8de
     *x))
      else
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .lt. (gyuq8dex / 4.0))then
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = dlog(gyuq8dex / 4.0)
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = cpz4fgkx
      else
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = kwvo4ury
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = cpz4fgkx
      endif
      endif
23400 continue
23401 continue
      else
      cpz4fgkx = dlog(y7sdgtqi(3+hj3ftvzu))
      do23406 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = dlog(tlgduey8(ayfnwr1v,hj3ftvzu)
     * + 0.03125d0)
      m0ibglfx(2*hj3ftvzu, ayfnwr1v) = cpz4fgkx
23406 continue
23407 continue
      endif
      endif
      endif
      endif
      if(qfx3vhct .eq. 8)then
      do23410 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = tlgduey8(ayfnwr1v,hj3ftvzu)
23410 continue
23411 continue
      endif
      return
      end
      subroutine dlgpwe0c(tlgduey8, wmat, m0ibglfx, t8hwvalr, ghz9vuba, 
     *rbne6ouj, wpuarq2m, rsynp1go, n3iasxug, uaf2xgqy, kuzxj1lo, wy1vqf
     *zu, afpc0kns, br5ovgcj, dimu, hj3ftvzu, qfx3vhct, zjkrtol8, unhycz
     *0e, vm4xjosb)
      implicit logical (a-z)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, hj3ftvzu, zj
     *krtol8, unhycz0e
      double precision tlgduey8(kuzxj1lo,afpc0kns), wmat(kuzxj1lo,*), m0
     *ibglfx(wy1vqfzu,kuzxj1lo), t8hwvalr(afpc0kns,kuzxj1lo), vm4xjosb(k
     *uzxj1lo), ghz9vuba(kuzxj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy1vqfzu)
     *, wpuarq2m(dimu,kuzxj1lo), rsynp1go, n3iasxug, uaf2xgqy
      integer ayfnwr1v, qfx3vhct
      double precision xd4mybgja, xd4mybgjb, xd4mybgjc, anopu9vi
      logical lbgwvp3q
      double precision hdqsx7bk, dkdeta, dldk, ux3nadiw, ed2ldk2, n2kers
     *mx, bzmd6ftvmat(1,1), kkmat(1,1), nm0eljqk(1,1)
      integer hbsl0gto, dvhw1ulq, sguwj9ty
      double precision jtnbu2hz, uqnkc6zgd, uqnkc6zgt, dldshape
      double precision fvn3iasxug, xk7dnvei
      integer okobr6tcex
      br5ovgcj = 1
      hbsl0gto = 1
      n2kersmx = 0.990d0
      n2kersmx = 0.995d0
      if(qfx3vhct .eq. 1)then
      do23414 ayfnwr1v=1,kuzxj1lo 
      xd4mybgja = t8hwvalr(hj3ftvzu,ayfnwr1v) * (1.0d0 - t8hwvalr(hj3ftv
     *zu,ayfnwr1v))
      xd4mybgjb = xd4mybgja * wmat(ayfnwr1v,1)
      if(xd4mybgja .lt. n3iasxug)then
      xd4mybgja = n3iasxug
      endif
      if(xd4mybgjb .lt. n3iasxug)then
      xd4mybgjb = n3iasxug
      wpuarq2m(hj3ftvzu,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(hj3ftvzu,ayfnwr1v) = dsqrt(xd4mybgjb)
      endif
      rbne6ouj(ayfnwr1v,hj3ftvzu) = xd4mybgjb
      ghz9vuba(ayfnwr1v,hj3ftvzu) = m0ibglfx(hj3ftvzu,ayfnwr1v) + (tlgdu
     *ey8(ayfnwr1v,hj3ftvzu)-t8hwvalr(hj3ftvzu,ayfnwr1v)) / xd4mybgja
23414 continue
23415 continue
      endif
      if(qfx3vhct .eq. 2)then
      do23422 ayfnwr1v=1,kuzxj1lo 
      xd4mybgja = t8hwvalr(hj3ftvzu,ayfnwr1v)
      xd4mybgjb = xd4mybgja * wmat(ayfnwr1v,1)
      if(xd4mybgjb .lt. n3iasxug)then
      xd4mybgjb = n3iasxug
      wpuarq2m(hj3ftvzu,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(hj3ftvzu,ayfnwr1v) = dsqrt(xd4mybgjb)
      endif
      rbne6ouj(ayfnwr1v,hj3ftvzu) = xd4mybgjb
      if(tlgduey8(ayfnwr1v,hj3ftvzu) .gt. 0.0d0)then
      xd4mybgjc = xd4mybgja
      if(xd4mybgjc .lt. n3iasxug)then
      xd4mybgjc = n3iasxug
      endif
      ghz9vuba(ayfnwr1v,hj3ftvzu) = m0ibglfx(hj3ftvzu,ayfnwr1v) + (tlgdu
     *ey8(ayfnwr1v,hj3ftvzu)-xd4mybgjc)/xd4mybgjc
      else
      ghz9vuba(ayfnwr1v,hj3ftvzu) = m0ibglfx(hj3ftvzu,ayfnwr1v) - 1.0d0
      endif
23422 continue
23423 continue
      endif
      if(qfx3vhct .eq. 4)then
      do23432 ayfnwr1v=1,kuzxj1lo 
      if((t8hwvalr(hj3ftvzu,ayfnwr1v) .lt. n3iasxug) .or. (t8hwvalr(hj3f
     *tvzu,ayfnwr1v) .gt. 1.0d0 - n3iasxug))then
      xd4mybgja = n3iasxug
      xd4mybgjb = xd4mybgja * wmat(ayfnwr1v,1)
      if(xd4mybgjb .lt. n3iasxug)then
      xd4mybgjb = n3iasxug
      wpuarq2m(hj3ftvzu,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(hj3ftvzu,ayfnwr1v) = dsqrt(xd4mybgjb)
      endif
      rbne6ouj(ayfnwr1v,hj3ftvzu) = xd4mybgjb
      ghz9vuba(ayfnwr1v,hj3ftvzu) = m0ibglfx(hj3ftvzu,ayfnwr1v) + (tlgdu
     *ey8(ayfnwr1v,hj3ftvzu)-t8hwvalr(hj3ftvzu,ayfnwr1v)) / xd4mybgja
      else
      xd4mybgja = -(1.0d0 - t8hwvalr(hj3ftvzu,ayfnwr1v)) * dlog(1.0d0 - 
     *t8hwvalr(hj3ftvzu,ayfnwr1v))
      if(xd4mybgja .lt. n3iasxug)then
      xd4mybgja = n3iasxug
      endif
      xd4mybgjb = -xd4mybgja * wmat(ayfnwr1v,1) * dlog(1.0d0 - t8hwvalr(
     *hj3ftvzu,ayfnwr1v)) / t8hwvalr(hj3ftvzu,ayfnwr1v)
      if(xd4mybgjb .lt. n3iasxug)then
      xd4mybgjb = n3iasxug
      endif
      rbne6ouj(ayfnwr1v,hj3ftvzu) = xd4mybgjb
      wpuarq2m(hj3ftvzu,ayfnwr1v) = dsqrt(xd4mybgjb)
      ghz9vuba(ayfnwr1v,hj3ftvzu) = m0ibglfx(hj3ftvzu,ayfnwr1v) + (tlgdu
     *ey8(ayfnwr1v,hj3ftvzu)-t8hwvalr(hj3ftvzu,ayfnwr1v)) / xd4mybgja
      endif
23432 continue
23433 continue
      endif
      if(qfx3vhct .eq. 5)then
      fvn3iasxug = 1.0d-20
      anopu9vi = 34.0d0
      do23444 ayfnwr1v=1,kuzxj1lo 
      if(m0ibglfx(2*hj3ftvzu,ayfnwr1v) .gt. anopu9vi)then
      jtnbu2hz = dexp(anopu9vi)
      lbgwvp3q = .true.
      else
      if(m0ibglfx(2*hj3ftvzu,ayfnwr1v) .lt. -anopu9vi)then
      jtnbu2hz = dexp(-anopu9vi)
      lbgwvp3q = .true.
      else
      jtnbu2hz = dexp(m0ibglfx(2*hj3ftvzu,ayfnwr1v))
      lbgwvp3q = .false.
      endif
      endif
      call vdgam1(jtnbu2hz, uqnkc6zgd, okobr6tcex)
      if(okobr6tcex .ne. 1)then
      call intpr("error in dlgpwe0c okobr6tcex 1: ",-1,okobr6tcex,1)
      endif
      xk7dnvei = t8hwvalr(hj3ftvzu,ayfnwr1v)
      if(xk7dnvei .lt. fvn3iasxug)then
      xk7dnvei = fvn3iasxug
      endif
      dldshape = dlog(tlgduey8(ayfnwr1v,hj3ftvzu)) + dlog(jtnbu2hz) - dl
     *og(xk7dnvei) + 1.0d0 - uqnkc6zgd - tlgduey8(ayfnwr1v,hj3ftvzu) / x
     *k7dnvei
      call vtgam1(jtnbu2hz, uqnkc6zgt, okobr6tcex)
      if(okobr6tcex .ne. 1)then
      call intpr("error in dlgpwe0c okobr6tcex 2: ",-1,okobr6tcex,1)
      endif
      rbne6ouj(ayfnwr1v,2*hj3ftvzu-1) = wmat(ayfnwr1v,1) * jtnbu2hz
      xd4mybgja = jtnbu2hz * uqnkc6zgt - 1.0d0
      rbne6ouj(ayfnwr1v,2*hj3ftvzu ) = wmat(ayfnwr1v,1) * jtnbu2hz * xd4
     *mybgja
      if(rbne6ouj(ayfnwr1v,2*hj3ftvzu-1) .lt. n3iasxug)then
      rbne6ouj(ayfnwr1v,2*hj3ftvzu-1) = n3iasxug
      wpuarq2m(2*hj3ftvzu-1,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(2*hj3ftvzu-1,ayfnwr1v) = dsqrt(rbne6ouj(ayfnwr1v,2*hj3ftv
     *zu-1))
      endif
      if(rbne6ouj(ayfnwr1v,2*hj3ftvzu) .lt. n3iasxug)then
      rbne6ouj(ayfnwr1v,2*hj3ftvzu) = n3iasxug
      wpuarq2m(2*hj3ftvzu,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(2*hj3ftvzu,ayfnwr1v) = dsqrt(rbne6ouj(ayfnwr1v,2*hj3ftvzu
     *))
      endif
      if(xd4mybgja .lt. fvn3iasxug)then
      xd4mybgja = fvn3iasxug
      endif
      ghz9vuba(ayfnwr1v,2*hj3ftvzu-1) = m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) 
     *+ tlgduey8(ayfnwr1v,hj3ftvzu) / xk7dnvei - 1.0d0
      ghz9vuba(ayfnwr1v,2*hj3ftvzu ) = m0ibglfx(2*hj3ftvzu ,ayfnwr1v) + 
     *dldshape / xd4mybgja
23444 continue
23445 continue
      endif
      if(qfx3vhct .eq. 3)then
      anopu9vi = 34.0d0
      fvn3iasxug = 1.0d-20
      do23464 ayfnwr1v=1,kuzxj1lo 
      if(m0ibglfx(2*hj3ftvzu,ayfnwr1v) .gt. anopu9vi)then
      hdqsx7bk = dexp(anopu9vi)
      lbgwvp3q = .true.
      else
      if(m0ibglfx(2*hj3ftvzu,ayfnwr1v) .lt. -anopu9vi)then
      hdqsx7bk = dexp(-anopu9vi)
      lbgwvp3q = .true.
      else
      hdqsx7bk = dexp(m0ibglfx(2*hj3ftvzu,ayfnwr1v))
      lbgwvp3q = .false.
      endif
      endif
      xk7dnvei = t8hwvalr(hj3ftvzu,ayfnwr1v)
      if(xk7dnvei .lt. fvn3iasxug)then
      xk7dnvei = fvn3iasxug
      endif
      call vdgam1(tlgduey8(ayfnwr1v,hj3ftvzu) + hdqsx7bk, xd4mybgja, oko
     *br6tcex)
      if(okobr6tcex .ne. 1)then
      endif
      call vdgam1(hdqsx7bk, xd4mybgjb, okobr6tcex)
      if(okobr6tcex .ne. 1)then
      endif
      dldk = xd4mybgja - xd4mybgjb - (tlgduey8(ayfnwr1v,hj3ftvzu) + hdqs
     *x7bk) / (xk7dnvei + hdqsx7bk) + 1.0d0 + dlog(hdqsx7bk / (xk7dnvei 
     *+ hdqsx7bk))
      dkdeta = hdqsx7bk
      kkmat(1,1) = hdqsx7bk
      nm0eljqk(1,1) = xk7dnvei
      sguwj9ty = 5000
      call enbin9(bzmd6ftvmat, kkmat, nm0eljqk, n2kersmx, hbsl0gto, dvhw
     *1ulq, hbsl0gto, ux3nadiw, rsynp1go, sguwj9ty)
      if(dvhw1ulq .ne. 1)then
      zjkrtol8 = 5
      return
      endif
      ed2ldk2 = -bzmd6ftvmat(1,1) - 1.0d0 / hdqsx7bk + 1.0d0 / (hdqsx7bk
     * + xk7dnvei)
      rbne6ouj(ayfnwr1v,2*hj3ftvzu-1) = wmat(ayfnwr1v,1) * xk7dnvei * hd
     *qsx7bk / (xk7dnvei + hdqsx7bk)
      rbne6ouj(ayfnwr1v,2*hj3ftvzu ) = wmat(ayfnwr1v,1) * hdqsx7bk * (-b
     *zmd6ftvmat(1,1)*hdqsx7bk - 1.0d0 + hdqsx7bk / (hdqsx7bk + xk7dnvei
     *))
      if(rbne6ouj(ayfnwr1v,2*hj3ftvzu-1) .lt. n3iasxug)then
      rbne6ouj(ayfnwr1v,2*hj3ftvzu-1) = n3iasxug
      wpuarq2m(2*hj3ftvzu-1,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(2*hj3ftvzu-1,ayfnwr1v) = dsqrt(rbne6ouj(ayfnwr1v,2*hj3ftv
     *zu-1))
      endif
      if(rbne6ouj(ayfnwr1v,2*hj3ftvzu) .lt. n3iasxug)then
      rbne6ouj(ayfnwr1v,2*hj3ftvzu) = n3iasxug
      wpuarq2m(2*hj3ftvzu,ayfnwr1v) = uaf2xgqy
      else
      wpuarq2m(2*hj3ftvzu,ayfnwr1v) = dsqrt(rbne6ouj(ayfnwr1v,2*hj3ftvzu
     *))
      endif
      ghz9vuba(ayfnwr1v,2*hj3ftvzu-1) = m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) 
     *+ tlgduey8(ayfnwr1v,hj3ftvzu) / xk7dnvei - 1.0d0
      ghz9vuba(ayfnwr1v,2*hj3ftvzu ) = m0ibglfx(2*hj3ftvzu ,ayfnwr1v) + 
     *dldk / (dkdeta * ed2ldk2)
23464 continue
23465 continue
      endif
      if(qfx3vhct .eq. 8)then
      do23484 ayfnwr1v=1,kuzxj1lo 
      rbne6ouj(ayfnwr1v,hj3ftvzu) = wmat(ayfnwr1v,1)
      wpuarq2m(hj3ftvzu,ayfnwr1v) = dsqrt(rbne6ouj(ayfnwr1v,hj3ftvzu))
      ghz9vuba(ayfnwr1v,hj3ftvzu) = tlgduey8(ayfnwr1v,hj3ftvzu)
23484 continue
23485 continue
      endif
      if(unhycz0e .eq. 1)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23490 ayfnwr1v=1,kuzxj1lo 
      ghz9vuba(ayfnwr1v,2*hj3ftvzu-1) = ghz9vuba(ayfnwr1v,2*hj3ftvzu-1) 
     *- vm4xjosb(ayfnwr1v)
23490 continue
23491 continue
      else
      do23492 ayfnwr1v=1,kuzxj1lo 
      ghz9vuba(ayfnwr1v,hj3ftvzu) = ghz9vuba(ayfnwr1v,hj3ftvzu) - vm4xjo
     *sb(ayfnwr1v)
23492 continue
23493 continue
      endif
      endif
      return
      end
      subroutine cqo2f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4
     *xjosb, t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj,
     * fasrkub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu,
     * zjkrtol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm
     *, y7sdgtqi)
      implicit logical (a-z)
      integer xui7hqwl(18), tgiyxdw1(*), dufozmt7(*)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, ge
     *s1xpkr(*)
      double precision lncwkfq7(kuzxj1lo,*), tlgduey8(kuzxj1lo,afpc0kns)
     *, kifxa0he(kuzxj1lo,*), wmat(kuzxj1lo,*), m0ibglfx(wy1vqfzu,kuzxj1
     *lo), vm4xjosb(kuzxj1lo), t8hwvalr(afpc0kns,kuzxj1lo)
      double precision ghz9vuba(kuzxj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy
     *1vqfzu), wpuarq2m(dimu,kuzxj1lo), w8znmyce(br5ovgcj,*)
      double precision vc6hatuj(br5ovgcj,*), fasrkub3(*), tlq9wpes, beta
     *(*), y7sdgtqi(*)
      double precision twk(wy1vqfzu,kuzxj1lo,2), wkmm(wy1vqfzu * (wy1vqf
     *zu + 1))
      integer ayfnwr1v, yq6lorbx, gp1jxzuh, hyqwtp6i, ptr, i1loc, i2, iz
     *ero0, iter, fmzq7aob, xwdf5ltg, dimw, f7svlajr, qfx3vhct, c5aesxku
     *l
      integer job, info, qemj9asg, xlpjcg3s, eu3oxvyb, vtsou9pz, unhycz0
     *e, zaupqv9b
      integer hbsl0gto, wr0lbopv
      double precision rpto5qwb, n3iasxug, pvofyg8z, wiptsjx8, uylxqtc7,
     * bh2vgiay, uaf2xgqy, vsoihn1r, rsynp1go
      double precision hmayv1xt1, hmayv1xt2
      integer x1jrewny
      hbsl0gto = 1
      x1jrewny = 0
      kifxa0he(1,1) = 1
      wkmm(1) = 0.0d0
      xwdf5ltg = xui7hqwl(1)
      fmzq7aob = xui7hqwl(2)
      xlpjcg3s = xui7hqwl(3)
      dimw = xui7hqwl(4)
      f7svlajr = xui7hqwl(5)
      qfx3vhct = xui7hqwl(6)
      c5aesxkul = xui7hqwl(7)
      xui7hqwl(9) = 0
      eu3oxvyb = xui7hqwl(11)
      vtsou9pz = xui7hqwl(12)
      unhycz0e = xui7hqwl(14)
      zaupqv9b = xui7hqwl(15)
      wr0lbopv = xui7hqwl(18)
      n3iasxug = y7sdgtqi(1)
      uaf2xgqy = dsqrt(n3iasxug)
      if((qfx3vhct .eq. 1) .or. (qfx3vhct .eq. 4))then
      vsoihn1r = dlog(n3iasxug)
      endif
      bh2vgiay = y7sdgtqi(2)
      rsynp1go = y7sdgtqi(3)
      uylxqtc7 = 0.0d0
      izero0 = 0
      zjkrtol8 = 1
      call qpsedg8xf(tgiyxdw1, dufozmt7, xwdf5ltg)
      hyqwtp6i = xwdf5ltg * (xwdf5ltg+1) / 2
      call flncwkfq72(lncwkfq7, w8znmyce, kuzxj1lo, wy1vqfzu, br5ovgcj, 
     *xwdf5ltg, qfx3vhct, afpc0kns, fmzq7aob, eu3oxvyb, hyqwtp6i, tgiyxd
     *w1, dufozmt7, unhycz0e, vm4xjosb)
653   hmayv1xt2 = 1.0d0
      if(f7svlajr .eq. 0)then
      do23498 yq6lorbx=1,afpc0kns 
      call ietam6(tlgduey8, m0ibglfx, y7sdgtqi, kuzxj1lo, wy1vqfzu, afpc
     *0kns, qfx3vhct, yq6lorbx, wmat, wr0lbopv)
23498 continue
23499 continue
      else
      if(f7svlajr .eq. 2)then
      call pkc4ejib(w8znmyce, beta, m0ibglfx, kuzxj1lo, wy1vqfzu, br5ovg
     *cj, xlpjcg3s, vtsou9pz, izero0, qfx3vhct, unhycz0e, vm4xjosb)
      endif
      endif
      call nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0kns, qf
     *x3vhct, izero0)
      if(f7svlajr .eq. 2)then
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, rpto5qwb, izero0, n3iasxug, vsoihn1r
     *, hbsl0gto)
      else
      rpto5qwb = -1.0d0
      endif
      do23504 iter=1,c5aesxkul 
      do23506 yq6lorbx=1,afpc0kns 
      call dlgpwe0c(tlgduey8, wmat, m0ibglfx, t8hwvalr, ghz9vuba, rbne6o
     *uj, wpuarq2m, rsynp1go, n3iasxug, uaf2xgqy, kuzxj1lo, wy1vqfzu, af
     *pc0kns, br5ovgcj, dimu, yq6lorbx, qfx3vhct, zjkrtol8, unhycz0e, vm
     *4xjosb)
23506 continue
23507 continue
      do23508 yq6lorbx=1,xlpjcg3s 
      do23510 ayfnwr1v=1,br5ovgcj 
      vc6hatuj(ayfnwr1v,yq6lorbx) = w8znmyce(ayfnwr1v,yq6lorbx)
23510 continue
23511 continue
23508 continue
23509 continue
      do23512 yq6lorbx=1,xlpjcg3s 
      ptr = 1
      do23514 i1loc=1,kuzxj1lo 
      do23516 i2=1,wy1vqfzu 
      vc6hatuj(ptr,yq6lorbx) = wpuarq2m(i2,i1loc) * vc6hatuj(ptr,yq6lorb
     *x)
      ptr = ptr + 1
23516 continue
23517 continue
23514 continue
23515 continue
23512 continue
23513 continue
      do23518 gp1jxzuh=1,xlpjcg3s 
      ges1xpkr(gp1jxzuh) = gp1jxzuh
23518 continue
23519 continue
      pvofyg8z = 1.0d-7
      call vqrdca(vc6hatuj,br5ovgcj,br5ovgcj,xlpjcg3s,fasrkub3,ges1xpkr,
     *twk,qemj9asg,pvofyg8z)
      if(qemj9asg .ne. xlpjcg3s)then
      zjkrtol8 = 2
      return
      endif
      do23522 ayfnwr1v=1,kuzxj1lo 
      do23524 yq6lorbx=1,wy1vqfzu 
      twk(yq6lorbx,ayfnwr1v,1) = wpuarq2m(yq6lorbx,ayfnwr1v) * ghz9vuba(
     *ayfnwr1v,yq6lorbx)
23524 continue
23525 continue
23522 continue
23523 continue
      job = 101
      call vdqrsl(vc6hatuj,br5ovgcj,br5ovgcj,qemj9asg,fasrkub3, twk, uyl
     *xqtc7, twk(1,1,2), beta, uylxqtc7,m0ibglfx,job,info)
      do23526 ayfnwr1v=1,kuzxj1lo 
      do23528 yq6lorbx=1,wy1vqfzu 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) / wpuarq
     *2m(yq6lorbx,ayfnwr1v)
23528 continue
23529 continue
23526 continue
23527 continue
      if(unhycz0e .eq. 1)then
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23534 ayfnwr1v=1,kuzxj1lo 
      do23536 yq6lorbx=1,afpc0kns 
      m0ibglfx(2*yq6lorbx-1,ayfnwr1v) = m0ibglfx(2*yq6lorbx-1,ayfnwr1v) 
     *+ vm4xjosb(ayfnwr1v)
23536 continue
23537 continue
23534 continue
23535 continue
      else
      do23538 ayfnwr1v=1,kuzxj1lo 
      do23540 yq6lorbx=1,wy1vqfzu 
      m0ibglfx(yq6lorbx,ayfnwr1v) = m0ibglfx(yq6lorbx,ayfnwr1v) + vm4xjo
     *sb(ayfnwr1v)
23540 continue
23541 continue
23538 continue
23539 continue
      endif
      endif
      call nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0kns, qf
     *x3vhct, izero0)
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, tlq9wpes,izero0,n3iasxug,vsoihn1r, h
     *bsl0gto)
      wiptsjx8 = dabs(tlq9wpes - rpto5qwb) / (1.0d0 + dabs(tlq9wpes))
      if(wiptsjx8 .lt. bh2vgiay)then
      zjkrtol8 = 0
      xui7hqwl(8) = iter
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, tlq9wpes,izero0,n3iasxug,vsoihn1r, i
     *zero0)
      endif
      x1jrewny = 1
      goto 20097
      else
      rpto5qwb = tlq9wpes
      x1jrewny = 0
      endif
23504 continue
23505 continue
20097 hmayv1xt1 = 0.0d0
      if(x1jrewny .eq. 1)then
      return
      endif
      if(f7svlajr .eq. 1 .or. f7svlajr .eq. 2)then
      f7svlajr = 0
      xui7hqwl(9) = 1
      goto 653
      endif
      zjkrtol8 = 3
      return
      end
      subroutine cqo1f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4
     *xjosb, t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj,
     * fasrkub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu,
     * zjkrtol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm
     *, y7sdgtqi)
      implicit logical (a-z)
      integer xui7hqwl(18), tgiyxdw1(*), dufozmt7(*)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, ge
     *s1xpkr(*)
      double precision lncwkfq7(kuzxj1lo,*), tlgduey8(kuzxj1lo,afpc0kns)
     *, wmat(kuzxj1lo,*), m0ibglfx(wy1vqfzu,kuzxj1lo), vm4xjosb(kuzxj1lo
     *), t8hwvalr(afpc0kns,kuzxj1lo), kifxa0he(kuzxj1lo,*), ghz9vuba(kuz
     *xj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy1vqfzu), wpuarq2m(dimu,kuzxj1
     *lo), w8znmyce(br5ovgcj,*)
      double precision vc6hatuj(br5ovgcj,*), fasrkub3(*), tlq9wpes, beta
     *(*), y7sdgtqi(*)
      double precision twk(br5ovgcj,3), wkmm(wy1vqfzu*(wy1vqfzu+1))
      integer ayfnwr1v, hj3ftvzu, hyqwtp6i, izero0, iter, fmzq7aob, unhy
     *cz0e, xwdf5ltg, dimw, f7svlajr, qfx3vhct, c5aesxkul
      integer job, info, qemj9asg, xlpjcg3s, vtsou9pz, zaupqv9b
      integer hbsl0gto, p1, wr0lbopv
      double precision rpto5qwb, n3iasxug, pvofyg8z, wiptsjx8, uylxqtc7,
     * bh2vgiay, uaf2xgqy, vsoihn1r, rsynp1go
      integer gp1jxzuh
      double precision aqg1vdmo, hmayv1xt
      aqg1vdmo = 0.0d0
      hbsl0gto = 1
      wkmm(1) = 1.0d0
      call intpr("entering cqo1f hbsl0gto ------------------------------
     *-: ",-1,hbsl0gto,1)
      call intpr("in cqo1f afpc0kns: ",-1,afpc0kns,1)
      xwdf5ltg = xui7hqwl(1)
      fmzq7aob = xui7hqwl(2)
      xlpjcg3s = xui7hqwl(3)
      dimw = xui7hqwl(4)
      f7svlajr = xui7hqwl(5)
      qfx3vhct = xui7hqwl(6)
      c5aesxkul = xui7hqwl(7)
      xui7hqwl(9) = 0
      vtsou9pz = xui7hqwl(12)
      if(vtsou9pz .ne. 1)then
      zjkrtol8 = 4
      return
      endif
      unhycz0e = xui7hqwl(14)
      zaupqv9b = xui7hqwl(15)
      p1 = xui7hqwl(16)
      wr0lbopv = xui7hqwl(18)
      call intpr("Entry to cqo1f: f7svlajr ",-1,f7svlajr,1)
      n3iasxug = y7sdgtqi(1)
      uaf2xgqy = dsqrt(n3iasxug)
      if((qfx3vhct .eq. 1) .or. (qfx3vhct .eq. 4))then
      vsoihn1r = dlog(n3iasxug)
      endif
      bh2vgiay = y7sdgtqi(2)
      rsynp1go = y7sdgtqi(3)
      uylxqtc7 = 0.0d0
      izero0 = 0
      zjkrtol8 = 1
      call qpsedg8xf(tgiyxdw1, dufozmt7, xwdf5ltg)
      hyqwtp6i = xwdf5ltg * (xwdf5ltg+1) / 2
      call flncwkfq71(lncwkfq7, w8znmyce, kuzxj1lo, xwdf5ltg, qfx3vhct, 
     *vm4xjosb, br5ovgcj, xlpjcg3s, hyqwtp6i, tgiyxdw1, dufozmt7, kifxa0
     *he, p1, unhycz0e)
      call dblepr("cqo1f: vm4xjosb()",-1,vm4xjosb,kuzxj1lo)
      call dblepr("cqo1f: w8znmyce(,)",-1,w8znmyce,br5ovgcj*xlpjcg3s)
      call dblepr("cqo1f: wmat(,1)",-1,wmat(1,1),kuzxj1lo)
      do23554 hj3ftvzu=1,afpc0kns 
      call intpr("cqo1f: hj3ftvzu======================: ",-1,hj3ftvzu,1
     *)
653   hmayv1xt = 1.0d0
      if(f7svlajr .eq. 0)then
      call intpr("cqo1f: calling ietam6 ",-1,hj3ftvzu,1)
      call ietam6(tlgduey8, m0ibglfx, y7sdgtqi, kuzxj1lo, wy1vqfzu, afpc
     *0kns, qfx3vhct, hj3ftvzu, wmat, wr0lbopv)
      else
      if(f7svlajr .eq. 2)then
      call intpr("cqo1f: calling pkc4ejib; vtsou9pz== ",-1,vtsou9pz,1)
      call pkc4ejib(w8znmyce, beta(1+(hj3ftvzu-1)*xlpjcg3s), m0ibglfx, k
     *uzxj1lo, wy1vqfzu, br5ovgcj, xlpjcg3s, vtsou9pz, hj3ftvzu, qfx3vhc
     *t, unhycz0e, vm4xjosb)
      endif
      endif
      call nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0kns, qf
     *x3vhct, hj3ftvzu)
      if(f7svlajr .eq. 2)then
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, rpto5qwb, hj3ftvzu, n3iasxug, vsoihn
     *1r, hbsl0gto)
      else
      rpto5qwb = -1.0d0
      endif
      do23562 iter=1,c5aesxkul 
      call intpr("iter: ",-1,iter,1)
      call intpr("posn 7: ",-1,hbsl0gto,1)
      call intpr("qfx3vhct: ",-1,qfx3vhct,1)
      call dblepr("rpto5qwb",-1,rpto5qwb,1)
      call dlgpwe0c(tlgduey8, wmat, m0ibglfx, t8hwvalr, ghz9vuba, rbne6o
     *uj, wpuarq2m, rsynp1go, n3iasxug, uaf2xgqy, kuzxj1lo, wy1vqfzu, af
     *pc0kns, br5ovgcj, dimu, hj3ftvzu, qfx3vhct, zjkrtol8, unhycz0e, vm
     *4xjosb)
      call dblepr("cqo1f: m0ibglfx",-1,m0ibglfx,wy1vqfzu*kuzxj1lo)
      call dblepr("cqo1f: wpuarq2m",-1,wpuarq2m,dimu*kuzxj1lo)
      call dblepr("cqo1f: ghz9vuba",-1,ghz9vuba,kuzxj1lo*wy1vqfzu)
      call dblepr("cqo1f: rbne6ouj",-1,rbne6ouj,kuzxj1lo*wy1vqfzu)
      do23564 gp1jxzuh=1,xlpjcg3s 
      do23566 ayfnwr1v=1,br5ovgcj 
      vc6hatuj(ayfnwr1v,gp1jxzuh) = w8znmyce(ayfnwr1v,gp1jxzuh)
23566 continue
23567 continue
23564 continue
23565 continue
      call intpr("posn 3: ",-1,hbsl0gto,1)
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23570 gp1jxzuh=1,xlpjcg3s 
      do23572 ayfnwr1v=1,kuzxj1lo 
      vc6hatuj(2*ayfnwr1v-1,gp1jxzuh) = wpuarq2m(2*hj3ftvzu-1,ayfnwr1v) 
     ** vc6hatuj(2*ayfnwr1v-1,gp1jxzuh)
      vc6hatuj(2*ayfnwr1v ,gp1jxzuh) = wpuarq2m(2*hj3ftvzu ,ayfnwr1v) * 
     *vc6hatuj(2*ayfnwr1v ,gp1jxzuh)
23572 continue
23573 continue
23570 continue
23571 continue
      else
      do23574 gp1jxzuh=1,xlpjcg3s 
      do23576 ayfnwr1v=1,kuzxj1lo 
      vc6hatuj(ayfnwr1v,gp1jxzuh) = wpuarq2m(hj3ftvzu,ayfnwr1v) * vc6hat
     *uj(ayfnwr1v,gp1jxzuh)
23576 continue
23577 continue
23574 continue
23575 continue
      endif
      call intpr("posn 4: ",-1,hbsl0gto,1)
      do23578 gp1jxzuh=1,xlpjcg3s 
      ges1xpkr(gp1jxzuh) = gp1jxzuh
23578 continue
23579 continue
      call dblepr("cqo1f: vc6hatuj",-1,vc6hatuj,br5ovgcj*xlpjcg3s)
      call intpr("iter: ",-1,iter,1)
      pvofyg8z = 1.0d-7
      call vqrdca(vc6hatuj,br5ovgcj,br5ovgcj,xlpjcg3s,fasrkub3,ges1xpkr,
     *twk,qemj9asg,pvofyg8z)
      call intpr("ges1xpkr: ",-1,ges1xpkr,xlpjcg3s)
      if(qemj9asg .ne. xlpjcg3s)then
      zjkrtol8 = 2
      return
      endif
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23584 ayfnwr1v=1,kuzxj1lo 
      twk(2*ayfnwr1v-1,1) = wpuarq2m(2*hj3ftvzu-1,ayfnwr1v) * ghz9vuba(a
     *yfnwr1v,2*hj3ftvzu-1)
      twk(2*ayfnwr1v ,1) = wpuarq2m(2*hj3ftvzu ,ayfnwr1v) * ghz9vuba(ayf
     *nwr1v,2*hj3ftvzu )
23584 continue
23585 continue
      else
      do23586 ayfnwr1v=1,kuzxj1lo
      twk(ayfnwr1v,1) = wpuarq2m(hj3ftvzu,ayfnwr1v) * ghz9vuba(ayfnwr1v,
     *hj3ftvzu)
23586 continue
23587 continue
      endif
      call intpr("posn 5: ",-1,hbsl0gto,1)
      job = 101
      call intpr("posn 6: ",-1,hbsl0gto,1)
      call vdqrsl(vc6hatuj,br5ovgcj,br5ovgcj,qemj9asg,fasrkub3, twk(1,1)
     *, uylxqtc7, twk(1,2), beta(1+(hj3ftvzu-1)*xlpjcg3s), uylxqtc7,twk(
     *1,3),job,info)
      call dblepr("beta(1+(hj3ftvzu-1)*xlpjcg3s)",-1,beta(1+(hj3ftvzu-1)
     **xlpjcg3s),xlpjcg3s)
      if(zaupqv9b .gt. 1)then
      endif
      do23590 gp1jxzuh=1,xlpjcg3s 
      twk(gp1jxzuh,1) = beta((hj3ftvzu-1)*xlpjcg3s + gp1jxzuh)
23590 continue
23591 continue
      do23592 gp1jxzuh=1,xlpjcg3s 
      beta((hj3ftvzu-1)*xlpjcg3s + ges1xpkr(gp1jxzuh)) = twk(gp1jxzuh,1)
23592 continue
23593 continue
      call intpr("posn 7: ",-1,hbsl0gto,1)
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      do23596 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = twk(2*ayfnwr1v-1,3) / wpuarq2m(2
     **hj3ftvzu-1,ayfnwr1v)
      m0ibglfx(2*hj3ftvzu ,ayfnwr1v) = twk(2*ayfnwr1v ,3) / wpuarq2m(2*h
     *j3ftvzu ,ayfnwr1v)
23596 continue
23597 continue
      if(unhycz0e .eq. 1)then
      do23600 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) = m0ibglfx(2*hj3ftvzu-1,ayfnwr1v) 
     *+ vm4xjosb(ayfnwr1v)
23600 continue
23601 continue
      endif
      else
      do23602 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = twk(ayfnwr1v,3) / wpuarq2m(hj3ftvzu,
     *ayfnwr1v)
23602 continue
23603 continue
      if(unhycz0e .eq. 1)then
      do23606 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(hj3ftvzu,ayfnwr1v) = m0ibglfx(hj3ftvzu,ayfnwr1v) + vm4xjo
     *sb(ayfnwr1v)
23606 continue
23607 continue
      endif
      endif
      call intpr("posn 8: ",-1,hbsl0gto,1)
      call nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0kns, qf
     *x3vhct, hj3ftvzu)
      call intpr("posn 8b: ",-1,hbsl0gto,1)
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, tlq9wpes,hj3ftvzu,n3iasxug,vsoihn1r,
     *hbsl0gto)
      call intpr("posn 8c: ",-1,hbsl0gto,1)
      wiptsjx8 = dabs(tlq9wpes - rpto5qwb) / (1.0d0 + dabs(tlq9wpes))
      call intpr("cqo1f: iter -------------",-1,iter,1)
      call dblepr("cqo1f: wiptsjx8",-1,wiptsjx8,1)
      if(wiptsjx8 .lt. bh2vgiay)then
      zjkrtol8 = 0
      xui7hqwl(8)=iter
      call intpr("cqo1f xui7hqwl(8): ",-1,xui7hqwl(8),1)
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, tlq9wpes,hj3ftvzu,n3iasxug,vsoihn1r,
     * izero0)
      endif
      aqg1vdmo = aqg1vdmo + tlq9wpes
      goto 1011
      else
      rpto5qwb = tlq9wpes
      endif
      call intpr("posn 9: ",-1,hbsl0gto,1)
23562 continue
23563 continue
      call intpr("cqo1f; unsuccessful convergence: ",-1,hbsl0gto,1)
      if(f7svlajr .eq. 1)then
      f7svlajr = 0
      xui7hqwl(9) = 1
      goto 653
      endif
      zjkrtol8 = 3
1011  hmayv1xt = 1.0d0
23554 continue
23555 continue
      call intpr("exiting cqo1f hbsl0gto ============================ : 
     *",-1,hbsl0gto,1)
      tlq9wpes = aqg1vdmo
      return
      end
      subroutine vcao6f(lncwkfq7, tlgduey8, wmat, m0ibglfx, t8hwvalr, gh
     *z9vuba, rbne6ouj, wpuarq2m, vc6hatuj, fasrkub3, ges1xpkr, kuzxj1lo
     *, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, xui7hqwl, tlq9wpes
     *, beta, twk, wkmm, y7sdgtqi, psdvgce3,qfozcl5b, kiye1wjz, ezlgm2up
     *, nef, which, ub4xioar,kispwgx3,s0, zyodca3j, lxyst1eb, mbvnaor6, 
     *hjm2ktyr, jnxpuym2, hnpt1zym, fzm1ihwj, iz2nbfjc, work1, wk2, wwkm
     *m, work3, sgdub, bmb, ifys6woa, mwk, ttwk, rpyis2kc, zv2xfhei, nbz
     *jkpi3, acpios9q, itwk, jwbkl9fp)
      implicit logical (a-z)
      integer xui7hqwl(19)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, ge
     *s1xpkr(*)
      double precision lncwkfq7(kuzxj1lo,*), tlgduey8(kuzxj1lo,afpc0kns)
     *, wmat(kuzxj1lo,*), m0ibglfx(wy1vqfzu,kuzxj1lo), t8hwvalr(afpc0kns
     *,kuzxj1lo)
      double precision ghz9vuba(kuzxj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy
     *1vqfzu), wpuarq2m(dimu,kuzxj1lo)
      double precision vc6hatuj(br5ovgcj,2), fasrkub3(*), tlq9wpes, beta
     *(*), y7sdgtqi(*)
      double precision twk(br5ovgcj,3), wkmm(wy1vqfzu*(wy1vqfzu+1))
      integer hj3ftvzu, ehtjigf4, izero0, iter, xwdf5ltg, dimw, f7svlajr
     *, qfx3vhct, c5aesxkul
      integer vtsou9pz, zaupqv9b, xlpjcg3s
      integer hbsl0gto, sedf7mxb
      double precision rpto5qwb, n3iasxug, wiptsjx8, uylxqtc7, bh2vgiay,
     * uaf2xgqy, vsoihn1r, rsynp1go
      double precision aqg1vdmo, hmayv1xt
      integer psdvgce3(15), qfozcl5b, ezlgm2up(*),nef(*),which(*), jnxpu
     *ym2(*), hnpt1zym(*), fzm1ihwj(*), iz2nbfjc(*)
      integer wr0lbopv, acpios9q(*), itwk(*), jwbkl9fp(*)
      integer nbzjkpi3(*)
      double precision kiye1wjz(*)
      double precision ub4xioar(qfozcl5b,kuzxj1lo), kispwgx3(kuzxj1lo,*)
     *,s0(wy1vqfzu), zyodca3j(qfozcl5b,kuzxj1lo), lxyst1eb(qfozcl5b,kuzx
     *j1lo), mbvnaor6(kuzxj1lo,*), hjm2ktyr(qfozcl5b,*), work1(*), wk2(k
     *uzxj1lo,qfozcl5b), work3(*), sgdub(*), bmb(*), ifys6woa(*), mwk(*)
     *, rpyis2kc(*), zv2xfhei(*)
      integer qes4mujl
      integer ayfnwr1v, kij0gwer, xumj5dnk
      integer irhm4cfa, lyma1kwc
      double precision xbignn(2), lncrw8mg, ufkq9rpg, r3eoxkzp, wld4qctn
      double precision zpcqv3uj, resss
      double precision vm4xjosb(2)
      lncrw8mg=0.0d0
      ufkq9rpg=0.0d0
      r3eoxkzp=0.0d0
      wld4qctn=0.0d0
      irhm4cfa = xui7hqwl(19)
      aqg1vdmo = 0.0d0
      hbsl0gto = 1
      wkmm(1) = 1.0d0
      twk(1,1) = 1.0d0
      xwdf5ltg = xui7hqwl(1)
      xlpjcg3s = xui7hqwl(3)
      dimw = xui7hqwl(4)
      f7svlajr = xui7hqwl(5)
      qfx3vhct = xui7hqwl(6)
      c5aesxkul = xui7hqwl(7)
      xui7hqwl(9) = 0
      lyma1kwc = xui7hqwl(11)
      vtsou9pz = xui7hqwl(12)
      if((vtsou9pz .ne. 1) .or. (lyma1kwc .ne. xwdf5ltg))then
      zjkrtol8 = 4
      return
      endif
      zaupqv9b = xui7hqwl(15)
      wr0lbopv = xui7hqwl(18)
      zpcqv3uj = y7sdgtqi(3+afpc0kns+afpc0kns+2)
      n3iasxug = y7sdgtqi(1)
      uaf2xgqy = dsqrt(n3iasxug)
      if((qfx3vhct .eq. 1) .or. (qfx3vhct .eq. 4))then
      vsoihn1r = dlog(n3iasxug)
      endif
      bh2vgiay = y7sdgtqi(2)
      rsynp1go = y7sdgtqi(3)
      uylxqtc7 = 0.0d0
      izero0 = 0
      zjkrtol8 = 1
      do23618 hj3ftvzu=1,afpc0kns 
653   hmayv1xt = 1.0d0
      if(f7svlajr .eq. 0)then
      call ietam6(tlgduey8, m0ibglfx, y7sdgtqi, kuzxj1lo, wy1vqfzu, afpc
     *0kns, qfx3vhct, hj3ftvzu, wmat, wr0lbopv)
      else
      if(f7svlajr .ne. 1)then
      zjkrtol8 = 6
      return
      endif
      endif
      call nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0kns, qf
     *x3vhct, hj3ftvzu)
      if(f7svlajr .eq. 2)then
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, rpto5qwb, hj3ftvzu, n3iasxug, vsoihn
     *1r, hbsl0gto)
      else
      rpto5qwb = -1.0d0
      endif
      do23626 iter=1,c5aesxkul 
      call flncwkfq76(lncwkfq7, vc6hatuj, kuzxj1lo, br5ovgcj, xwdf5ltg, 
     *qfx3vhct)
      psdvgce3(7) = 0
      call dlgpwe0c(tlgduey8, wmat, m0ibglfx, t8hwvalr, ghz9vuba, rbne6o
     *uj, wpuarq2m, rsynp1go, n3iasxug, uaf2xgqy, kuzxj1lo, wy1vqfzu, af
     *pc0kns, br5ovgcj, dimu, hj3ftvzu, qfx3vhct, zjkrtol8, izero0, vm4x
     *josb)
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      qes4mujl = 2*hj3ftvzu-1
      else
      qes4mujl = hj3ftvzu
      endif
      do23630 kij0gwer=1,qfozcl5b 
      do23632 ayfnwr1v=1,kuzxj1lo 
      zyodca3j(kij0gwer,ayfnwr1v) = wpuarq2m(qes4mujl-1+kij0gwer,ayfnwr1
     *v)
      lxyst1eb(kij0gwer,ayfnwr1v) = m0ibglfx(qes4mujl-1+kij0gwer,ayfnwr1
     *v)
23632 continue
23633 continue
23630 continue
23631 continue
      sedf7mxb = lyma1kwc * afpc0kns
      ehtjigf4 = xwdf5ltg * (hj3ftvzu-1)
      if(iter .eq. 1)then
      lncrw8mg = kiye1wjz( ehtjigf4 + hnpt1zym(1))
      ufkq9rpg = kiye1wjz(sedf7mxb + ehtjigf4 + hnpt1zym(1))
      if(xwdf5ltg .eq. 2)then
      r3eoxkzp = kiye1wjz( ehtjigf4 + hnpt1zym(2))
      wld4qctn = kiye1wjz(sedf7mxb + ehtjigf4 + hnpt1zym(2))
      endif
      do23638 kij0gwer=1,lyma1kwc 
      do23640 ayfnwr1v=1,kuzxj1lo 
      kispwgx3(ayfnwr1v,ehtjigf4 + hnpt1zym(kij0gwer)) = 0.0d0
23640 continue
23641 continue
23638 continue
23639 continue
      else
      kiye1wjz( ehtjigf4 + hnpt1zym(1)) = lncrw8mg
      kiye1wjz(sedf7mxb + ehtjigf4 + hnpt1zym(1)) = ufkq9rpg
      if(xwdf5ltg .eq. 2)then
      kiye1wjz( ehtjigf4 + hnpt1zym(2)) = r3eoxkzp
      kiye1wjz(sedf7mxb + ehtjigf4 + hnpt1zym(2)) = wld4qctn
      endif
      endif
      call vbfa(irhm4cfa,kuzxj1lo,qfozcl5b,psdvgce3, mbvnaor6, ghz9vuba(
     *1,qes4mujl), rbne6ouj(1,qes4mujl), kiye1wjz( ehtjigf4 + hnpt1zym(1
     *)), kiye1wjz(sedf7mxb + ehtjigf4 + hnpt1zym(1)), ezlgm2up,nef,whic
     *h, ub4xioar,kispwgx3(1,ehtjigf4 + hnpt1zym(1)), lxyst1eb,s0, beta(
     *1+(hj3ftvzu-1)*xlpjcg3s), cov,zpcqv3uj, vc6hatuj,fasrkub3, ges1xpk
     *r, xbignn, zyodca3j, hjm2ktyr, jnxpuym2, hnpt1zym, fzm1ihwj, iz2nb
     *fjc, work1, wk2, wwkmm, work3, sgdub, bmb, ifys6woa, mwk, ttwk, rp
     *yis2kc(1+(hj3ftvzu-1)*(nbzjkpi3(1+xwdf5ltg)-1)), zv2xfhei, resss, 
     *nbzjkpi3, acpios9q, itwk, jwbkl9fp)
      y7sdgtqi(3+afpc0kns+afpc0kns+1) = resss
      xumj5dnk = psdvgce3(14)
      if(xumj5dnk .ne. 0)then
      call intpr("vcao6f: exiting because of an error",-1,xumj5dnk,1)
      zjkrtol8 = 8
      return
      endif
      do23646 kij0gwer=1,qfozcl5b 
      do23648 ayfnwr1v=1,kuzxj1lo 
      m0ibglfx(qes4mujl-1+kij0gwer,ayfnwr1v) = lxyst1eb(kij0gwer,ayfnwr1
     *v)
23648 continue
23649 continue
23646 continue
23647 continue
      call nipyajc1(m0ibglfx, t8hwvalr, kuzxj1lo, wy1vqfzu, afpc0kns, qf
     *x3vhct, hj3ftvzu)
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, tlq9wpes, hj3ftvzu, n3iasxug, vsoihn
     *1r, hbsl0gto)
      wiptsjx8 = dabs(tlq9wpes - rpto5qwb) / (1.0d0 + dabs(tlq9wpes))
      if(wiptsjx8 .lt. bh2vgiay)then
      zjkrtol8 = 0
      xui7hqwl(8) = iter
      if((qfx3vhct .eq. 3) .or. (qfx3vhct .eq. 5))then
      call shjlwft5(qfx3vhct, tlgduey8, wmat, t8hwvalr, kuzxj1lo, wy1vqf
     *zu, afpc0kns, dimw, m0ibglfx, tlq9wpes,hj3ftvzu,n3iasxug,vsoihn1r,
     * izero0)
      endif
      aqg1vdmo = aqg1vdmo + tlq9wpes
      goto 1011
      else
      rpto5qwb = tlq9wpes
      endif
23626 continue
23627 continue
      if(f7svlajr .eq. 1)then
      f7svlajr = 0
      xui7hqwl(9) = 1
      goto 653
      endif
      zjkrtol8 = 3
1011  hmayv1xt = 1.0d0
23618 continue
23619 continue
      tlq9wpes = aqg1vdmo
      return
      end
      subroutine dcqof(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4
     *xjosb, t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj,
     * fasrkub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu,
     * zjkrtol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm
     *, y7sdgtqi, atujnxb8, yxiwebc5, k7hulceq, p2, kpzavbj3, ydcnh9xl, 
     *ajul8wkv)
      implicit logical (a-z)
      integer xui7hqwl(19), tgiyxdw1(*), dufozmt7(*)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, ge
     *s1xpkr(*)
      integer vtsou9pz
      double precision lncwkfq7(kuzxj1lo,*), tlgduey8(kuzxj1lo,afpc0kns)
     *, kifxa0he(kuzxj1lo,*), wmat(kuzxj1lo,*), m0ibglfx(wy1vqfzu,kuzxj1
     *lo), vm4xjosb(kuzxj1lo), t8hwvalr(afpc0kns,kuzxj1lo), ghz9vuba(kuz
     *xj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy1vqfzu), wpuarq2m(dimu,kuzxj1
     *lo), w8znmyce(br5ovgcj,*)
      double precision vc6hatuj(br5ovgcj,*), fasrkub3(*), tlq9wpes, beta
     *(*), y7sdgtqi(*)
      double precision twk(wy1vqfzu,kuzxj1lo,*), wkmm(wy1vqfzu*(wy1vqfzu
     *+1))
      integer p2
      double precision atujnxb8(kuzxj1lo,p2), yxiwebc5(kuzxj1lo,*), k7hu
     *lceq(p2,*), kpzavbj3(p2,*), ydcnh9xl, ajul8wkv(*)
      integer ayfnwr1v, xvr7bonh, hpmwnav2, xwdf5ltg, idlosrw8, gp1jxzuh
     *, exrkcn5d, wr0lbopv
      double precision summ, dev0
      xwdf5ltg = xui7hqwl(1)
      idlosrw8 = xui7hqwl(5)
      vtsou9pz = xui7hqwl(12)
      exrkcn5d = xui7hqwl(13)
      wr0lbopv = xui7hqwl(18)
      do23656 hpmwnav2=1,xwdf5ltg 
      do23658 ayfnwr1v=1,kuzxj1lo 
      summ = 0.0d0
      do23660 xvr7bonh=1,p2 
      summ = summ + atujnxb8(ayfnwr1v,xvr7bonh) * k7hulceq(xvr7bonh,hpmw
     *nav2)
23660 continue
23661 continue
      yxiwebc5(ayfnwr1v,hpmwnav2) = summ
      lncwkfq7(ayfnwr1v,hpmwnav2) = summ
23658 continue
23659 continue
23656 continue
23657 continue
      if(vtsou9pz.eq.1)then
      call cqo1f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, dev0, ajul8wkv, twk, wkmm, y7sd
     *gtqi)
      else
      call cqo2f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, dev0, ajul8wkv, twk, wkmm, y7sd
     *gtqi)
      endif
      do23664 xvr7bonh=1,p2 
      do23666 ayfnwr1v=1,kuzxj1lo 
      atujnxb8(ayfnwr1v,xvr7bonh) = ydcnh9xl * atujnxb8(ayfnwr1v,xvr7bon
     *h)
23666 continue
23667 continue
23664 continue
23665 continue
      do23668 hpmwnav2=1,xwdf5ltg 
      do23670 xvr7bonh=1,p2 
      do23672 ayfnwr1v=1,kuzxj1lo 
      lncwkfq7(ayfnwr1v,hpmwnav2)=yxiwebc5(ayfnwr1v,hpmwnav2)+atujnxb8(a
     *yfnwr1v,xvr7bonh)
23672 continue
23673 continue
      xui7hqwl(5) = 2
      do23674 gp1jxzuh=1,exrkcn5d 
      beta(gp1jxzuh) = ajul8wkv(gp1jxzuh)
23674 continue
23675 continue
      if(vtsou9pz.eq.1)then
      call cqo1f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm, y7sd
     *gtqi)
      else
      call cqo2f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm, y7sd
     *gtqi)
      endif
      if(zjkrtol8 .ne. 0)then
      return
      endif
      kpzavbj3(xvr7bonh,hpmwnav2) = (tlq9wpes - dev0) / ydcnh9xl
23670 continue
23671 continue
      if(xwdf5ltg .gt. 1)then
      do23682 ayfnwr1v=1,kuzxj1lo 
      lncwkfq7(ayfnwr1v,hpmwnav2) = yxiwebc5(ayfnwr1v,hpmwnav2)
23682 continue
23683 continue
      endif
23668 continue
23669 continue
      xui7hqwl(5) = idlosrw8
      return
      end
      subroutine vdcaof(lncwkfq7, tlgduey8, wmat, m0ibglfx, t8hwvalr, gh
     *z9vuba, rbne6ouj, wpuarq2m, vc6hatuj, fasrkub3, ges1xpkr, kuzxj1lo
     *, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, xui7hqwl, tlq9wpes
     *, beta, twk, wkmm, y7sdgtqi, atujnxb8, yxiwebc5, k7hulceq, p2, kpz
     *avbj3, ajul8wkv, psdvgce3,qfozcl5b, kiye1wjz, ezlgm2up, nef, which
     *, ub4xioar,kispwgx3,s0, zyodca3j, lxyst1eb, mbvnaor6, hjm2ktyr, jn
     *xpuym2, hnpt1zym, fzm1ihwj, iz2nbfjc, work1, wk2, wwkmm, work3, sg
     *dub, bmb, ifys6woa, mwk, ttwk, rpyis2kc, zv2xfhei, nbzjkpi3, acpio
     *s9q, itwk, jwbkl9fp)
      implicit logical (a-z)
      integer xui7hqwl(19)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, ge
     *s1xpkr(*)
      integer vtsou9pz
      double precision lncwkfq7(kuzxj1lo,*), tlgduey8(kuzxj1lo,afpc0kns)
     *, wmat(kuzxj1lo,*), m0ibglfx(wy1vqfzu,kuzxj1lo), t8hwvalr(afpc0kns
     *,kuzxj1lo), ghz9vuba(kuzxj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy1vqfz
     *u), wpuarq2m(dimu,kuzxj1lo)
      double precision vc6hatuj(br5ovgcj,*), fasrkub3(*), tlq9wpes, beta
     *(*), y7sdgtqi(*)
      double precision twk(wy1vqfzu,kuzxj1lo,*)
      double precision wkmm(wy1vqfzu*(wy1vqfzu+1))
      integer p2
      double precision atujnxb8(kuzxj1lo,p2), yxiwebc5(kuzxj1lo,*), k7hu
     *lceq(p2,*), kpzavbj3(p2,*), ydcnh9xl, ajul8wkv(*)
      integer ayfnwr1v, pp, hpmwnav2, xwdf5ltg, idlosrw8, exrkcn5d, wr0l
     *bopv
      double precision summ, dev0
      integer psdvgce3(15), qfozcl5b, ezlgm2up(*),nef(*),which(*), jnxpu
     *ym2(*), hnpt1zym(*), fzm1ihwj(*), iz2nbfjc(*), nbzjkpi3(2), acpios
     *9q(*), itwk(*), jwbkl9fp(2)
      double precision kiye1wjz(*)
      double precision ub4xioar(qfozcl5b,kuzxj1lo), kispwgx3(kuzxj1lo,*)
     *,s0(wy1vqfzu), zyodca3j(qfozcl5b,kuzxj1lo)
      double precision lxyst1eb(qfozcl5b,kuzxj1lo), mbvnaor6(kuzxj1lo,*)
     *, hjm2ktyr(qfozcl5b,*), work1(*), wk2(kuzxj1lo,qfozcl5b), work3(*)
     *, sgdub(*), bmb(*), ifys6woa(*), mwk(*), rpyis2kc(*), zv2xfhei(*),
     * resss
      integer irhm4cfa
      double precision zpcqv3uj
      resss = 0.0d0
      irhm4cfa = 0
      xwdf5ltg = xui7hqwl(1)
      idlosrw8 = xui7hqwl(5)
      vtsou9pz = xui7hqwl(12)
      exrkcn5d = xui7hqwl(13)
      wr0lbopv = xui7hqwl(18)
      zpcqv3uj = y7sdgtqi(3+afpc0kns+afpc0kns+2)
      ydcnh9xl = y7sdgtqi(3+afpc0kns+afpc0kns+3)
      do23684 hpmwnav2=1,xwdf5ltg 
      do23686 ayfnwr1v=1,kuzxj1lo 
      summ = 0.0d0
      do23688 pp=1,p2 
      summ = summ + atujnxb8(ayfnwr1v,pp) * k7hulceq(pp,hpmwnav2)
23688 continue
23689 continue
      yxiwebc5(ayfnwr1v,hpmwnav2) = summ
      lncwkfq7(ayfnwr1v,hpmwnav2) = summ
23686 continue
23687 continue
23684 continue
23685 continue
      if(vtsou9pz.eq.1)then
      call vcao6f(lncwkfq7, tlgduey8, wmat, m0ibglfx, t8hwvalr, ghz9vuba
     *, rbne6ouj, wpuarq2m, vc6hatuj, fasrkub3, ges1xpkr, kuzxj1lo, wy1v
     *qfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, xui7hqwl, dev0, ajul8wkv
     *, twk, wkmm, y7sdgtqi, psdvgce3,qfozcl5b, kiye1wjz, ezlgm2up, nef,
     * which, ub4xioar,kispwgx3,s0, zyodca3j, lxyst1eb, mbvnaor6, hjm2kt
     *yr, jnxpuym2, hnpt1zym, fzm1ihwj, iz2nbfjc, work1, wk2, wwkmm, wor
     *k3, sgdub, bmb, ifys6woa, mwk, ttwk, rpyis2kc, zv2xfhei, nbzjkpi3,
     * acpios9q, itwk, jwbkl9fp)
      y7sdgtqi(3+afpc0kns+afpc0kns+1) = resss
      else
      endif
      do23692 pp=1,p2 
      do23694 ayfnwr1v=1,kuzxj1lo 
      atujnxb8(ayfnwr1v,pp) = ydcnh9xl * atujnxb8(ayfnwr1v,pp)
23694 continue
23695 continue
23692 continue
23693 continue
      do23696 hpmwnav2=1,xwdf5ltg 
      do23698 pp=1,p2 
      do23700 ayfnwr1v=1,kuzxj1lo 
      lncwkfq7(ayfnwr1v,hpmwnav2) = yxiwebc5(ayfnwr1v,hpmwnav2) + atujnx
     *b8(ayfnwr1v,pp)
23700 continue
23701 continue
      xui7hqwl(5) = 0
      if(vtsou9pz.eq.1)then
      call vcao6f(lncwkfq7, tlgduey8, wmat, m0ibglfx, t8hwvalr, ghz9vuba
     *, rbne6ouj, wpuarq2m, vc6hatuj, fasrkub3, ges1xpkr, kuzxj1lo, wy1v
     *qfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, xui7hqwl, tlq9wpes, beta
     *, twk, wkmm, y7sdgtqi, psdvgce3,qfozcl5b, kiye1wjz, ezlgm2up, nef,
     * which, ub4xioar,kispwgx3,s0, zyodca3j, lxyst1eb, mbvnaor6, hjm2kt
     *yr, jnxpuym2, hnpt1zym, fzm1ihwj, iz2nbfjc, work1, wk2, wwkmm, wor
     *k3, sgdub, bmb, ifys6woa, mwk, ttwk, rpyis2kc, zv2xfhei, nbzjkpi3,
     * acpios9q, itwk, jwbkl9fp)
      y7sdgtqi(3+afpc0kns+afpc0kns+1) = resss
      else
      endif
      if(zjkrtol8 .ne. 0)then
      return
      endif
      kpzavbj3(pp,hpmwnav2) = (tlq9wpes - dev0) / ydcnh9xl
23698 continue
23699 continue
      if(xwdf5ltg .gt. 1)then
      do23708 ayfnwr1v=1,kuzxj1lo 
      lncwkfq7(ayfnwr1v,hpmwnav2) = yxiwebc5(ayfnwr1v,hpmwnav2)
23708 continue
23709 continue
      endif
23696 continue
23697 continue
      xui7hqwl(5) = idlosrw8
      return
      end
      subroutine duqof(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4
     *xjosb, t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj,
     * fasrkub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu,
     * zjkrtol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm
     *, y7sdgtqi, yxiwebc5, kpzavbj3, ydcnh9xl, ajul8wkv)
      implicit logical (a-z)
      integer xui7hqwl(19), tgiyxdw1(*), dufozmt7(*)
      integer kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrtol8, ge
     *s1xpkr(*)
      integer vtsou9pz
      double precision lncwkfq7(kuzxj1lo,*), tlgduey8(kuzxj1lo,afpc0kns)
     *, kifxa0he(kuzxj1lo,*), wmat(kuzxj1lo,*), m0ibglfx(wy1vqfzu,kuzxj1
     *lo), vm4xjosb(kuzxj1lo), t8hwvalr(afpc0kns,kuzxj1lo), ghz9vuba(kuz
     *xj1lo,wy1vqfzu), rbne6ouj(kuzxj1lo,wy1vqfzu), wpuarq2m(dimu,kuzxj1
     *lo), w8znmyce(br5ovgcj,*)
      double precision vc6hatuj(br5ovgcj,*), fasrkub3(*), tlq9wpes, beta
     *(*), y7sdgtqi(*)
      double precision twk(wy1vqfzu,kuzxj1lo,*), wkmm(wy1vqfzu*(wy1vqfzu
     *+1))
      double precision yxiwebc5(kuzxj1lo,*), kpzavbj3(kuzxj1lo,*), ydcnh
     *9xl, ajul8wkv(*)
      integer ayfnwr1v, hpmwnav2, xwdf5ltg, idlosrw8, gp1jxzuh, exrkcn5d
      double precision dev0
      xwdf5ltg = xui7hqwl(1)
      idlosrw8 = xui7hqwl(5)
      vtsou9pz = xui7hqwl(12)
      exrkcn5d = xui7hqwl(13)
      if(vtsou9pz.eq.1)then
      call cqo1f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, dev0, ajul8wkv, twk, wkmm, y7sd
     *gtqi)
      else
      call cqo2f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, dev0, ajul8wkv, twk, wkmm, y7sd
     *gtqi)
      endif
      do23712 hpmwnav2=1,xwdf5ltg 
      do23714 ayfnwr1v=1,kuzxj1lo 
      lncwkfq7(ayfnwr1v,hpmwnav2) = yxiwebc5(ayfnwr1v,hpmwnav2) + ydcnh9
     *xl
      xui7hqwl(5) = 2
      do23716 gp1jxzuh=1,exrkcn5d 
      beta(gp1jxzuh) = ajul8wkv(gp1jxzuh)
23716 continue
23717 continue
      if(vtsou9pz.eq.1)then
      call cqo1f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm, y7sd
     *gtqi)
      else
      call cqo2f(lncwkfq7, tlgduey8, kifxa0he, wmat, m0ibglfx, vm4xjosb,
     * t8hwvalr, ghz9vuba, rbne6ouj, wpuarq2m, w8znmyce, vc6hatuj, fasrk
     *ub3, ges1xpkr, kuzxj1lo, wy1vqfzu, afpc0kns, br5ovgcj, dimu, zjkrt
     *ol8, xui7hqwl, tgiyxdw1, dufozmt7, tlq9wpes, beta, twk, wkmm, y7sd
     *gtqi)
      endif
      if(zjkrtol8 .ne. 0)then
      return
      endif
      kpzavbj3(ayfnwr1v,hpmwnav2) = (tlq9wpes - dev0) / ydcnh9xl
      lncwkfq7(ayfnwr1v,hpmwnav2) = yxiwebc5(ayfnwr1v,hpmwnav2)
23714 continue
23715 continue
23712 continue
23713 continue
      xui7hqwl(5) = idlosrw8
      return
      end
