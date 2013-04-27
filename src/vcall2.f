C Output from Public domain Ratfor, version 1.01
      subroutine vcall2(onemor,w,y,eta,beta,u)
      logical onemor
      double precision w(*), y(*), eta(*), beta(*), u(*)
      onemor = .true.
      w(1) = 1.0d0
      y(1) = 1.0d0
      eta(1) = 1.0d0
      beta(1) = 1.0d0
      u(1) = 1.0d0
      return
      end
      subroutine vcall1(onemor,y,eta,beta,u,xbig,cpxbig)
      logical onemor, cpxbig
      double precision y(*), eta(*), beta(*), u(*), xbig(*)
      onemor = .true.
      y(1) = 1.0d0
      eta(1) = 1.0d0
      beta(1) = 1.0d0
      u(1) = 1.0d0
      xbig(1) = 1.0d0
      cpxbig = .true.
      return
      end
