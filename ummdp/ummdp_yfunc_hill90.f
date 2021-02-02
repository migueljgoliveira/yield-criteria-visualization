c---------------------------------------------------------------(hill90)
c
c     Hill-1990 anisotropic yield function and its dfferentials
c
      subroutine jancae_hill90 ( stress,ndyld,pryld,nreq,
     1                           se,dseds )
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension stress(6)
      dimension s(3),dseds(3),pryld(ndyld)
c
      dimension c(3,3),v(3)
      dimension A1(3),A2(3,3),A3(3,3),A4(3,3)
      dimension dfds1(3),dfds2(3),dfds3(3),dfds4(3)
      dimension dfds_t(3)
      dimension dxds1(3),dxds2(3),dxds3(3),dxds4(3)
c
c ----------------------------------------------------------------------
c
cf2py intent(in) stress,ndyld,pryld,nreq
cf2py intent(out) se,dseds
cf2py depend(ndyld) pryld
cf2py depend(3) dseds
c
c ----------------------------------------------------------------------
c
c                                                     ---- stress tensor
      s(1) = stress(1)
      s(2) = stress(2)
      s(3) = stress(4)
c                                       ---- define a1-matrix, a2-matrix
      data a1/ 1.0d0 , 1.0d0 , 0.0d0 /,
     &     a2/ 1.0d0 ,-1.0d0 , 0.0d0 ,
     &        -1.0d0 , 1.0d0 , 0.d00 ,
     &         0.0d0 , 0.0d0 , 4.0d0 /,
     &     a3/ 1.0d0 , 0.0d0 , 0.0d0 ,
     &         0.0d0 , 1.0d0 , 0.0d0 ,
     &         0.0d0 , 0.0d0 , 2.0d0 /
c                                        ---- set anisotropic parameters
      a    = pryld(1+1)
      b    = pryld(1+2)
      tau  = pryld(1+3)
      sigb = pryld(1+4)
      am   = pryld(1+5)
c
      syini = 1.0d0
      sigbtm = (sigb/tau)**am
      alarge = 1.0d0 + sigbtm - 2.0d0*a + b
c
c                               ---- coef. matrix of material parameters
c                                ---- define a4-matrix consists of a & b
      a4 = 0.0d0
      a4(1,1) = -2.0d0*a + b
      a4(2,2) = 2.0d0*a + b
      a4(1,2) = - b
      a4(2,1) = - b
c                                                              ---- fai1
      x1 = s(1) + s(2)
      fai1 = abs(x1)**am
c                                                              ---- fai2
      call jancae_mv  ( v,a2,s,3,3 )
      call jancae_vvs ( x2,s,v,3 )
      fai2 = sigbtm * (x2)**(am/2.0d0)
c                                                              ---- fai3
      call jancae_mv  ( v,a3,s,3,3 )
      call jancae_vvs ( x3,s,v,3 )
      fai3 = (x3)**(am/2.0d0-1.0d0)
c                                                              ---- fai4
      call jancae_mv  ( v,a4,s,3,3 )
      call jancae_vvs ( x4,s,v,3 )
      fai4 = x4
c                                             ---- yield fuction : fyild
      fyild = fai1 + fai2 + fai3*fai4
c
c                                                 ---- equivalent stress
      se = (fyild/alarge)**(1.0d0/am)
c                                             ---- 1st order differential
      if ( nreq .ge. 1 ) return
c                                                          ---- dfai1/ds
      dxds1(1) = 1.0
      dxds1(2) = 1.0
      dxds1(3) = 0.0
c
      wrk = am * (abs(x1)**(am-2))*x1
      do i = 1,3
        dfds1(i) = wrk * dxds1(i)
      end do
c                                                          ---- dfai2/ds
      wrk = sigbtm * (am/2.0) * (x2)**(am/2.0-1.0)
      call jancae_mv ( dxds2,a2,s,3,3 )
      do i = 1,3
        dxds2(i) = 2.0 * dxds2(i)
        dfds2(i) = wrk * dxds2(i)
      end do
c                                                          ---- dfai3/ds
      wrk = (am/2.0-1.0) * (x3)**(am/2.0-2.0)
      call jancae_mv ( dxds3,a3,s,3,3 )
      do i = 1,3
        dxds3(i) = 2.0 * dxds3(i)
        dfds3(i) = wrk * dxds3(i)
      end do
c                                                          ---- dfai4/ds
      call jancae_mv ( dxds4,a4,s,3,3 )
      do i = 1,3
        dxds4(i) = 2.0 * dxds4(i)
        dfds4(i) = dxds4(i)
      end do
c
c        ---- 1st order differential coefficient of yield fuction result
c                                     ---- dfai/ds()= result = dfds_t(i)
      do i = 1,3
        dfds_t(i) = dfds1(i) + dfds2(i) + dfds3(i)*fai4 + fai3*dfds4(i)
      end do
c           ---- 1st order differential coefficient of equivalent stress
      wrk = (abs(fyild/alarge))**(1.0/am-1.0) / (am*alarge)
      do i = 1,3
        dseds(i) = wrk * dfds_t(i)
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate scaler product of vectors
c     s={u}T{v}
c
      subroutine jancae_vvs ( s,u,v,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v(n),u(n)
c
      s = 0.0
      do i = 1,n
        s = s + u(i)*v(i)
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate multiplication of matrix and vector
c     {v}=[a]{u}
c
      subroutine jancae_mv (v,a,u,nv,nu)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v(nv),a(nv,nu),u(nu)
c
      v = 0.0d0
      do i = 1,nv
        do j = 1,nu
          v(i) = v(i) + a(i,j)*u(j)
        enddo
      enddo
c
      return
      end
c
c
c