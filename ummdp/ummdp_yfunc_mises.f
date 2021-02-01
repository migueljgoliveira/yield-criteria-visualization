c----------------------------------------------------------------(mises)
c
c     von Mises yield function
c
      subroutine jancae_mises ( stress,ndyld,pryld,nreq,
     1                          se,dseds )
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension stress(6),pryld(ndyld)
      dimension s(6),dseds(6)
      dimension c(6,6),v(6)
c
c ----------------------------------------------------------------------
c
cf2py intent(in) stress,ndyld,pryld,nreq
cf2py intent(out) se,dseds
cf2py depend(ndyld) pryld
cf2py depend(6) stress,dseds
c
c ----------------------------------------------------------------------
c
c                                                     ---- stress tensor
      s(1) = stress(1)
      s(2) = stress(2)
      s(3) = stress(3)
      s(4) = stress(4)
      s(5) = stress(5)
      s(6) = stress(6)
c                                                      ---- coef. matrix
      c = 0.0d0
      do i = 1,3
        do j = 1,3
          c(i,j) = -0.5d0
        enddo
        c(i,i) = 1.0d0
      enddo
      do i = 4,6
        c(i,i) = 3.0d0
      enddo
c
      call jancae_mv  ( v,c,s,6,6 )
      call jancae_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      se = sqrt(phi)
c
      return
      end
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