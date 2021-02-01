c---------------------------------------------------------------(hill48)
c
c     Hill 1948 yield function and its dfferentials
c
      subroutine jancae_hill_1948 ( stress,ndyld,pryld,nreq,
     1                              se,dseds )
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension stress(6)
      dimension s(6),dseds(6),pryld(ndyld)
      dimension c(6,6),v(6)
c
c ----------------------------------------------------------------------
c
cf2py intent(in) stress,ndyld,pryld,nreq
cf2py intent(out) se,dseds
cf2py depend(ndyld) pryld
cf2py depend(6) dseds
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
c                                            ---- anisotropic parameters
      pf = pryld(1)
      pg = pryld(2)
      ph = pryld(3)
      pl = pryld(4)
      pm = pryld(5)
      pn = pryld(6)
c                                                      ---- coef. matrix
      c = 0.0d0
      c(1,1) =  pg + ph
      c(1,2) = -ph
      c(1,3) = -pg
      c(2,1) = -ph
      c(2,2) =  pf + ph
      c(2,3) = -pf
      c(3,1) = -pg
      c(3,2) = -pf
      c(3,3) =  pf + pg
      c(4,4) = 2.0d0 * pn
      c(5,5) = 2.0d0 * pm
      c(6,6) = 2.0d0 * pl
      do i = 1,6
        do j = 1,6
          c(i,j) = c(i,j) / (pg+ph)
        enddo
      enddo
c
      call jancae_mv  ( v,c,s,6,6 )
      call jancae_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      if ( phi .le. 0.0 ) phi = 0.0
      se = sqrt(phi)
c                                            ---- 1st order differential
      if ( nreq .ge. 1 ) then
        do i = 1,6
          dseds(i) = v(i) / se
        end do
      end if
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
      subroutine jancae_mv ( v,a,u,nv,nu )
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
