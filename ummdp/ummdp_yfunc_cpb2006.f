c--------------------------------------------------------------(cpb2006)
c
c     CPB 2006 yield function and its dfferentials
c
      subroutine jancae_cpb2006 ( stress,ndyld,pryld,nreq,
     1                            se,dseds )
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension stress(6)
      dimension s(6),dseds(6),pryld(ndyld)
      dimension sigma(6),psigma(3)
      dimension phi(3),psi(3),omega(3)
      dimension c(6,6),ct(6,6)
      dimension DFDH(3),DFDpsigma(3)
      dimension DpsigmaDH(3,3),DHdsigma(3,6),DsigmaDs(6,6)
      dimension DFDs(6)
      dimension dummat(3,6)
      parameter ( pi=3.141592653589793d0 )
      parameter ( eps=1.0d-5 )
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
c
c                                        ---- set anisotropic parameters
c
      c = 0.0d0
c
      c(1,1) = pryld(1)            ! C11 = 1.0
      c(1,2) = pryld(2)            ! C12 = 0.0
      c(1,3) = pryld(3)            ! C13 = 0.0
      c(2,1) = pryld(4)            ! C21 = 0.0
      c(2,2) = pryld(5)            ! C22 = 1.0
      c(2,3) = pryld(6)            ! C23 = 0.0
      c(3,1) = pryld(7)            ! C31 = 0.0
      c(3,2) = pryld(8)            ! C32 = 0.0
      c(3,3) = pryld(9)            ! C33 = 1.0
      c(4,4) = pryld(10)           ! C44 = 1.0    ! tau_xy
      c(5,5) = pryld(11)           ! C55 = 1.0    ! tau_xz
      c(6,6) = pryld(12)           ! C66 = 1.0    ! tau_yz
      a      = pryld(13)           ! a = 2.0
      ck     = pryld(14)           ! k = 0.0
c
      ai = 1.0d0 / a
c
c                                         ---- Calculate phi, psi, omega
      phi(1) = (2.0d0*c(1,1) - c(1,2) - c(1,3)) / 3.0d0
      phi(2) = (2.0d0*c(1,2) - c(2,2) - c(2,3)) / 3.0d0
      phi(3) = (2.0d0*c(1,3) - c(2,3) - c(3,3)) / 3.0d0
c
      psi(1) = (-c(1,1) + 2.0d0*c(1,2) - c(1,3)) / 3.0d0
      psi(2) = (-c(1,2) + 2.0d0*c(2,2) - c(2,3)) / 3.0d0
      psi(3) = (-c(1,3) + 2.0d0*c(2,3) - c(3,3)) / 3.0d0
c
      omega(1) = (c(1,1) + c(1,2) - 2.0d0*c(1,3)) / 3.0d0
      omega(2) = (c(1,2) + c(2,2) - 2.0d0*c(2,3)) / 3.0d0
      omega(3) = (c(1,3) + c(2,3) - 2.0d0*c(3,3)) / 3.0d0
c
c      ---- Calculate 4th order orthotropic tensor "L" ( named ct here )
      ct = 0.0d0
c
      ct(1,1) =    phi(1)
      ct(1,2) =    psi(1)
      ct(1,3) = -omega(1)
      ct(2,1) =    phi(2)
      ct(2,2) =    psi(2)
      ct(2,3) = -omega(2)
      ct(3,1) =    phi(3)
      ct(3,2) =    psi(3)
      ct(3,3) = -omega(3)
      ct(4,4) =    c(4,4)
      ct(5,5) =    c(5,5)
      ct(6,6) =    c(6,6)
c
c                               ---- Calculate linear transformed stress
      call jancae_mv ( sigma,ct,s,6,6 )
c
c            ---- Calculate principal values of transformed stress sigma
c                                                     by Cardan's method
c
c                                       ---- 1st, 2nd and 3rd invariants
      H1 = (sigma(1) + sigma(2) + sigma(3)) / 3.0d0
      H2 = (sigma(5)**2.0d0 + sigma(6)**2.0d0 + sigma(4)**2.0d0 -
     1      sigma(2)*sigma(3) - sigma(3)*sigma(1) - 
     2      sigma(1)*sigma(2)) / 3.0d0
      H3 = (2.0d0*sigma(5)*sigma(6)*sigma(4) + 
     1      sigma(1)*sigma(2)*sigma(3) - sigma(1)*sigma(6)**2.0d0 -
     2      sigma(2)*sigma(5)**2.0d0 -
     3      sigma(3)*sigma(4)**2.0d0) / 2.0d0 ! sigma(5) <-> sigma(6)
c
      p = H1**2.0d0 + H2
      q = (2.0d0*H1**3.0d0 + 3.0d0*H1*H2 + 2.0d0*H3) / 2.0d0
      if ( abs(p) .ge. 1.0d-16 ) then
        theta = q / (p**1.5d0)
        if ( theta >  1.0d0 ) theta =  1.0d0
        if ( theta < -1.0d0 ) theta = -1.0d0
        theta = acos(theta)
      else
        theta = 0.0
      end if
c
c                               ---- calculate principal values of sigma
      psigma(1) = 2.0d0*sqrt(p)*cos(theta/3.0d0) + H1
      psigma(2) = 2.0d0*sqrt(p)*cos((theta+4.0d0*pi)/3.0d0) + H1
      psigma(3) = 2.0d0*sqrt(p)*cos((theta+2.0d0*pi)/3.0d0) + H1
c
c                                                 ---- equivalent stress
c
c                                          ---- calculate yield function
      F = (abs(psigma(1)) - ck*psigma(1))**a + 
     1    (abs(psigma(2)) - ck*psigma(2))**a + 
     2    (abs(psigma(3)) - ck*psigma(3))**a
c                                           ---- denominator coefficient
      D = (abs(phi(1)) - ck*phi(1))**a + (abs(phi(2)) - ck*phi(2))**a +
     1    (abs(phi(3)) - ck*phi(3))**a
c
      se = (F/D) ** ai
c
c                                            ---- 1st order differential
      if ( nreq .ge. 1 ) then
c                                              ---- D(se)/D(F) -> Scalar
        DseDF = (1.0d0/D)**ai * ai * F**(ai-1.0d0)
c                                    ---- D(F)/D(psigma) -> 1 x 3 Vector
        do i = 1,3
          DFDpsigma(i) = a * (psigma(i)/abs(psigma(i))-ck) *
     1                   (abs(psigma(i))-ck*psigma(i))**(a-1.0d0)
        end do
c
c                                 ---- D(F)/D(H) by using D(psigma)/D(H)
c                                         D(F)/D(H)      -> 1 x 3 Vector
c                                         D(psigma)/D(H) -> 3 x 3 Matrix
        DFDH = 0.0d0
        DpsigmaDH = 0.0d0
c
        if ( abs(psigma(2)-psigma(3)) / se > eps .and.
     1       abs(psigma(2)-psigma(1)) / se > eps ) then
c                                                 ---- not Singular case
          do i = 1,3
            denom = psigma(i)**2.0d0 - 2.0d0*H1*psigma(i) - H2
            DpsigmaDH(i,1) = psigma(i)**2.0d0/denom
            DpsigmaDH(i,2) = psigma(i)/denom
            DpsigmaDH(i,3) = 2.0d0/3.0d0/denom
          end do
c
          do i = 1,3
            do j = 1,3
              DFDH(i) = DFDH(i) + DFDpsigma(j)*DpsigmaDH(j,i)
            end do
          end do
        else
c                                                     ---- singular case
          if ( abs(psigma(2)-psigma(3)) / se .le. eps) then
c                                           ---- Case1 S2=S3 ( theta=0 )
            denom = psigma(1)**2.0d0 - 2.0d0*H1*psigma(1) - H2
            DpsigmaDH(1,1) = psigma(1)**2.0d0/denom
            DpsigmaDH(1,2) = psigma(1) / denom
            DpsigmaDH(1,3) = 2.0d0/3.0d0/denom
c
            DFDH(1) = DpsigmaDH(1,1)*(DFDpsigma(1)-DFDpsigma(2)) +
     1                3.0d0*DFDpsigma(2)
            DFDH(2) = DpsigmaDH(1,2) * (DFDpsigma(1)-DFDpsigma(2))
            DFDH(3) = DpsigmaDH(1,3) * (DFDpsigma(1)-DFDpsigma(2))
          elseif ( abs(psigma(2)-psigma(1)) / se .le. eps ) then
c
c                                          ---- Case2 S2=S1 ( theta=pi )
            denom = psigma(3)**2.0d0 - 2.0d0*H1*psigma(3) - H2
            DpsigmaDH(3,1) = psigma(3)**2.0d0/denom
            DpsigmaDH(3,2) = psigma(3) / denom
            DpsigmaDH(3,3) = 2.0d0/3.0d0/denom
c
            DFDH(1) = DpsigmaDH(3,1)*(DFDpsigma(3)-DFDpsigma(2)) +
     1                3.0d0*DFDpsigma(2)
            DFDH(2) = DpsigmaDH(3,2) * (DFDpsigma(3)-DFDpsigma(2))
            DFDH(3) = DpsigmaDH(3,3) * (DFDpsigma(3)-DFDpsigma(2))
          endif
        endif
c
c                                     ---- D(H)/D(sigma) -> 3 x 6 Matrix
        DHDsigma = 0.0d0
c
        DHDsigma(1,1) = 1.0d0 / 3.0d0
        DHDsigma(1,2) = 1.0d0 / 3.0d0
        DHDsigma(1,3) = 1.0d0 / 3.0d0
c
        DHDsigma(2,1) = -1.0d0/3.0d0*(sigma(2)+sigma(3))
        DHDsigma(2,2) = -1.0d0/3.0d0*(sigma(3)+sigma(1))
        DHDsigma(2,3) = -1.0d0/3.0d0*(sigma(1)+sigma(2))
        DHDsigma(2,4) = 2.0d0/3.0d0*sigma(4)
        DHDsigma(2,5) = 2.0d0/3.0d0*sigma(5)
        DHDsigma(2,6) = 2.0d0/3.0d0*sigma(6)
c
c                                            !!-sigma(5)**2.0d0)
        DHDsigma(3,1) = 0.5d0 * (sigma(2)*sigma(3)-sigma(6)**2.0d0)
c                                            !!-sigma(6)**2.0d0)
        DHDsigma(3,2) = 0.5d0 * (sigma(3)*sigma(1)-sigma(5)**2.0d0)
        DHDsigma(3,3) = 0.5d0 * (sigma(1)*sigma(2)-sigma(4)**2.0d0)
        DHDsigma(3,4) = sigma(5)*sigma(6) - sigma(3)*sigma(4)
        DHDsigma(3,6) = sigma(6)*sigma(4) - sigma(1)*sigma(5) !!...(3,5)=
        DHDsigma(3,5) = sigma(4)*sigma(5) - sigma(2)*sigma(6) !!...(3,6)=
c
c                                     ---- D(sigma)/D(s) -> 6 x 6 Matrix
        do i = 1,6
          do j = 1,6
            DsigmaDs(i,j) = ct(i,j)
          end do
        end do
c
c                                        ---- D(se)/D(s) -> 1 x 3 Vector
        DFDs = 0.0d0
        dummat = 0.0d0
c
        do i = 1,6
          do j = 1,3
            do k = 1,6
              dummat(j,i) = dummat(j,i) + DHDsigma(j,k)*DsigmaDs(k,i)
            end do
            DFDs(i) = DFDs(i) + DFDH(j)*dummat(j,i)
          end do
        end do
      end if
c
      do i = 1,6
        dseds(i) = DseDF * DFDs(i)
      end do
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