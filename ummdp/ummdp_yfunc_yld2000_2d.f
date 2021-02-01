c-----------------------------------------------------------(yld2000_2d)
c
c     Yld2000-2d yield function and its dfferentials
c
      subroutine jancae_yld2000_2d ( stress,ndyld,pryld,nreq,
     1                               se,dseds )
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension stress(6)
      dimension s(3),dseds(3),pryld(ndyld)
      dimension a(8),am(2,3,3),xx(2,2),yy(2,3),phi(2)
      dimension dsedphi(2),dphidx(2,2)
      dimension dxdy(2,2,3),dyds(2,3,3)
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
c                                            ---- anisotropic parameters
      do i = 1,8
        a(i) = pryld(i)
      enddo
      em = pryld(9)
c                                         ---- set linear transf. matrix
      call jancae_yld2000_2d_am ( a,am )
c                                                 ---- equivalent stress
      call jancae_yld2000_2d_xyphi ( s,em,am,xx,yy,phi )
      q = phi(1) + phi(2)
      if ( q .le. 0.0 ) q = 0.
c
      se = (0.5d0*q)**(1.0d0/em)
c                                            ---- 1st order differential
      if ( nreq .ge. 1 ) then
        call jancae_yld2000_2d_ds1 ( em,am,xx,yy,phi,
     1                               dsedphi,dphidx,
     2                               dxdy,dyds,se )
        dseds = 0.0d0
        do nd = 1,2
          do m = 1,2
            do k = 1,3
              do i = 1,3
                dseds(i) = dseds(i) + dsedphi(nd)*dphidx(nd,m)*
     1                                dxdy(nd,m,k)*dyds(nd,k,i)
              end do
            end do
          end do
        end do
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------(yld2000_2d)
c     set barlat-yld2k linear transformation matrix am
c
      subroutine jancae_yld2000_2d_am ( a,am )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(8),am(2,3,3)
c
c                                      ---- linear transformation matrix
      am(1,1,1) =  2.0d0*a(1)
      am(1,1,2) = -1.0d0*a(1)
      am(1,1,3) =  0.0
      am(1,2,1) = -1.0d0*a(2)
      am(1,2,2) =  2.0d0*a(2)
      am(1,2,3) =  0.0
      am(1,3,1) =  0.0
      am(1,3,2) =  0.0
      am(1,3,3) =  3.0d0*a(7)
c
      am(2,1,1) = -2.0d0*a(3) + 2.0d0*a(4) + 8.0d0*a(5) - 2.0d0*a(6)
      am(2,1,2) =        a(3) - 4.0d0*a(4) - 4.0d0*a(5) + 4.0d0*a(6)
      am(2,1,3) =  0.0
      am(2,2,1) =  4.0d0*a(3) - 4.0d0*a(4) - 4.0d0*a(5) +       a(6)
      am(2,2,2) = -2.0d0*a(3) + 8.0d0*a(4) + 2.0d0*a(5) - 2.0d0*a(6)
      am(2,2,3) =  0.0
      am(2,3,1) =  0.0
      am(2,3,2) =  0.0
      am(2,3,3) =  9.0d0*a(8)
c
      do i = 1,3
        do j = 1,3
          am(1,i,j) = am(1,i,j) / 3.0d0
          am(2,i,j) = am(2,i,j) / 9.0d0
        end do
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------(yld2000_2d)
c     calc. barlat-yld2k function xx,yy,phi
c
      subroutine jancae_yld2000_2d_xyphi ( s,em,am,xx,yy,phi )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),am(2,3,3),xx(2,2),yy(2,3),phi(2)
      dimension p(2)
c
      p(1) =  1.0d0
      p(2) = -1.0d0
c                                                      ---- {yy}=[am]{s}
      yy = 0.0d0
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            yy(nd,i) = yy(nd,i) + am(nd,i,j)*s(j)
          end do
        end do
      end do
c                                      ---- {xx}=principle value of {yy}
      do nd = 1,2
        a = (yy(nd,1)-yy(nd,2))**2.d0 + 4.0d0*yy(nd,3)**2.d0
        a = sqrt(a)
        do i = 1,2
          xx(nd,i) = 0.5d0*(yy(nd,1)+yy(nd,2)+p(i)*a)
        end do
      end do
c                                                 ---- phi(1) and phi(2)
      nd = 1
      phi(nd) = abs(xx(nd,1)-xx(nd,2))**em
      nd = 2
      phi(nd) = abs(2.0d0*xx(nd,2)+xx(nd,1))**em +
     1          abs(2.0d0*xx(nd,1)+xx(nd,2))**em
c
      return
      end
c
c
c
c-----------------------------------------------------------(yld2000_2d)
c     set 1st order differential of parameters
c
      subroutine jancae_yld2000_2d_ds1 ( em,am,xx,yy,phi,
     1                                   dsedphi,dphidx,
     2                                   dxdy,dyds,se )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension am(2,3,3),xx(2,2),yy(2,3),phi(2),
     1          dsedphi(2),dphidx(2,2),
     2          dxdy(2,2,3),dyds(2,3,3)
c
      dimension p(2)
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                          ---- dse/dphi
      q = phi(1) + phi(2)
      if ( q .le. 0.0 ) q = 0.0
      do nd = 1,2
        dsedphi(nd) = (0.5d0**emi)*emi*q**(emi-1.0d0)
      end do
c                                                           ---- dphi/dx
      nd = 1
      a0 = xx(nd,1) - xx(nd,2)
      b0 = abs(a0)
      sgn0 = 0
      if ( b0 .ge. eps*se ) sgn0 = a0/b0
      dphidx(nd,1) =  em*b0**(em-1.0d0) * sgn0
      dphidx(nd,2) = -em*b0**(em-1.0d0) * sgn0
c
      nd = 2
      a1 = 2.0d0*xx(nd,1) +       xx(nd,2)
      a2 =       xx(nd,1) + 2.0d0*xx(nd,2)
      b1 = abs(a1)
      b2 = abs(a2)
      sgn1 = 0.0
      sgn2 = 0.0
      if ( b1 .ge. eps*se ) sgn1 = a1 / b1
      if ( b2 .ge. eps*se ) sgn2 = a2 / b2
      dphidx(nd,1) = em*(2.0d0*b1**(em-1.0d0)*sgn1 + 
     1                         b2**(em-1.0d0)*sgn2)
      dphidx(nd,2) = em*(      b1**(em-1.0d0)*sgn1 + 
     1                   2.0d0*b2**(em-1.0d0)*sgn2)
c
      do nd = 1,2
        a = (yy(nd,1)-yy(nd,2))*(yy(nd,1)-yy(nd,2)) + 
     1      4.0d0*yy(nd,3)*yy(nd,3)
        a = sqrt(a)
        if ( a .gt. eps*se ) then
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0*(1.0d0+p(j)*(yy(nd,1)-yy(nd,2))/a)
            dxdy(nd,j,2) = 0.5d0*(1.0d0-p(j)*(yy(nd,1)-yy(nd,2))/a)
            dxdy(nd,j,3) = 2.0d0*       p(j)* yy(nd,3)         /a
          end do
        else
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0*(1.0d0+0.0)
            dxdy(nd,j,2) = 0.5d0*(1.0d0-0.0)
            dxdy(nd,j,3) = 2.0d0*       0.0
          end do
        end if
      end do
c
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            dyds(nd,i,j) = am(nd,i,j)
          end do
        end do
      end do
c
      return
      end
c
c
c
