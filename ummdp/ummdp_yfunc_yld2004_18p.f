c----------------------------------------------------------(yld2004-18p)
c
c     Yld2004-18p yield function and its differentials
c
      subroutine jancae_yld2004_18p ( stress,ndyld,pryld,nreq,
     1                                se,dseds )
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension stress(6)
      dimension s(6),dseds(6),pryld(ndyld)
      dimension sp1(6),sp2(6),cp1(6,6),cp2(6,6),cl(6,6),ctp1(6,6)
      dimension ctp2(6,6),psp1(3),psp2(3),hp1(3),hp2(3)
      dimension dfadpsp1(3),dfadpsp2(3),dpsdhp1(3,3),dpsdhp2(3,3)
      dimension dfadhp1(3),dfadhp2(3)
      dimension dhdsp1(3,6),dhdsp2(3,6),dsdsp1(6,6), dsdsp2(6,6)
      dimension dfads(6)
      dimension delta(3,3)
      dimension xx1(3,6),xx2(3,6)
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
      pi = acos(-1.0d0)
      eps2 = 1.0d-15
      eps3 = 1.0d-8
      del = 1.0d-4
c                                                   ---- Kronecker Delta
      delta = 0.0d0
      do i = 1,3
        delta(i,i) = 1.0d0
      end do
c                                        ---- set anisotropic parameters
      cp1 = 0.0d0
      cp2 = 0.0d0
      cp1(1,2) = -pryld(1)
      cp1(1,3) = -pryld(2)
      cp1(2,1) = -pryld(3)
      cp1(2,3) = -pryld(4)
      cp1(3,1) = -pryld(5)
      cp1(3,2) = -pryld(6)
      cp1(4,4) =  pryld(9)  ! for tau_xy (c'66 in original paper)
      cp1(5,5) =  pryld(8)  ! for tau_xz (c'55 in original paper)
      cp1(6,6) =  pryld(7)  ! for tau_zx (c'44 in original paper)
      cp2(1,2) = -pryld(10)
      cp2(1,3) = -pryld(11)
      cp2(2,1) = -pryld(12)
      cp2(2,3) = -pryld(13)
      cp2(3,1) = -pryld(14)
      cp2(3,2) = -pryld(15)
      cp2(4,4) =  pryld(18) ! for tau_xy (c"66 in original paper)
      cp2(5,5) =  pryld(17) ! for tau_xz (c"55 in original paper)
      cp2(6,6) =  pryld(16) ! for tau_zx (c"44 in original paper)
      am       =  pryld(19)
      ami = 1.0d0 / am
c
      dc = 4.0d0
c
c             ---- set matrix for transforming Cauchy stress to deviator
      cl = 0.0d0
      do i = 1,3
        do j = 1,3
          if ( i == j ) then
            cl(i,j) = 2.0d0
          else
            cl(i,j) = -1.0d0
          end if
        end do
      end do
      do i = 4,6
        cl(i,i) = 3.0d0
      end do
      do i = 1,6
        do j = 1,6
          cl(i,j) = cl(i,j) / 3.0d0
        end do
      end do
c
c                  ---- matrix for transforming Cauchy stress to sp1,sp2
      call jancae_mm ( ctp1,cp1,cl,6,6,6 )
      call jancae_mm ( ctp2,cp2,cl,6,6,6 )
c                                  ---- calculation of equivalent stress
      call jancae_yld2004_18p_yf ( ctp1,ctp2,s,am,ami,dc,pi,
     &                             sp1,sp2,psp1,psp2,hp1,hp2,
     &                             cetpq1,cetpq2,fai,se )
c
c                                            ---- 1st order differential
      if ( nreq .ge. 1 ) then
c                                                     ---- d(fai)/d(psp)
        do i = 1,3
          dfadpsp1(i) = 0.0d0
          dfadpsp2(i) = 0.0d0
          do j = 1,3
            dfadpsp1(i) = dfadpsp1(i) + (psp1(i)-psp2(j)) *
     &                    abs(psp1(i)-psp2(j))**(am-2.0d0)
            dfadpsp2(i) = dfadpsp2(i) + (psp1(j)-psp2(i)) *
     &                    abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
          dfadpsp1(i) = dfadpsp1(i) * am
          dfadpsp2(i) = dfadpsp2(i) * (-am)
        end do
c                                         ---- d(psp)/d(hp)&d(fai)/d(hp)
        dpsdhp1 = 0.0d0
        dpsdhp2 = 0.0d0
        dfadhp1 = 0.0d0
        dfadhp2 = 0.0d0
c                                                  ---- theta'<>0 & <>pi
        if ( abs(cetpq1-1.0d0) >= eps2 .and. 
     &       abs(cetpq1+1.0d0) >= eps2 ) then
          do i = 1,3
            call jancae_yld2004_18p_dpsdhp ( i,psp1,hp1,dpsdhp1 )
          end do
c                                                          ---- theta'=0
        else if ( abs(cetpq1-1.0d0) < eps2 ) then
          i = 1
          call jancae_yld2004_18p_dpsdhp ( i,psp1,hp1,dpsdhp1 )
          do i = 2,3
            do j = 1,3
              dpsdhp1(i,j) = -0.5d0 * (dpsdhp1(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                         ---- theta'=pi
        else
          i = 3
          call jancae_yld2004_18p_dpsdhp ( i,psp1,hp1,dpsdhp1 )
          do i = 1,2
            do j = 1,3
              dpsdhp1(i,j) = -0.5d0 * (dpsdhp1(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c                                                 ---- theta''<>0 & <>pi
        if ( abs(cetpq2-1.0d0) >= eps2 .and. 
     &       abs(cetpq2+1.0d0) >= eps2 ) then
          do i = 1,3
            call jancae_yld2004_18p_dpsdhp ( i,psp2,hp2,dpsdhp2 ) 
          end do
c                                                         ---- theta''=0
        else if ( abs(cetpq2-1.0d0) < eps2 ) then
          i = 1
          call jancae_yld2004_18p_dpsdhp ( i,psp2,hp2,dpsdhp2 )
          do i = 2,3
            do j = 1,3
              dpsdhp2(i,j) = -0.5d0 * (dpsdhp2(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                        ---- theta''=pi
        else
          i = 3
          call jancae_yld2004_18p_dpsdhp ( i,psp2,hp2,dpsdhp2 )
          do i = 1,2
            do j = 1,3
              dpsdhp2(i,j) = -0.5d0 * (dpsdhp2(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c
        do i = 1,3
          do j = 1,3
            dfadhp1(i) = dfadhp1(i) + dfadpsp1(j)*dpsdhp1(j,i)
            dfadhp2(i) = dfadhp2(i) + dfadpsp2(j)*dpsdhp2(j,i)
          end do
        end do
c                                                       ---- d(hp)/d(sp)
        dhdsp1 = 0.0d0
        dhdsp2 = 0.0d0
        do i = 1,3
          j = mod(i,  3) + 1
          k = mod(i+1,3) + 1
          l = mod(i,  3) + 4
          if ( i == 1 ) l = 6
          if ( i == 2 ) l = 5
          dhdsp1(1,i) =  1.0d0/3.0d0
          dhdsp2(1,i) =  1.0d0/3.0d0
          dhdsp1(2,i) = -1.0d0/3.0d0 * (sp1(j)+sp1(k))
          dhdsp2(2,i) = -1.0d0/3.0d0 * (sp2(j)+sp2(k))
          dhdsp1(3,i) =  1.0d0/2.0d0 * (sp1(j)*sp1(k)-sp1(l)**2)
          dhdsp2(3,i) =  1.0d0/2.0d0 * (sp2(j)*sp2(k)-sp2(l)**2)
        end do
        do i = 4,6
          k = mod(i+1,3) + 1
          l = mod(i,  3) + 4
          m = mod(i+1,3) + 4
          if ( i == 5 ) k = 2
          if ( i == 6 ) k = 1
          dhdsp1(2,i) = 2.0d0/3.0d0 * sp1(i)
          dhdsp2(2,i) = 2.0d0/3.0d0 * sp2(i)
          dhdsp1(3,i) = sp1(l)*sp1(m) - sp1(k)*sp1(i)
          dhdsp2(3,i) = sp2(l)*sp2(m) - sp2(k)*sp2(i)
        end do
c                                                        ---- d(sp)/d(s)
        do i = 1,6
          do j = 1,6
            dsdsp1(i,j) = ctp1(i,j)
            dsdsp2(i,j) = ctp2(i,j)
          end do
        end do
c                                                       ---- d(fai)/d(s)
        dfads = 0.0d0
        xx1 = 0.0d0
        xx2 = 0.0d0
        do l = 1,6
          do j = 1,3
            do k = 1,6
              xx1(j,l) = xx1(j,l) + dhdsp1(j,k)*dsdsp1(k,l)
              xx2(j,l) = xx2(j,l) + dhdsp2(j,k)*dsdsp2(k,l)
            end do
            dfads(l) = dfads(l) +
     &                 dfadhp1(j)*xx1(j,l) + 
     &                 dfadhp2(j)*xx2(j,l)
          end do
        end do
c                                                        ---- d(se)/d(s)
        dsedfa = fai**(ami-1.0d0) / am / dc**ami
        do i = 1,6
          dseds(i) = dsedfa * dfads(i)
        end do
c
      end if
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate yield function
c
      subroutine jancae_yld2004_18p_yf ( ctp1,ctp2,s,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,se )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ctp1(6,6),ctp2(6,6),s(6)
      dimension sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3)
c
      call jancae_yld2004_18p_yfsub ( ctp1,s,pi,sp1,psp1,hp1,cetpq1 )
      call jancae_yld2004_18p_yfsub ( ctp2,s,pi,sp2,psp2,hp2,cetpq2 )
c                                                    ---- yield function
      fai = 0.0d0
      do i = 1,3
        do j = 1,3
          fai = fai + (abs(psp1(i)-psp2(j)))**am
        end do
      end do
c                                                 ---- equivalent stress
      se = (fai/dc) ** ami
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate yield function2
c
      subroutine jancae_yld2004_18p_yfsub ( ctp,s,pi,sp,psp,hp,cetpq )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ctp(6,6),s(6)
      dimension sp(6),psp(3),hp(3)
c
c                                         ---- linear-transformed stress
      call jancae_mv ( sp,ctp,s,6,6 )
c                                                  ---- invariants of sp
      hp(1) = (sp(1)+sp(2)+sp(3)) / 3.0d0
      hp(2) = (sp(5)**2+sp(6)**2+sp(4)**2 - 
     1         sp(2)*sp(3)-sp(3)*sp(1)-sp(1)*sp(2)) / 3.0d0
      hp(3) = (2.0d0*sp(5)*sp(6)*sp(4)+sp(1)*sp(2)*sp(3) -
     1         sp(1)*sp(6)**2-sp(2)*sp(5)**2-sp(3)*sp(4)**2) / 2.0d0
c                           ---- coefficients of characteristic equation
      hpq = sqrt(hp(1)**2 + hp(2)**2 + hp(3)**2)
      if ( hpq .gt. 1.0e-16 ) then
        cep = hp(1)**2 + hp(2)
        ceq = (2.0d0*hp(1)**3+3.0d0*hp(1)*hp(2)+2.0d0*hp(3)) / 2.0d0
        cetpq = ceq / cep**(3.0d0/2.0d0)
        if ( cetpq >  1.0d0 ) cetpq =  1.0d0
        if ( cetpq < -1.0d0 ) cetpq = -1.0d0
        cet = acos(cetpq)
c                                           ---- principal values of sp1
        psp(1) = 2.0d0*sqrt(cep)*cos(cet/3.0d0) + hp(1)
        psp(2) = 2.0d0*sqrt(cep)*cos((cet+4.0d0*pi)/3.0d0) + hp(1)
        psp(3) = 2.0d0*sqrt(cep)*cos((cet+2.0d0*pi)/3.0d0) + hp(1)
      else
        cetpq = 0.0
        do i = 1,3
          psp(i) = 0.0
        end do
      end if
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate d(psp)/d(hp)
c
      subroutine jancae_yld2004_18p_dpsdhp ( i,psp,hp,dpsdhp )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension psp(3),hp(3),dpsdhp(3,3)
c
      dummy = psp(i)**2-2.0d0*hp(1)*psp(i) - hp(2)
      dpsdhp(i,1) = psp(i)**2/dummy
      dpsdhp(i,2) = psp(i) / dummy
      dpsdhp(i,3) = 2.0d0/3.0d0/dummy
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate multiplication of matrix and matrix
c     [a]=[b][c]
c
      subroutine jancae_mm (a,b,c,na1,na2,nbc)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(na1,na2),b(na1,nbc),c(nbc,na2)
c
      a = 0.0d0
      do i = 1,na1
        do j = 1,na2
          do k = 1,nbc
            a(i,j) = a(i,j) + b(i,k)*c(k,j)
          end do
        end do
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