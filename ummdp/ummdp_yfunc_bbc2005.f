c--------------------------------------------------------------(bbc2005)
c
c     BBC 2005 yield function and its differentials
c
	subroutine jancae_bbc2005 ( stress,ndyld,pryld,nreq,
     1                            se,dseds )
c
c-----------------------------------------------------------------------
c
	implicit real*8 (a-h,o-z)
	dimension stress(6)
	dimension s(3),dseds(3),pryld(ndyld)
c
	real*8 a,b,L,M,N,P,Q,R,nn,fir,las,oo
	real*8 th(3),lth(3),pp
	real*8 phi,Al,AA,BB,CC,DD,kk
	real*8 lth1_2,lth12,lth2_3,lth23
	real*8 dlthds(3,3),dphidlth(3),dsedphi
	integer i,j,d,e,k,mm,fact,ii,jj
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
	k = pryld(1+1)
	a = pryld(1+2)
	b = pryld(1+3)
	L = pryld(1+4)
	M = pryld(1+5)
	N = pryld(1+6)
	P = pryld(1+7)
	Q = pryld(1+8)
	R = pryld(1+9)
c
c                                                 ---- equivalent stress
	th(1) = L*s(1) + M*s(2)
	th(2) = sqrt((N*s(1)-P*s(2))**2+s(3)**2)
	th(3) = sqrt((Q*s(1)-R*s(2))**2+s(3)**2)
c
	AA = th(2)**2 + th(1)**2
	BB = 2 * th(2) * th(1)
	CC = th(2)**2 + th(3)**2
	DD = 2 * th(2) * th(3)
c
	mm = int(k/2)
	nn = 0.0d0
	oo = 0.0d0
c
	do i = 0,mm
	  nn = nn + fact(k)/(fact(k-2*i)*fact(2*i))*
     &       AA**(k-2*i)*BB**(2*i)
	  oo = oo + fact(k)/(fact(k-2*i)*fact(2*i))*
     &       CC**(k-2*i)*DD**(2*i)
	end do
c
	fir = 2.0d0 * a * nn
	las = 2.0d0 * b * oo
c
	phi = fir + las
c
	kk = dble(k)
c
	Al = (a*(N+L)**(2*k)+a*(N-L)**(2*k) +
     &   b*(N+Q)**(2*k)+b*(N-Q)**(2*k))**(1.0d0/(2.0d0*kk))
c
	se = phi**(1/(2*kk)) / Al
c
	dseds(:) = 0
c
c                                           ----  1st order differential
	if ( nreq .ge. 1 ) then
c
	  dsedphi = (1.0d0/(2.0d0*kk))*(phi**(1.0d0/(2.0d0*kk)-1.0d0)/Al)
c
	  lth(1) = (L*s(1)+M*s(2))**2
	  lth(2) = (N*s(1)-P*s(2))**2+s(3)**2
	  lth(3) = (Q*s(1)-R*s(2))**2+s(3)**2
c
c                    ----to avoid division by zero or zero of zero power
c 
	  if ( lth(1) < 1e-15 * se**2) then
	    lth(1) = 1e-15 * se**2
	  end if
c
	  if ( lth(2) < 1e-15 * se**2) then
	    lth(2) = 1e-15 * se**2
	  end if
c
	  if ( lth(3) < 1e-15 * se**2) then
	    lth(3) = 1e-15 * se**2
	  end if
c
	  lth1_2 = lth(1)+lth(2)
	  lth12  = lth(1)*lth(2)
	  lth2_3 = lth(2)+lth(3)
	  lth23  = lth(2)*lth(3)
c
	  if (lth1_2 < 1e-15 * se**2) then
	    lth1_2 = 1e-15 * se**2
	  end if
c
	  if (lth2_3 < 1e-15 * se**2) then
	    lth2_3 = 1e-15 * se**2
	  end if
c
	  dphidlth(:) = 0.0d0
c
	  do i = 0,mm  
	    dphidlth(1) = dphidlth(1) + fact(k)/(fact(k-2*i)*fact(2*i))*
     &     (2*a*(i*4**i*lth(2)**i*lth(1)**(i-1)*lth1_2
     &     **(k-2*i)+(k-2*i)*4**i*lth12**i*lth1_2**(-2*i+k-1)))
c
	    dphidlth(2) = dphidlth(2) + fact(k)/(fact(k-2*i)*fact(2*i))*
     &    (2*a*(i*4**i*lth(2)**(i-1)*lth(1)**i*lth1_2
     &    **(k-2*i)+(k-2*i)*4**i*lth12**i*
     &    lth1_2**(-2*i+k-1))
     &    +2*b*(i*4**i*lth(2)**(i-1)*lth(3)**i*lth2_3
     &    **(k-2*i)+(k-2*i)*4**i*lth23**i*
     &    lth2_3**(-2*i+k-1)))
c
	    dphidlth(3) = dphidlth(3) + fact(k)/(fact(k-2*i)*fact(2*i))*
     &     (2*b*(i*4**i*lth(2)**i*lth(3)**(i-1)*lth2_3
     &     **(k-2*i)+(k-2*i)*4**i*lth23**i*lth2_3**(-2*i+k-1)))
	  end do
c
	  dlthds(1,1) =  2 * L * (M*s(2)+L*s(1))
	  dlthds(1,2) =  2 * M * (M*s(2)+L*s(1))
	  dlthds(1,3) =  0.0d0
	  dlthds(2,1) =  2 * N * (N*s(1)-P*s(2))
	  dlthds(2,2) = -2 * P * (N*s(1)-P*s(2))
	  dlthds(2,3) =  2 * s(3)
	  dlthds(3,1) =  2 * Q * (Q*s(1)-R*s(2))
	  dlthds(3,2) = -2 * R * (Q*s(1)-R*s(2))
	  dlthds(3,3) =  2 * s(3)
c
	  do i = 1,3
	    dseds(i) = 0.0d0
	    do j = 1,3
				dseds(i) = dseds(i) + dsedphi*dphidlth(j)*dlthds(j,i)
	    end do
	  end do
c
	end if
c
	return
	end
c
c
c
c----  function to solve factorial ----------------------------(bbc2005)
	integer function fact(n) result(m)
	integer, intent(in) :: n
	m = 1
	if ( n .ge. 1 ) then
	    do i = 1,n
		m = m * i
	    end do
	end if
c
	end function fact
c
c
c