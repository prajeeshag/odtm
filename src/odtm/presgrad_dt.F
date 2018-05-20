        module presgrad_mod
        use advec_mod, only : sum_2pt
        implicit none
        private

        public :: pressure_integral, density_nonlinear
        contains

        subroutine pressure_integral (loop_ind)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	subroutine to calculated presure gradient for each layer
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only: i, j, k, km, taun, pres_gradu, pres_gradv, h, eta
            use size_mod, only: rkmh, dxu, dyv, rkmu, rkmv  
            use param_mod, only : denss, g
        implicit none
            integer, intent(in) :: loop_ind 
            integer :: ip, im, jp, jm, kk, ii, jj
            real :: rsum1, rsum2, rsum, rsumeta, rsumeta_init, rho_km, rho_kk
            real :: rtemp1, rho_k, rwall

	pres_gradu(k) = 0.0
	pres_gradv(k) = 0.0
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1
	
!===================== TERM-I ======================================
	rsum1 = 0.0
	rsum2 = 0.0
	rsum = 0.0
	rsumeta = 0.0
	rsumeta_init = 0.0

	do kk=1,km-1
	rho_km = rho_2pt (i,j,km)
	rho_kk = rho_2pt (i,j,kk)
	rsumeta = rsumeta + (rho_km - rho_kk )*h(i,kk,j,taun) /rho_km
	rsum1 = ((rho_km - rho_kk )/rho_km)*h(i,kk,j,taun)
	rho_km = rho_2pt (i,j,km)
	rho_kk = rho_2pt (i,j,kk)
	rsum2 = ((rho_km - rho_kk )/rho_km)*h(im,kk,j,taun)
	rsum = rsum + (rsum1 - rsum2)/dxu(i,j)
	enddo
	eta(i,k,j,taun) = rsumeta*rkmh(i,j)
	rtemp1 = sum_2pt (3,i,j,im,j) /sum_2pt (31,i,j,im,j)
	rho_k = rho_2pt (i,j,k)
	rsum1 = rsum*rho_k*rtemp1
	pres_gradu(k) = rsum1*g 

!==================== TERM-I =======================================
	
!==================== TERM-II ======================================
	rsum1 = 0.0
	rsum2 = 0.0
	rsum = 0.0
	do kk=1,k-1,1
	rho_k = rho_2pt (i,j,k)
	rho_kk = rho_2pt(i,j,kk)
	rsum1 = ( rho_k - rho_kk)*h(i,kk,j,taun)
	rho_k = rho_2pt (i,j,k)
	rho_kk = rho_2pt(i,j,kk)
	rsum2 = ( rho_k - rho_kk)*h(im,kk,j,taun)
	rsum = rsum + (rsum1 - rsum2)/dxu(i,j)
	enddo
	rtemp1 = sum_2pt (3,im,j,i,j) /sum_2pt (31,im,j,i,j)
	rsum2 = rsum*rtemp1
	rsum2 = rsum2*g
	pres_gradu(k) = pres_gradu(k) - rsum2 

!==================== TERM-II ======================================
	
!===================== TERM-I ======================================
	rsum1 = 0.0
	rsum2 = 0.0
	rsum = 0.0
	rsumeta = 0.0
	rsumeta_init = 0.0

	do kk=1,km-1
	rho_km = rho_2pt (i,j,km)
	rho_kk = rho_2pt (i,j,kk)
c        if (rho_km .le. rho_kk) rho_km = rho_kk+0.1
	rsum1 = ((rho_km - rho_kk )/rho_km)*h(i,kk,j,taun)
	rho_km = rho_2pt (i,j,km)
	rho_kk = rho_2pt (i,j,kk)
c        if (rho_km .le. rho_kk) rho_km = rho_kk+0.1
	rsum2 = ((rho_km - rho_kk )/rho_km)*h(i,kk,jm,taun)
	rsum = rsum + (rsum1 - rsum2)/dyv(i,j)
	enddo
	rtemp1 = sum_2pt (3,i,j,i,jm) /sum_2pt (31,i,j,i,jm)
	rho_k = rho_2pt (i,j,k)
	rsum1 = rsum*rho_k*rtemp1
	pres_gradv(k) = rsum1*g 
	

!==================== TERM-I =======================================
	
!==================== TERM-II ======================================
	rsum1 = 0.0
	rsum2 = 0.0
	rsum = 0.0
	do kk=1,k-1,1
	rho_k = rho_2pt (i,j,k)
	rho_kk = rho_2pt(i,j,kk)
	rsum1 = ( rho_k - rho_kk)*h(i,kk,j,taun)
	rho_k = rho_2pt (i,j,k)
	rho_kk = rho_2pt(i,j,kk)
	rsum2 = ( rho_k - rho_kk)*h(i,kk,jm,taun)
	rsum = rsum + (rsum1 - rsum2)/dyv(i,j)
	enddo
	rtemp1 = sum_2pt (3,i,j,i,jm) /sum_2pt (31,i,j,i,jm)
	rsum2 = rsum*rtemp1
	rsum2 = rsum2*g
	pres_gradv(k) =  pres_gradv(k) - rsum2
	
!=================== TERM-II ==========================================

! Thermodynamic term treatement

        rsum = 0.0
        rwall = 1.0
        do ii=im,ip
        do jj=jm,jp
        rwall = rkmu(ii,jj) * rkmv(ii,jj) * rkmh(ii,jj)
        enddo
        enddo
	do kk=k,km-1
        rho_km = rho_2pt (i,j,km)
        rho_kk = rho_2pt (i,j,kk)
        if (rho_km .lt. 1000.0 .or. rho_kk .lt. 1000.0) rwall = 0.0
        rsum1 = ((rho_km - rho_kk )/rho_km) 
        rho_km = rho_2pt (im,j,km)
        rho_kk = rho_2pt (im,j,kk)
        if (rho_km .lt. 1000.0 .or. rho_kk .lt. 1000.0) rwall = 0.0
        rsum2 = ((rho_km - rho_kk )/rho_km) 
        rsum = rsum + (rsum1 - rsum2)/dxu(i,j)
        enddo
        rtemp1 = sum_2pt (3,i,j,im,j) /sum_2pt (31,i,j,im,j)
        rho_k = rho_2pt (i,j,k)
        rsum1 = rsum*rho_k*rtemp1*rtemp1  ! rtemp1 is twice because
                                          ! depth integrated form
        pres_gradu(k) = pres_gradu(k)  + rsum1*g * rwall 
     &  * rkmu(i,j) * rkmv(i,j) * rkmh(i,j)


        rsum = 0.0
        do kk=k,km-1
        rho_km = rho_2pt (i,j,km)
        rho_kk = rho_2pt (i,j,kk)
        if (rho_km .lt. 1000.0 .or. rho_kk .lt. 1000.0) rwall = 0.0
        rsum1 = ((rho_km - rho_kk )/rho_km)
        rho_km = rho_2pt (i,jm,km)
        rho_kk = rho_2pt (i,jm,kk)
        if (rho_km .lt. 1000.0 .or. rho_kk .lt. 1000.0) rwall = 0.0
        rsum2 = ((rho_km - rho_kk )/rho_km)
        rsum = rsum + (rsum1 - rsum2)/dyv(i,j)
        enddo
        rtemp1 = sum_2pt (3,i,j,i,jm) /sum_2pt (31,i,j,i,jm)
        rho_k = rho_2pt (i,j,k)
        rsum1 = rsum*rho_k*rtemp1*rtemp1  ! rtemp1 is twice because
                                          ! depth integrated form
        pres_gradv(k) = pres_gradv(k)  + rsum1*g * rwall


c        pvort(i,k,j) = rsum1*g*rkmv(i,j)
! End of Thermodynamic term
        
	
c         if (rho_2pt(i,j,k+1) .gt. rho_2pt (i,j,k))  then
	 pres_gradu(k) = pres_gradu(k)*rkmu(i,j)
	 pres_gradv(k) = pres_gradv(k)*rkmv(i,j) 
c         else
c         pres_gradu(k) = 0.0
c         pres_gradv(k) = 0.0
c         endif
c	write(*,*) pres_gradu(k), pres_gradv(k)
	
	return
	end subroutine pressure_integral 


	function rho_2pt (i1,j1,k1)

        use size_mod, only : taun, t, zdz, rkmh
        use param_mod, only : denss
        implicit none
        
        double precision :: rho_2pt
        integer, intent(in) :: i1, j1, k1
        real :: VolExpan, SalContr, rkappa_vol, rTemp, rSalt, rDep

	rho_2pt = 0.0
	VolExpan = -0.15
        SalContr =  0.78
        rkappa_vol = 0.0045
        rTemp = t(i1,k1,j1,1,taun)
        rSalt = t(i1,k1,j1,2,taun)
        rDep = zdz(k1)
	rho_2pt = density_nonlinear (rTemp, rSalt, rDep) !rho(km) + rho_2pt

        denss(i1,k1,j1) = rho_2pt * rkmh(i1,j1)
        
	end function rho_2pt


        function density_nonlinear (rTemp, rSalt, rDepth)

        implicit none

        real, intent(in) :: rTemp, rSalt, rDepth
        real :: density_nonlinear, rKtsp

        real :: rT1, rT2, rT3, rT4, rT5, rS1, rS2, rS1b2, rP1, rP2
        real :: rTempt, rSaltt, rho_0, rT, rS

         rT1 = rTemp**1
         rT2 = rTempt**2 !rTempt is not assigned !!??? Prajeesh  
         rT3 = rTempt**3
         rT4 = rTempt**4
         rT5 = rTempt**5
         rS1 = rSalt**1
         rS2 = rSaltt**2 !rSaltt is not assigned !!?? Prajeesh
         rS1b2 = rSalt**1.5
         rP1 = rDepth*0.1  ! in bars approximately
         rP2 = rP1**2  
! bar (1 millibar = 0.1 kPa)
! http://docs.bluerobotics.com/calc/pressure-depth/
! 10m sea water depth = 100.45 kPa = 1004.5 millibar atm press = 1 bar)
! 100m sea water depth = 1004.54 kPa = 10045.5 millibar atm press =
! 10.045 bar
! 1000m sea water depth = 10045.44 kPa = 100454.4 millibar atm press =
! 100.4544 bar
        
        
         rho_0 = (
     1          +999.842594           +  6.793952e-2*rT1
     2          -9.095290e-3*rT2      +  1.001685e-4*rT3
     3          -1.120083e-6*rT4      +  6.536332e-9*rT5
     4          +8.24493e-1*rS1       -  4.0899e-3*rT1*rS1
     5          +7.64380e-5*rT2*rS1   -  8.2467e-7*rT3*rS1
     6          +5.38750e-9*rT4*rS1   -  5.72466e-3*rS1b2
     7          +1.02270e-4*rT1*rS1b2 -  1.6546e-6*rT2*rS1b2
     8          +4.83140e-4*rS2 )


        rKtsp = +19652.21
     1          +148.4206    *rT        -  2.327105*rT2
     2          +1.360477e-2 *rT3       -  5.155288e-5*rT4
     3          +3.239908    *rP1       +  1.43713e-3*rT1*rP1    
     4          +1.16092e-4  *rT2*rP1   -  5.77905e-7*rT3*rP1
     5          +8.50935e-5  *rP2       -  6.12293e-6*rT1*rP2
     6          +5.2787e-8   *rT2*rP2
     7          +54.6746*rS1            -  0.603459*rT1*rS1
     8          +1.09987e-2*rT2*rS      -  6.1670e-5*rT3*rS1
     9          +7.944e-2*rS1b2         +  1.6483e-2*rT*rS1b2
     1          -5.3009e-4*rT2*rS1b2    +  2.2838e-3*rP1*rS1
     2          -1.0981e-5*rT1*rP1*rS1  -  1.6078e-6*rT2*rP1*rS1
     3          +1.91075e-4*rP1*rS1b2   -  9.9348e-7*rP2*rS1         
     4          +2.0816e-8*rT1*rP2*rS1  +  9.1697e-10*rT2*rP2*rS1


        density_nonlinear = rho_0*(1.0-rP1/rKtsp)**(-1) 
        
        end function density_nonlinear
        end module presgrad_mod
