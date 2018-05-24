        module diffuse_two_mod
        implicit none
        private

        public :: laplacian

        contains
    
#ifdef constant_laplacian
	subroutine laplacian
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	subroutine to find diffusion term
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : imt, jmt, i, j, k, gdx, gdy, rdx, rdy, u, v, h, taum
        use param_mod, only : diffuse, diffuse_uu, diffuse_vv, dphi, dthe
        use param_mod, only : re, theta
        implicit none
	real r2cos2theta, r2
	real rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6
	real rtemp7, rtemp8, rtemp9, rtemp10, rtemp11
	real rtemp12, rtemp13, rtemp14, rtemp15, rtemp16
	real rtemp17, rtemp18, rtemp20, rtemp21, rtemp22
        integer :: ip, im, jp, jm, ip2, im2, jp2, jm2
        real :: Ahx, Ahy, rH1, rH2, rH3, rH4

	 r2cos2theta = Re*Re*cos(theta(j))*cos(theta(j))
	 r2cos2theta = 1.0/r2cos2theta
	 r2 = Re*Re
	 r2 = 1.0/r2
	
	 ip = i + 1
	 im = i - 1
        if (im<1) im = 1 !Prajeesh  
	 jp = j + 1
	 jm = j - 1
        if (jm<1) jm = 1 !Prajeesh  
	 ip2 = i + 2
	 im2 = i - 2
        if (im2<1) im2 = 1 !Prajeesh  
	 jp2 = j + 2
	 jm2 = j - 2
        if (jm2<1) jm2 = 1 !Prajeesh  
        
c =========================== diffuse_uu ============================
	 rtemp1 = (gdx(i+1) - gdx(i))
	 rtemp2 = (gdy(j+1) - gdy(j))
	 Ahx = diffuse *max(rtemp1, rtemp2) ! taking Ah here
	 Ahy = diffuse *max(rtemp1, rtemp2) ! taking Ah here


	 dphi = (rdx(i+1) + rdx(i))/2.0
	 rH1 = sum_2ptd (3,ip,j,i,j) /sum_2ptd (31,ip,j,i,j)
	 rH2 = sum_2ptd (3,im,j,i,j) /sum_2ptd (31,im,j,i,j)
	 rtemp1 = (u(i+1,j,k,taum) - u(i,j,k,taum))/dphi
	 dphi = (rdx(i-1) + rdx(i))/2.0
	 rH3 = sum_2ptd (3,im,j,i,j) /sum_2ptd (31,im,j,i,j)
	 rH4 = sum_2ptd (3,im2,j,im,j) /sum_2ptd (31,im2,j,im,j)
	 rtemp2 = (u(i,j,k,taum) - u(i-1,j,k,taum))/dphi
	 dphi = rdx(i)
	 rtemp3 = (rtemp1 - rtemp2)*Ahx*rH2/dphi  ! giving Ah here
	 rtemp4 = r2cos2theta*rtemp3  ! reserve

	 dthe = (rdy(j+1) + rdy(j))/2.0
	 rH1 = sum_2ptd (3,i,jp,im,jp) /sum_2ptd (31,i,jp,im,jp)
	 rH2 = sum_2ptd (3,i,j,im,j) /sum_2ptd (31,i,j,im,j)
	 rtemp5 = (u(i,j+1,k,taum) - u(i,j,k,taum))/dthe
	 dthe = (rdy(j-1) + rdy(j))/2.0
	 rH3 = sum_2ptd (3,i,j,im,j) /sum_2ptd (31,i,j,im,j)
	 rH4 = sum_2ptd (31,i,jm,im,jm) /sum_2ptd (31,i,jm,im,jm)
	 rtemp6 = (u(i,j,k,taum) - u(i,j-1,k,taum))/dthe
	 dthe = rdy(j)
	 rtemp7 = (rtemp5 - rtemp6)*Ahy*rH2/dthe ! taking Ah here
	 rtemp8 = r2*rtemp7 ! reserve

	 rH1 = sum_2ptd (3,i,j,im,j) /sum_2ptd (31,i,j,im,j)
	 rtemp9 = u(i,j,k,taum)*rH1*(1-2.0*cos(theta(j))**2) ! reserve
	 dphi = rdx(i)

	 rtemp10 = sum_2ptd(2,i,jp,i,j)/
     &                           sum_2ptd(21,i,jp,i,j)
	 rtemp11 = sum_2ptd(2,im,jp,im,j)/
     &                             sum_2ptd(21,im,jp,im,j)
	 rH1 = sum_2ptd(3,i,j,im,j)/sum_2ptd(31,i,j,im,j)
	 rtemp12 = (rtemp10 - rtemp11)*2.0*sin(theta(j))*rH1/dphi  ! reserve

	 diffuse_uu = rtemp4 + rtemp8 - r2cos2theta*(rtemp9
     &   + rtemp12) * Ahx

	
	
c ==================  diffuse_vv ==========================================
	
         rtemp1 = (gdx(i+1) - gdx(i))
         rtemp2 = (gdy(j+1) - gdy(j))
         Ahx = diffuse *max(rtemp1, rtemp2) ! taking Ah here
         Ahy = diffuse *max(rtemp1, rtemp2) ! taking Ah here

         dphi = (rdx(i+1) + rdx(i))/2.0
         rH1 = sum_2ptd(3,ip,j,ip,jm) /sum_2ptd(31,ip,j,ip,jm)
         rH2 = sum_2ptd(3,i,j,i,jm) /sum_2ptd(31,i,j,i,jm)
         rtemp1 = (v(i+1,j,k,taum) - v(i,j,k,taum))/dphi
         dphi = (rdx(i-1) + rdx(i))/2.0
         rH3 = sum_2ptd(3,i,j,i,jm) /sum_2ptd(31,i,j,i,jm) 
         rH4 = sum_2ptd(3,im,j,im,jm) /sum_2ptd(31,im,j,im,jm)
         rtemp2 = (v(i,j,k,taum) - v(i-1,j,k,taum))/dphi
         dphi = rdx(i)
         rtemp3 = (rtemp1 - rtemp2)*Ahx*rH2/dphi  ! giving Ah here
         rtemp4 = r2cos2theta*rtemp3  ! reserve

         dthe = (rdy(j+1) + rdy(j))/2.0
         rH1 = sum_2ptd (3,i,jp,i,j) /sum_2ptd (31,i,jp,i,j)
         rH2 = sum_2ptd (3,i,j,i,jm) /sum_2ptd (31,i,j,i,jm)
         rtemp5 = (v(i,j+1,k,taum) - v(i,j,k,taum))/dthe
         dthe = (rdy(j-1) + rdy(j))/2.0
         rH3 = sum_2ptd (3,i,j,i,jm) /sum_2ptd (31,i,j,i,jm)
         rH4 = sum_2ptd (3,i,jm,i,jm2) /sum_2ptd (31,i,jm,i,jm2)
         rtemp6 = (v(i,j,k,taum) - v(i,j-1,k,taum))/dthe
         dthe = rdy(j)
         rtemp7 = (rtemp5 - rtemp6)*Ahy*rH2/dthe ! taking Ah here
         rtemp8 = r2*rtemp7 ! reserve

         rH1 = sum_2ptd (3,i,j,i,jm) /sum_2ptd (31,i,j,i,jm)
         rtemp9 = v(i,j,k,taum)*rH1*(1-2.0*cos(theta(j))**2) ! reserve
         dphi = rdx(i)

         rtemp10 = h(i,j,k,taum)*sum_2ptd(1,i,j,ip,j)/
     &                           sum_2ptd(11,i,j,ip,j)
         rtemp11 = h(i,j-1,k,taum)*sum_2ptd(1,i,jm,ip,jm)/
     &                             sum_2ptd(11,i,jm,ip,jm)
	 rH1 = sum_2ptd(3,i,jm,i,j)/sum_2ptd(31,i,jm,i,j)
         rtemp12 = (rtemp10 - rtemp11)*2.0*sin(theta(j))*rH1/dphi  ! reserve

         diffuse_vv = rtemp4 + rtemp8 - r2cos2theta*(rtemp9
     &   - rtemp12) * Ahy


	return
	end subroutine laplacian
#endif



	function sum_4ptd (index, i1,j1, i2,j2, i3,j3, i4,j4)

        use size_mod, only : u, v, h, rkmu, rkmv, rkmh, k, taum

        implicit none

        real :: sum_4ptd
        integer, intent(in) :: index, i1, j1, i2, j2, i3, j3, i4, j4

!#include "size.h"
!#include "polar.h"
!#include "param.h"

        
        if (index .eq. 1) then
        sum_4ptd = (u(i1,j1,k,taum) + u(i2,j2,k,taum) + 
     &             u(i3,j3,k,taum) + u(i4,j4,k,taum))
        endif

        if (index .eq. 2) then
        sum_4ptd = (v(i1,j1,k,taum) + v(i2,j2,k,taum) + 
     &             v(i3,j3,k,taum) + v(i4,j4,k,taum))
        endif

        if (index .eq. 3) then
        sum_4ptd = (h(i1,j1,k,taum) + h(i2,j2,k,taum) + 
     &             h(i3,j3,k,taum) + h(i4,j4,k,taum))
        endif

        if (index .eq. 11) then
        sum_4ptd = (rkmu(i1,j1) + rkmu(i2,j2) + rkmu(i3,j3) + rkmu(i4,j4))
        sum_4ptd = max(1.0,sum_4ptd)
        endif

        if (index .eq. 21) then
        sum_4ptd = (rkmv(i1,j1) + rkmv(i2,j2) + rkmv(i3,j3) + rkmv(i4,j4))
        sum_4ptd = max(1.0,sum_4ptd)
        endif

        if (index .eq. 31) then
        sum_4ptd = (rkmh(i1,j1) + rkmh(i2,j2) + rkmh(i3,j3) + rkmh(i4,j4))
        sum_4ptd = max(1.0,sum_4ptd)
        endif


        end function sum_4ptd


        function sum_2ptd (index, i1,j1, i2,j2)
            use size_mod, only: u, v, h, taum, k, rkmu, rkmv, rkmh
        implicit none
        
        real :: sum_2ptd
        integer, intent(in) :: index, i1, j1, i2, j2

        if (index .eq. 1) then
        sum_2ptd = (u(i1,j1,k,taum) + u(i2,j2,k,taum))
        endif

        if (index .eq. 2) then
        sum_2ptd = (v(i1,j1,k,taum) + v(i2,j2,k,taum))
        endif

        if (index .eq. 3) then
        sum_2ptd = (h(i1,j1,k,taum) + h(i2,j2,k,taum))
        endif

        if (index .eq. 11) then
        sum_2ptd = (rkmu(i1,j1) + rkmu(i2,j2))
        sum_2ptd = max(1.0,sum_2ptd)
        endif

        if (index .eq. 21) then
        sum_2ptd = (rkmv(i1,j1) + rkmv(i2,j2))
        sum_2ptd = max(1.0,sum_2ptd)
        endif

        if (index .eq. 31) then
        sum_2ptd = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptd = max(1.0,sum_2ptd)
        endif


        end function sum_2ptd

        end module diffuse_two_mod

