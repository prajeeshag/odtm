	
        module tracer_mod
        use advec_mod, only : sum_2pt
        implicit none
        private
        
        public :: tracer
    
        contains

        subroutine tracer
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	subroutine to solve baroclinic pressure gradient
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : dxu, dyu, gdx, gdy, h, i, j, k, loop, month, rdx
        use size_mod, only : rdy, rkmh, salt_read, t, taum, taum1, taun, taun1
        use size_mod, only : taup1, temp_read, dyv, taup
        use param_mod, only : day2sec, diffuse_h, diffuse_my, diffuse_tr, dphi
        use param_mod, only : dthe, dtts, nn, theta, dt

        implicit none
        real relax_tr(nn), tempin(201), saltin(201)
        
        integer :: nt, ip, im, jp, jm
        real :: tracedvbdy, tracedubdx, rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6
        real :: rtheta_vu, rtheta_vd, tempout, saltout, hstart, hend, zero_out_tradv
        real :: relax_on_off, alpha_rgm, alpha_rgt, beta_mldm, beta_mldt, gama, rdrag
        real :: relax_South, widS, relax_East, widE
        integer :: irepeat, kmax, kk

        do nt = 1,nn
	taum1 = 1
	taun1 = 2
	taup1 = 3
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1

	tracedubdx = 0.0
	tracedvbdy = 0.0
	
	rtemp1 = sum_2ptr (1,ip,j,i,j,nt) /
     &                           sum_2ptr (11,ip,j,i,j,nt)
	rtemp2 = sum_2ptr (1,im,j,i,j,nt) /
     &                           sum_2ptr (11,im,j,i,j,nt)
	dphi = (rdx(i+1)+rdx(i))/2.0

	tracedubdx = (rtemp1 - rtemp2)* (sum_2pt (1,ip,j,i,j)/
     &                              sum_2pt (11,ip,j,i,j) ) /
     &              ( (dxu(i,j)+dxu(ip,j))/2  )


	rtheta_vu = cos((theta(j)+theta(j+1))*0.5)
        rtheta_vd = cos((theta(j)+theta(j-1))*0.5)
	dthe = (rdy(j) + rdy(j+1))/2.0

	rtemp3 = sum_2ptr (1,i,j,i,jp,nt) /
     &                            sum_2ptr (11,i,j,i,jp,nt)
	rtemp4 = sum_2ptr (1,i,j,i,jm,nt) /
     &                           sum_2ptr (11,i,j,i,jm,nt)
	tracedvbdy = ( rtemp3 - rtemp4 )* (sum_2pt (2,i,jp,i,j)/
     &                              sum_2pt (21,i,jp,i,j) ) /
     &                ( (dyv(i,jp)+dyv(i,j))/2 )



	
#ifdef restrt
        if ( loop .eq. 1) dtts = dtts
#else
        if ( loop .eq. 1) dtts = dt
#endif

	rtemp1 = rkmh(ip,j)*
     &    (t(ip,k,j,nt,taum) - t(i,k,j,nt,taum) )/dxu(ip,j)
	rtemp2 = rkmh(im,j)*
     &    (t(i,k,j,nt,taum) - t(im,k,j,nt,taum) )/dxu(i,j)
	rtemp3 = (rtemp1 - rtemp2)/((dxu(ip,j)+dxu(i,j))*0.5) !reserve

	rtemp4 = rkmh(i,jp)*
     &   (t(i,k,jp,nt,taun) - t(i,k,j,nt,taun))/dyv(i,jp)
	dthe = (rdy(j-1) + rdy(j))/2.0
	rtemp5 = rkmh(i,jm)*
     &   (t(i,k,j,nt,taun) - t(i,k,jm,nt,taun))/dyv(i,j)
	rtemp6 = (rtemp4 - rtemp5)/((dyv(i,jp)+dyv(i,j))*0.5)  ! reserve
	
	diffuse_h = rtemp3 + rtemp6

	rtemp1 = (gdx(i+1) - gdx(i))
        rtemp2 = (gdy(j+1) - gdy(j))

        kmax = 201
        do kk=1,kmax
         tempin(kk) = temp_read(i,kk,j,month)
         saltin(kk) = salt_read(i,kk,j,month)
         enddo
        call interp_extrap_initial (i,j,k,kmax,tempin,
     &   saltin,tempout,saltout)


        relax_tr(1) = (tempout - t(i,k,j,1,taun))/
     &      (30.0*day2sec)
        relax_tr(2) = (saltout - t(i,k,j,2,taun))/
     &      (30.0*day2sec)

        hstart = 0
        kmax = 51
        do irepeat = 1,k-1
        hstart = hstart + h(i,irepeat,j,taun)
        enddo
        hend = h(i,k,j,taun) + hstart
! upadting MYM 
        if (hend .le. 5.0*kmax) then 
        zero_out_tradv = 0.0
        else
        zero_out_tradv = 1.0
        endif

        relax_on_off = 1.0
#include "cparam.h"  ! <== This file overrides above parameter setup
        
c        pvort(i,k,j) = zero_out_tradv

	 rtemp1 = (gdx(i+1) - gdx(i))
         rtemp2 = (gdy(j+1) - gdy(j))


        t(i,k,j,nt,taup) = t(i,k,j,nt,taum)*rkmh(i,j) + (
     &  -(
     &     tracedubdx * zero_out_tradv
     &     + tracedvbdy * zero_out_tradv
     &     - diffuse_h * diffuse_tr * max(rtemp1, rtemp2)
     &     - relax_tr(nt) * relax_on_off
     &      )
     &  ) *dtts*rkmh(i,j)

!     &     + max(0.0, we(i,k,j)*0.0)*
!     &          (t(i,k,j,nt,taun) - t(i,k+1,j,nt,taun))/
!     &       (he(k))
	
	enddo
        
	return
	end subroutine tracer



         function sum_2ptr (index, i1,j1, i2,j2,nt)
            use size_mod, only : k, t, rkmh, taun
        implicit none
        double precision :: sum_2ptr
        integer, intent(in) :: index, i1,j1, i2,j2,nt

        if (index .eq. 1) then
        sum_2ptr = (t(i1,k,j1,nt,taun) + t(i2,k,j2,nt,taun))
        endif

        if (index .eq. 11) then
        sum_2ptr = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptr = max(1.0,sum_2ptr)
        endif

        end function sum_2ptr 

        end module tracer_mod

