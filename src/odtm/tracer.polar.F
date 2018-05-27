        module tracer_mod

        use advec_mod, only : sum_2pt
        implicit none

        logical, public :: rgm_zero_tracer_adv = .true.
        logical :: relax_tracer = .false.


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
        use size_mod, only : gdx, gdy, h, i, j, k, loop, month, rdx
        use size_mod, only : rdy, rkmh, salt_read, t, taum, taum1, taun, taun1
        use size_mod, only : taup1, temp_read, taup, kclim, theta, nn
        use size_mod, only : rdxu, rdyv, rdxh, rdyh, kmaxMYM
        use param_mod, only : day2sec, diffuse_h, diffuse_my, diffuse_tr, dphi
        use param_mod, only : dtts, dt, rad2deg

        use interp_extrap_initial_mod, only : interp_extrap_initial 

        implicit none
        real relax_tr(nn), tempin(201), saltin(201)
        
        integer :: nt, ip, im, jp, jm
        real :: tracedvbdy, tracedubdx, rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6
        real :: tempout, saltout, hstart, hend, zero_out_tradv
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

        zero_out_tradv = 1.0

        if (rgm_zero_tracer_adv) then
            hstart = 0
            kmax = kmaxMYM
            do irepeat = 1,k-1
                hstart = hstart + h(i,j,irepeat,taun)
            enddo

            hend = h(i,j,k,taun) + hstart

            if (hend .le. 5.0*kmax) then 
                zero_out_tradv = 0.0
            endif
        endif

	tracedubdx = 0.0
	tracedvbdy = 0.0
	
        if (zero_out_tradv == 1.0) then
	rtemp1 = sum_2ptr (1,ip,j,i,j,nt) /
     &                           sum_2ptr (11,ip,j,i,j,nt)
	rtemp2 = sum_2ptr (1,im,j,i,j,nt) /
     &                           sum_2ptr (11,im,j,i,j,nt)

	tracedubdx = (rtemp1 - rtemp2)* (sum_2pt (1,ip,j,i,j)/
     &                              sum_2pt (11,ip,j,i,j) ) * rdxh(i,j)


	rtemp3 = sum_2ptr (1,i,j,i,jp,nt) /
     &                            sum_2ptr (11,i,j,i,jp,nt)
	rtemp4 = sum_2ptr (1,i,j,i,jm,nt) /
     &                           sum_2ptr (11,i,j,i,jm,nt)
	tracedvbdy = ( rtemp3 - rtemp4 )* (sum_2pt (2,i,jp,i,j)/
     &                              sum_2pt (21,i,jp,i,j) ) * rdyh(i,j) 

        endif



	rtemp1 = rkmh(ip,j)*(t(ip,j,k,nt,taum) - t(i,j,k,nt,taum) )*rdxu(ip,j)
	rtemp2 = rkmh(im,j)*(t(i,j,k,nt,taum) - t(im,j,k,nt,taum) )*rdxu(i,j)
	rtemp3 = (rtemp1 - rtemp2)*rdxh(i,j) !reserve

	rtemp4 = rkmh(i,jp)*(t(i,jp,k,nt,taun) - t(i,j,k,nt,taun))*rdyv(i,jp)
	rtemp5 = rkmh(i,jm)*(t(i,j,k,nt,taun) - t(i,jm,k,nt,taun))*rdyv(i,j)
	rtemp6 = (rtemp4 - rtemp5)*rdyh(i,j)
	
	diffuse_h = rtemp3 + rtemp6

        kmax = kclim
        
        relax_tr(:) = 0.

        if (relax_tracer) then
        do kk=1,kmax
         tempin(kk) = temp_read(i,j,kk,month)
         saltin(kk) = salt_read(i,j,kk,month)
         enddo

        call interp_extrap_initial (i,j,k,kmax,tempin,
     &   saltin,tempout,saltout)

        relax_tr(1) = (tempout - t(i,j,k,1,taun))/
     &      (30.0*day2sec)
        relax_tr(2) = (saltout - t(i,j,k,2,taun))/
     &      (30.0*day2sec)

        endif

! upadting MYM 

        relax_on_off = 1.0
#include "cparam.h"  ! <== This file overrides above parameter setup
        
        
	     rtemp1 = rdx(i) * rad2deg
         rtemp2 = rdy(j) * rad2deg


        t(i,j,k,nt,taup) = t(i,j,k,nt,taum)*rkmh(i,j) + (
     &  -(
     &     tracedubdx * zero_out_tradv
     &     + tracedvbdy * zero_out_tradv
     &     - diffuse_h * diffuse_tr * max(rtemp1, rtemp2)
     &     - relax_tr(nt) * relax_on_off
     &      )
     &  ) *dtts*rkmh(i,j)

	enddo
        
	return
	end subroutine tracer



         function sum_2ptr (index, i1,j1, i2,j2,nt)
            use size_mod, only : k, t, rkmh, taun
        implicit none
        real :: sum_2ptr
        integer, intent(in) :: index, i1,j1, i2,j2,nt

        if (index .eq. 1) then
        sum_2ptr = (t(i1,j1,k,nt,taun) + t(i2,j2,k,nt,taun))
        endif

        if (index .eq. 11) then
        sum_2ptr = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptr = max(1.0,sum_2ptr)
        endif

        end function sum_2ptr 

        end module tracer_mod

