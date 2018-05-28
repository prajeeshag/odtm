        module couple_mod
        use interp_extrap_mod, only : interp_extrap
        use mpp_mod, only : mpp_pe

        implicit none

        private

        public :: couple_rgmld

        contains

        subroutine couple_rgmld 

        ! linear interpolation for wind forcing

        use size_mod, only : isc, iec, jsc, jec
        use size_mod, only : taum, taun, taup, taus, tracedubdxmld, tracedvbdymld
        use size_mod, only : loop, km, i, j, k, rkmh, t, u, rkmu, v
        use size_mod, only : rkmv, h, nn, temp, rkmu, rkmv
        use size_mod, only : salt, uvel, vvel, rdx, rdxh, rdy, rdyh
        use size_mod, only : temp_read, salt_read, month, theta, kmaxMYM
        use param_mod, only : day2sec, dt, sec2day

        implicit none

        real :: tempout_taum, tempout_taup
        real :: saltout_taum, saltout_taup
        real :: uvelout_taum, uvelout_taup
        real :: vvelout_taum, vvelout_taup

        real :: tempout_taumC, tempout_taupC
        real :: saltout_taumC, saltout_taupC
        real :: uvelout_taumC, uvelout_taupC
        real :: vvelout_taumC, vvelout_taupC

        real :: tempout_taumL, tempout_taupL
        real :: saltout_taumL, saltout_taupL
        real :: uvelout_taumL, uvelout_taupL
        real :: vvelout_taumL, vvelout_taupL

        real :: tempout_taumB, tempout_taupB
        real :: saltout_taumB, saltout_taupB
        real :: uvelout_taumB, uvelout_taupB
        real :: vvelout_taumB, vvelout_taupB

        real :: alpha_rgm, beta_mldm
        real :: alpha_rgt, beta_mldt

        real :: diffuse_temp, diffuse_salt, rtheta_vu, rtheta_vd
        real :: relax_on_off, gama, diffuse_MY, diffuse_tr
        real :: relax_south, widS, relax_East, widE, rdrag
        real :: rtemp1, rtemp2, rtemp3, rtemp4
        integer :: hstart, hend, ii, jj, kk, tau_index
        integer :: im, jm, irepeat, nt, ip, jp, kmax

        tracedubdxmld(:,:) = 0.0
        tracedvbdymld(:,:) = 0.0
        
        !taum = 1
        !taun = 2
        !taup = 3
        !taus = 4

        alpha_rgm = 1.0  *      mod(loop,2)
        alpha_rgt = 1.0  * (1 - mod(loop,2))
        beta_mldm = 1.0  *      mod(loop,2)
        beta_mldt = 1.0  * (1 - mod(loop,2))

#include "cparam.h"  ! <== This file overrides above parameter setup

            kmax  = kmaxMYM

        do i=isc,iec
        do j=jsc,jec
        do k=1,km-1

        im = i-1
        jm = j-1

        tau_index = 1  
        call interp_extrap (i,j,k,kmax,tau_index,
     &  tempout_taumC,saltout_taumC,uvelout_taumC,vvelout_taumC)

        tau_index = 2
        call interp_extrap (i,j,k,kmax,tau_index,
     &  tempout_taupC,saltout_taupC,uvelout_taupC,vvelout_taupC)
        
        tau_index = 1  
        call interp_extrap (im,j,k,kmax,tau_index,
     &  tempout_taumL,saltout_taumL,uvelout_taumL,vvelout_taumL)

        tau_index = 2
        call interp_extrap (im,j,k,kmax,tau_index,
     &  tempout_taupL,saltout_taupL,uvelout_taupL,vvelout_taupL)
        
        tau_index = 1  
        call interp_extrap (i,jm,k,kmax,tau_index,
     &  tempout_taumB,saltout_taumB,uvelout_taumB,vvelout_taumB)

        tau_index = 2
        call interp_extrap (i,jm,k,kmax,tau_index,
     &  tempout_taupB,saltout_taupB,uvelout_taupB,vvelout_taupB)

        
        tempout_taup = tempout_taupC
        tempout_taum = tempout_taumC
        saltout_taup = saltout_taupC
        saltout_taum = saltout_taumC

        uvelout_taup = (uvelout_taupC + uvelout_taupL)/
     &                 max(1.0, rkmh(i,j) + rkmh(im,j) )
        uvelout_taum = (uvelout_taumC + uvelout_taumL)/
     &                 max(1.0, rkmh(i,j) + rkmh(im,j) )
       
      
        vvelout_taup = (vvelout_taupC + vvelout_taupB)/
     &                 max(1.0, rkmh(i,j) + rkmh(i,jm) )
        vvelout_taum = (vvelout_taumC + vvelout_taumB)/
     &                 max(1.0, rkmh(i,j) + rkmh(i,jm) )
       
        

! upadting RG model
! Tbar = (1-alpha_rg)T + alpha_rg (T_taup - T_taun)


        t(i,j,k,1,taup) =  (t(i,j,k,1,taup)
     &  + alpha_rgt*(tempout_taup - tempout_taum)) * rkmh(i,j)
        t(i,j,k,2,taup) =  (t(i,j,k,2,taup)  
     &  + alpha_rgt*(saltout_taup - saltout_taum)) * rkmh(i,j)
        u(i,j,k,taup) =  (u(i,j,k,taup)
     &   + alpha_rgm*(uvelout_taup - uvelout_taum )) * rkmu(i,j)
        v(i,j,k,taup) =  (v(i,j,k,taup)
     &   + alpha_rgm*(vvelout_taup - vvelout_taum )) * rkmv(i,j)

! T_taup = (1-beta_mld)T_taup + beta_mld (Tbar_taup - Tbar_taun)

        hstart = 0
        do irepeat = 1,k-1
        hstart = hstart + h(i,j,irepeat,taun)
        enddo
        hend = h(i,j,k,taun) + hstart
! upadting MYM 
        do kk = 1+hstart/5, hend/5
        if (kk .gt. kmax-1) goto 101
        
        do nt = 1,nn+2

        im = i-1
        ip = i+1
        jm = j-1
        jp = j+1


        rtemp1 = sum_2ptr_mld (1,ip,j,i,j,kk,nt) /
     &               sum_2ptr_mld (11,ip,j,i,j,kk,nt)

        rtemp2 = sum_2ptr_mld (1,im,j,i,j,kk,nt) /
     &               sum_2ptr_mld (11,im,j,i,j,kk,nt)

        tracedubdxmld(kk,nt) = (rtemp1 - rtemp2)*
     &         ( sum_2pt_mld (1,ip,j,i,j,k)/
     &          sum_2pt_mld (11,ip,j,i,j,k) ) * rdxh(i,j)

        rtemp3 = sum_2ptr_mld (1,i,j,i,jp,kk,nt) /
     &         sum_2ptr_mld (11,i,j,i,jp,kk,nt)

        rtemp4 = sum_2ptr_mld (1,i,j,i,jm,kk,nt) /
     &         sum_2ptr_mld (11,i,j,i,jm,kk,nt)

        tracedvbdymld(kk,nt) =(rtemp3-rtemp4 )* 
     &        ( sum_2pt_mld (2,i,jp,i,j,k)/
     &        sum_2pt_mld (12,i,jp,i,j,k) ) * rdyh(i,j)

        enddo

        call tsuv_diffuse (i,j,kk,diffuse_temp,diffuse_salt)

#include "cparam.h"

         temp(i,j,kk,2) = temp(i,j,kk,2) 
     &   - (tracedubdxmld(kk,1)
     &   + tracedvbdymld(kk,1)) * dt* beta_mldt * rkmh(i,j)
     &   + diffuse_temp * diffuse_MY * dt * rkmh(i,j)
     &   -  (temp(i,j,kk,2) -  temp_read(i,j,kk,month)) * dt * sec2day * 0.0

         salt(i,j,kk,2) = salt(i,j,kk,2) 
     &   - (tracedubdxmld(kk,2)
     &   + tracedvbdymld(kk,2)) *dt* beta_mldt * rkmh(i,j)
     &   + diffuse_salt * diffuse_MY * dt * rkmh(i,j)
     &   -  (salt(i,j,kk,2) -  salt_read(i,j,kk,month)) * dt * sec2day * 0.0

         uvel(i,j,kk,2) = uvel(i,j,kk,2) 
     &   - (tracedubdxmld(kk,3)
     &   + tracedvbdymld(kk,3)) *dt* beta_mldm  * rkmh(i,j)

         vvel(i,j,kk,2) = vvel(i,j,kk,2)
     &   - (tracedubdxmld(kk,4)
     &   + tracedvbdymld(kk,4)) *dt* beta_mldm  * rkmh(i,j)


        enddo

 101    continue

        enddo
        enddo
        enddo


        return
        end subroutine couple_rgmld 


         function sum_2ptr_mld (index, i1,j1, i2,j2,kk,nt)
                use size_mod, only: temp, salt, uvel, vvel, rkmh, rkmv, rkmu 
        implicit none

        real :: sum_2ptr_mld
        integer, intent (in) :: index, i1,j1, i2,j2,kk,nt

        real :: rtemp1, rtemp2

        if (index .eq. 1) then

         if (nt .eq. 1) rtemp1 = temp(i1,j1,kk,1)
         if (nt .eq. 1) rtemp2 = temp(i2,j2,kk,1)

         if (nt .eq. 2) rtemp1 = salt(i1,j1,kk,1)
         if (nt .eq. 2) rtemp2 = salt(i2,j2,kk,1)

         if (nt .eq. 3) rtemp1 = uvel(i1,j1,kk,1)
         if (nt .eq. 3) rtemp2 = uvel(i2,j2,kk,1)

         if (nt .eq. 4) rtemp1 = vvel(i1,j1,kk,1)
         if (nt .eq. 4) rtemp2 = vvel(i2,j2,kk,1)

        sum_2ptr_mld = (rtemp1 + rtemp2)

        endif

        if (index .eq. 11 .and. nt .eq. 1) then
        sum_2ptr_mld = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif
        if (index .eq. 11 .and. nt .eq. 2) then
        sum_2ptr_mld = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif
        if (index .eq. 11 .and. nt .eq. 3) then
        sum_2ptr_mld = (rkmu(i1,j1) + rkmu(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif
        if (index .eq. 11 .and. nt .eq. 4) then
        sum_2ptr_mld = (rkmv(i1,j1) + rkmv(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif

        end function sum_2ptr_mld


	function sum_2pt_mld (index, i1,j1, i2,j2, kk)
           use size_mod, only: u, v, taun, rkmu, rkmv
        implicit none
        real :: sum_2pt_mld
        integer, intent(in) :: index, i1,j1, i2,j2, kk
	
	if (index .eq. 1) then
	sum_2pt_mld = (u(i1,j1,kk,taun) + u(i2,j2,kk,taun))*1.0
	endif

	if (index .eq. 2) then
	sum_2pt_mld = (v(i1,j1,kk,taun) + v(i2,j2,kk,taun))*1.0
	endif

	if (index .eq. 11) then
	sum_2pt_mld = (rkmu(i1,j1) + rkmu(i2,j2))
	sum_2pt_mld = max(1.0,sum_2pt_mld)
	endif

	if (index .eq. 12) then
	sum_2pt_mld = (rkmv(i1,j1) + rkmv(i2,j2))
	sum_2pt_mld = max(1.0,sum_2pt_mld)
	endif

	end function sum_2pt_mld

        subroutine tsuv_diffuse (ii,jj,kk, diffuse_temp, diffuse_salt)
            use size_mod, only : temp, rdxu, rkmh, rdy, rdyv, salt, rdxh, rdyh
        implicit none

        integer :: ii, jj, kk
        real :: diffuse_temp, diffuse_salt, diiffuse_temp, diiffuse_salt

        integer :: iip, iim, jjp, jjm
        real :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6

        iip = ii+1
        iim = ii-1
        jjp = jj+1
        jjm = jj-1

        rtemp1 = rkmh(iip,jj)*
     &    (temp(iip,jj,kk,1) - temp(ii,jj,kk,1) )*rdxu(iip,jj)
        rtemp2 = rkmh(iim,jj)*
     &    (temp(ii,jj,kk,1) - temp(iim,jj,kk,1) )*rdxu(ii,jj)
        rtemp3 = (rtemp1 - rtemp1)*rdxh(ii,jj)

        rtemp4 = rkmh(ii,jjp)* (temp(ii,jjp,kk,1) - temp(ii,jj,kk,1))*rdyv(ii,jjp)
        rtemp5 = rkmh(ii,jjm)* (temp(ii,jj,kk,1) - temp(ii,jjm,kk,1))*rdyv(ii,jj)
        rtemp6 = (rtemp4 - rtemp5)*rdyh(ii,jj)  !


        diffuse_temp = rtemp3 + rtemp6
        if (rtemp3 .eq. 0.0 .or. rtemp6 .eq. 0.0) diffuse_temp = 0.0 !Prajeesh 


        rtemp1 = rkmh(iip,jj)*
     &    (salt(iip,jj,kk,1) - salt(ii,jj,kk,1) )*rdxu(iip,jj)
        rtemp2 = rkmh(iim,jj)*
     &    (salt(ii,jj,kk,1) - salt(iim,jj,kk,1) )*rdxu(ii,jj)
        rtemp3 = (rtemp1 - rtemp1)*rdxh(ii,jj)

        rtemp4 = rkmh(ii,jjp)*
     &   (salt(ii,jjp,kk,1) - salt(ii,jj,kk,1))*rdyv(ii,jjp)
        rtemp5 = rkmh(ii,jjm)*
     &   (salt(ii,jj,kk,1) - salt(ii,jjm,kk,1))*rdyv(ii,jj)
        rtemp6 = (rtemp4 - rtemp5)*rdyh(ii,jj)

        diffuse_salt = rtemp3 + rtemp6
        if (rtemp3 .eq. 0.0 .or. rtemp6 .eq. 0.0) diffuse_salt = 0.0 !Prajeesh 

        return
        end subroutine tsuv_diffuse

        end module couple_mod
