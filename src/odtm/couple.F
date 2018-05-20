        module couple_mod
        implicit none

        private

        public :: couple_rgmld

        contains

        subroutine couple_rgmld 

        ! linear interpolation for wind forcing

        use size_mod, only : taum, taun, taup, taus, tracedubdxmld, tracedvbdymld
        use size_mod, only : loop, imt, jmt, km, i, j, k, rkmh, t, u, rkmu, v
        use size_mod, only : rkmv, h, nn, tmask, smask, umask, vmask, temp
        use size_mod, only : salt, uvel, vvel, rdx, dxu, rdy, dyu, dyv
        use size_mod, only : temp_read, salt_read, month 
        use param_mod, only : day2sec, dphi, dt, dthe, sec2day, theta

        implicit none

        double precision :: tempout_taum, tempout_taup
        double precision :: saltout_taum, saltout_taup
        double precision :: uvelout_taum, uvelout_taup
        double precision :: vvelout_taum, vvelout_taup

        double precision :: tempout_taumC, tempout_taupC
        double precision :: saltout_taumC, saltout_taupC
        double precision :: uvelout_taumC, uvelout_taupC
        double precision :: vvelout_taumC, vvelout_taupC

        double precision :: tempout_taumL, tempout_taupL
        double precision :: saltout_taumL, saltout_taupL
        double precision :: uvelout_taumL, uvelout_taupL
        double precision :: vvelout_taumL, vvelout_taupL

        double precision :: tempout_taumB, tempout_taupB
        double precision :: saltout_taumB, saltout_taupB
        double precision :: uvelout_taumB, uvelout_taupB
        double precision :: vvelout_taumB, vvelout_taupB

        double precision :: alpha_rgm, beta_mldm
        double precision :: alpha_rgt, beta_mldt

        real :: diffuse_temp, diffuse_salt, rtheta_vu, rtheta_vd
        real :: relax_on_off, gama, diffuse_MY, diffuse_tr, zero_out_tradv
        real :: relax_south, widS, relax_East, widE, kmax, rdrag
        real :: rtemp1, rtemp2, rtemp3, rtemp4
        integer :: hstart, hend, ii, jj, kk, tau_index
        integer :: im, jm, irepeat, nt, ip, jp

        tracedubdxmld(:,:) = 0.0
        tracedvbdymld(:,:) = 0.0
        
        taum = 1
        taun = 2
        taup = 3
        taus = 4

        alpha_rgm = 1.0  *      mod(loop,2)
        alpha_rgt = 1.0  * (1 - mod(loop,2))
        beta_mldm = 1.0  *      mod(loop,2)
        beta_mldt = 1.0  * (1 - mod(loop,2))

#include "cparam.h"  ! <== This file overrides above parameter setup

        kmax  = 51

        do i=2,imt-1
        do j=2,jmt-1
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


        t(i,k,j,1,taup) =  (t(i,k,j,1,taup) + 
     &   alpha_rgt*(tempout_taup - tempout_taum)) * rkmh(i,j)
        t(i,k,j,2,taup) =  (t(i,k,j,2,taup) + 
     &   alpha_rgt*(saltout_taup - saltout_taum)) * rkmh(i,j)
        u(i,k,j,taup) =  (u(i,k,j,taup) +
     &   alpha_rgm*(uvelout_taup - uvelout_taum )) * rkmu(i,j)
        v(i,k,j,taup) =  (v(i,k,j,taup) +
     &   alpha_rgm*(vvelout_taup - vvelout_taum )) * rkmv(i,j)



! T_taup = (1-beta_mld)T_taup + beta_mld (Tbar_taup - Tbar_taun)

        hstart = 0
        do irepeat = 1,k-1
        hstart = hstart + h(i,irepeat,j,taun)
        enddo
        hend = h(i,k,j,taun) + hstart
! upadting MYM 
        do kk = 1+hstart/5, hend/5
        if (kk .gt. kmax-1) goto 101
        
        do nt = 1,nn+2

        im = i-1
        ip = i+1
        jm = j-1
        jp = j+1


        do ii=im,ip
        do jj=jm,jp
        tmask(ii,jj) = 1.0
        smask(ii,jj) = 1.0
        umask(ii,jj) = 1.0
        vmask(ii,jj) = 1.0
        if(temp(ii,kk,jj,1).eq. 0.0)tmask(ii,jj) = 0.0
        if(salt(ii,kk,jj,1).eq. 0.0)smask(ii,jj) = 0.0
        if(uvel(ii,kk,jj,1).eq. 0.0)umask(ii,jj) = 0.0
        if(vvel(ii,kk,jj,1).eq. 0.0)vmask(ii,jj) = 0.0
        enddo
        enddo
        

        rtemp1 = sum_2ptr_mld (1,ip,j,i,j,kk,nt) /
     &               sum_2ptr_mld (11,ip,j,i,j,kk,nt)
        rtemp2 = sum_2ptr_mld (1,im,j,i,j,kk,nt) /
     &               sum_2ptr_mld (11,im,j,i,j,kk,nt)
        dphi = (rdx(i+1)+rdx(i))/2.0

        tracedubdxmld(kk,nt) = (rtemp1 - rtemp2)*
     &         (sum_2pt_mld (1,ip,j,i,j,k)/
     &          sum_2pt_mld (11,ip,j,i,j,k) ) /
     &         ( (dxu(i,j)+dxu(ip,j))/2  )


        rtheta_vu = cos((theta(j)+theta(j+1))*0.5)
        rtheta_vd = cos((theta(j)+theta(j-1))*0.5)
        dthe = (rdy(j) + rdy(j+1))/2.0

c        if (i .eq. 47 .and. j .eq. 28) then
c        write(*,*)   tracedvbdymld(kk,nt), nt, kk
c        endif

        rtemp3 = sum_2ptr_mld (1,i,j,i,jp,kk,nt) /
     &         sum_2ptr_mld (11,i,j,i,jp,kk,nt)
        rtemp4 = sum_2ptr_mld (1,i,j,i,jm,kk,nt) /
     &         sum_2ptr_mld (11,i,j,i,jm,kk,nt)
        tracedvbdymld(kk,nt) =(rtemp3-rtemp4 )* 
     &        (sum_2pt_mld (2,i,jp,i,j,k)/
     &        sum_2pt_mld (12,i,jp,i,j,k) ) /
     &        ( (dyv(i,jp)+dyv(i,j))/2 )

        enddo

c         rmld_misc(i,kk,j) = tracedubdxmld(kk,3)

        call tsuv_diffuse (i,j,kk,diffuse_temp,diffuse_salt)

#include "cparam.h"

         temp(i,kk,j,1) = temp(i,kk,j,2) - (tracedubdxmld(kk,1)
     &   + tracedvbdymld(kk,1)) * dt* beta_mldt * tmask(i,j)
     &   + diffuse_temp * diffuse_MY * dt * tmask(i,j)
     &   -  (temp(i,kk,j,2) -  temp_read(i,kk,j,month)) * dt * sec2day * 0.0
         salt(i,kk,j,1) = salt(i,kk,j,2) - (tracedubdxmld(kk,2)
     &   + tracedvbdymld(kk,2)) *dt* beta_mldt * smask(i,j)
     &   + diffuse_salt * diffuse_MY * dt * smask(i,j)
     &   -  (salt(i,kk,j,2) -  salt_read(i,kk,j,month)) * dt * sec2day * 0.0
         uvel(i,kk,j,1) = uvel(i,kk,j,2) - (tracedubdxmld(kk,3)
     &   + tracedvbdymld(kk,3)) *dt* beta_mldm  * umask(i,j)
         vvel(i,kk,j,1) = vvel(i,kk,j,2) - (tracedubdxmld(kk,4)
     &   + tracedvbdymld(kk,4)) *dt* beta_mldm  * vmask(i,j)


        enddo
 101    continue

        enddo
        enddo
        enddo


        return
        end subroutine couple_rgmld 


         function sum_2ptr_mld (index, i1,j1, i2,j2,kk,nt)
                use size_mod, only: temp, salt, uvel, vvel, tmask, umask, smask, vmask
        implicit none

        double precision :: sum_2ptr_mld
        integer, intent (in) :: index, i1,j1, i2,j2,kk,nt

        real :: rtemp1, rtemp2

        if (index .eq. 1) then

         if (nt .eq. 1) rtemp1 = temp(i1,kk,j1,1)
         if (nt .eq. 1) rtemp2 = temp(i2,kk,j2,1)

         if (nt .eq. 2) rtemp1 = salt(i1,kk,j1,1)
         if (nt .eq. 2) rtemp2 = salt(i2,kk,j2,1)

         if (nt .eq. 3) rtemp1 = uvel(i1,kk,j1,1)
         if (nt .eq. 3) rtemp2 = uvel(i2,kk,j2,1)

         if (nt .eq. 4) rtemp1 = vvel(i1,kk,j1,1)
         if (nt .eq. 4) rtemp2 = vvel(i2,kk,j2,1)

        sum_2ptr_mld = (rtemp1 + rtemp2)

c        if (rtemp1 .eq. 0.0 .and. rtemp2 .ne. 0.0) sum_2ptr_mld = 0.0
c        if (rtemp2 .eq. 0.0 .and. rtemp1 .ne. 0.0) sum_2ptr_mld = 0.0
        
        endif

        if (index .eq. 11 .and. nt .eq. 1) then
        sum_2ptr_mld = (tmask(i1,j1) + tmask(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif
        if (index .eq. 11 .and. nt .eq. 2) then
        sum_2ptr_mld = (smask(i1,j1) + smask(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif
        if (index .eq. 11 .and. nt .eq. 3) then
        sum_2ptr_mld = (umask(i1,j1) + umask(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif
        if (index .eq. 11 .and. nt .eq. 4) then
        sum_2ptr_mld = (vmask(i1,j1) + vmask(i2,j2))
        sum_2ptr_mld = max(1.0,sum_2ptr_mld)
        endif

        end function sum_2ptr_mld


	function sum_2pt_mld (index, i1,j1, i2,j2, kk)
           use size_mod, only: u, uvel, v, vvel, taun, rkmu, rkmv
        implicit none
        double precision :: sum_2pt_mld
        integer, intent(in) :: index, i1,j1, i2,j2, kk
	
	if (index .eq. 1) then
	sum_2pt_mld = (u(i1,kk,j1,taun) + u(i2,kk,j2,taun))*1.0 + 
     & 	(uvel(i1,kk,j1,taun) + uvel(i2,kk,j2,taun)) * 0.0
	endif

	if (index .eq. 2) then
	sum_2pt_mld = (v(i1,kk,j1,taun) + v(i2,kk,j2,taun))*1.0 +
     & 	(vvel(i1,kk,j1,taun) + vvel(i2,kk,j2,taun))* 0.0
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
            use size_mod, only : temp, dxu, rkmh, rdy, dyv, salt
            use param_mod, only : dthe
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
     &    (temp(iip,kk,jj,2) - temp(ii,kk,jj,2) )/dxu(iip,jj)
        rtemp2 = rkmh(iim,jj)*
     &    (temp(ii,kk,jj,2) - temp(iim,kk,jj,2) )/dxu(ii,jj)
        rtemp3 = (rtemp1 - rtemp2)/((dxu(iip,jj)+dxu(ii,jj))*0.5)

        rtemp4 = rkmh(ii,jjp)*
     &   (temp(ii,kk,jjp,2) - temp(ii,kk,jj,2))/dyv(ii,jjp)
        dthe = (rdy(jj-1) + rdy(jj))/2.0
        rtemp5 = rkmh(ii,jjm)*
     &   (temp(ii,kk,jj,2) - temp(ii,kk,jjm,2))/dyv(ii,jj)
        rtemp6 = (rtemp4 - rtemp5)/((dyv(ii,jjp)+dyv(ii,jj))*0.5)  !


        diffuse_temp = rtemp3 + rtemp6
        if (rtemp3 .eq. 0.0 .or. rtemp6 .eq. 0.0) diiffuse_temp = 0.0 !Prajeesh 


        rtemp1 = rkmh(iip,jj)*
     &    (salt(iip,kk,jj,2) - salt(ii,kk,jj,2) )/dxu(iip,jj)
        rtemp2 = rkmh(iim,jj)*
     &    (salt(ii,kk,jj,2) - salt(iim,kk,jj,2) )/dxu(ii,jj)
        rtemp3 = (rtemp1 - rtemp2)/((dxu(iip,jj)+dxu(ii,jj))*0.5)

        rtemp4 = rkmh(ii,jjp)*
     &   (salt(ii,kk,jjp,2) - salt(ii,kk,jj,2))/dyv(ii,jjp)
        dthe = (rdy(jj-1) + rdy(jj))/2.0
        rtemp5 = rkmh(ii,jjm)*
     &   (salt(ii,kk,jj,2) - salt(ii,kk,jjm,2))/dyv(ii,jj)
        rtemp6 = (rtemp4 - rtemp5)/((dyv(ii,jjp)+dyv(ii,jj))*0.5)  !

        diffuse_salt = rtemp3 + rtemp6
        if (rtemp3 .eq. 0.0 .or. rtemp6 .eq. 0.0) diiffuse_salt = 0.0 !Prajeesh 

        return
        end subroutine tsuv_diffuse

        end module couple_mod
