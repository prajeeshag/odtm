
        module momentum_mod
        use advec_mod, only : advec        
        use heat_flux_mod, only : heat_flux 
        implicit none

        private

        public :: momentum

        contains

        subroutine momentum ()
        !	subroutine to solve momentum equations

        use size_mod, only : i, j, k, loop, taux_force, temp, h, rho, tauy_force
        use size_mod, only : cld, uwnd, vwnd, airt, sphm, taun, dz, u, v, taum, taup
        use size_mod, only : rkmv, rEnergy, rkmu
        use size_mod, only : fcoruh, pres_gradu, fcorvh, pres_gradv
        use param_mod, only : advec_uu, advec_uv, advec_vu, advec_vv, day2sec
        use param_mod, only : diffuse_my, diffuse_tr, diffuse_uu, diffuse_vv, dtts
        use param_mod, only : g, polar_u, polar_v, dt

        use diffuse_two_mod, only: laplacian

!| added by Vinu on 27-06-2018 for entrainment momentum account
!v
        use size_mod, only : we, wd, km
        use mpp_domains_mod, only : domain2d, mpp_update_domains
        use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, mpp_sum
!^
!| added by Vinu on 27-06-2018 for entrainment momentum account

        implicit none
        type(domain2d) :: domain
        real tx,ty, Lv
        integer :: ip, im, jp, jm
        real :: diffuse_additional, rtemp1, cldin, uwndin, vwndin, speed, airt_in
        real :: sphm_in, sst_in, qnet_in, qs, ql, qb, h_ekman
        real :: rtemp_relax_u1, rtemp_relax_v1, relaxterm_u1, relaxterm_v1
        real :: relaxterm_u2, relaxterm_u3, relaxterm_v2, relaxterm_v3
        real :: rH1, rtemp2, rKE, rpert, rPE
        real :: taux_1, tauy_1
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1

!| added by Vinu on 27-06-2018 for entrainment momentum account
!v
        call mpp_update_domains(we(:,:,1:km-1),domain)
        call mpp_update_domains(wd(:,:,1:km-1),domain)
!^
!| added by Vinu on 27-06-2018 for entrainment momentum account

#ifdef thermodynamic_forcing
	do kkk=1,km
        rho(kkk) = dens(i,kkk,j)
	enddo
       if (dens(i,1,j) .eq. 0.0) rho(1) = 1025.25
       if (dens(i,2,j) .eq. 0.0) rho(2) = 1026.59
       if (dens(i,3,j) .eq. 0.0) rho(3) = 1027.03
       if (dens(i,4,j) .eq. 0.0) rho(4) = 1027.30
       if (dens(i,5,j) .eq. 0.0) rho(5) = 1027.53
       if (dens(i,6,j) .eq. 0.0) rho(6) = 1027.77
#endif
	
	
c	 diffuse_additional = 1.0*(cos(theta(j)/deg2rad))**0.3
	 diffuse_additional = 1.0
	
c
c artificial exagerration of friction at lower left corner
c 
	
c	 if ( i .lt. 11 .and. j .lt. 26 ) then
c	  diffuse_additional = 3.0
c	 endif
	
	call advec
	
#ifdef density
	call nonbosinique
#endif
	
	call corf
	
#ifdef constant_laplacian	
	call laplacian
#endif
	
	
#ifdef smagorinsky_biharmonic
	call biharmon
#endif
	
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1
	rtemp1 = taux_force(i,j)  * (1 - mod(loop,2))
	taux_1 = rtemp1 / rho(1)  
	rtemp1 = tauy_force(i,j)  * (1 - mod(loop,2))
	tauy_1 = rtemp1 / rho(1)

#ifdef trace
         cldin= cld(i,j) * 0.0
         uwndin = uwnd(i,j)
         vwndin = vwnd(i,j)
         speed = (uwndin**2 + vwndin**2)**0.5
         airt_in = airt(i,j)
         sphm_in = sphm(i,j)
         sst_in = temp(i,j,1,1) ! temperature of first layer
         call heat_flux(sphm_in,airt_in,sst_in,uwndin,vwndin,cldin,qnet_in,
     &   qs, ql, qb, tx, ty, Lv)
	taux_1 = (tx / rho(1)  ) * (1 - mod(loop,2)) 
	tauy_1 = (ty / rho(1)  ) * (1 - mod(loop,2)) 
#endif


c        goto 183
	
	 h_ekman = 05.0

	 if (k .eq. 1) then
	 if (h(i,j,1,taun) .lt. h_ekman) then
	 taux_1 = taux_1 * h(i,j,1,taun)/h_ekman
	 tauy_1 = tauy_1 * h(i,j,1,taun)/h_ekman
	 endif
	 endif
	 if (k .eq. 2) then
	 if (h(i,j,1,taun) .lt. h_ekman) then
	 taux_1 = taux_1 * (h_ekman - h(i,j,1,taun))/h_ekman
	 tauy_1 = tauy_1 * (h_ekman - h(i,j,1,taun))/h_ekman
	 endif
	 endif
	
	 if (k .eq. 2 .and. h(i,j,1,taun) .ge. h_ekman) then
	 taux_1 = 0.0
	 tauy_1 = 0.0
	 endif

c 183    continue

	if (k .gt. 2) taux_1 = 0.0
	if (k .gt. 2) tauy_1 = 0.0
	
	
	
!#ifdef restrt
!	if ( loop .eq. 1) dtts = dtts
!#else
!	if ( loop .eq. 1) dtts = dt
!#endif
	
c	 if (rkmh(i,j) .ne. 0.0) then
c	  write(*,*) diffuse_uu, diffuse_vv, loop
c	 endif
	
#ifdef density
	 in_key = 1
	 call rgravity (in_key)
#endif
	rtemp_relax_u1 = dz(k)
	rtemp_relax_v1 = dz(k)
	if (u(i,j,k,taum) .gt. 1.0 .and. u(i,j,k,taum) .ne. 0.0) 
     &   relaxterm_u1 = rtemp_relax_u1*(1.0 - u(i,j,k,taum))/864000.0/3.0
	if (v(i,j,k,taum) .gt. 1.0 .and. v(i,j,k,taum) .ne. 0.0) 
     &   relaxterm_v1 = rtemp_relax_v1*(1.0 - v(i,j,k,taum))/864000.0/3.0
	if (u(i,j,k,taum) .lt. 1.0 .and. u(i,j,k,taum) .ne. 0.0) 
     &   relaxterm_u1 = rtemp_relax_u1*(1.0 - u(i,j,k,taum))/864000.0/3.0
	if (v(i,j,k,taum) .lt. 1.0 .and. v(i,j,k,taum) .ne. 0.0) 
     &   relaxterm_v1 = rtemp_relax_v1*(1.0 - v(i,j,k,taum))/864000.0/3.0
	
	relaxterm_u1 = 0.0
	relaxterm_u2 = 0.0
	relaxterm_u3 = 0.0
	relaxterm_v1 = 0.0
	relaxterm_v2 = 0.0
	relaxterm_v3 = 0.0
	
	rH1 = sum_2ptm(3,i,j,im,j)/ sum_2ptm(31,i,j,im,j)
	rH1 = max(1.0, rH1)
	
	
        u(i,j,k,taup) =
     &   u(i,j,k,taum)*rH1
     &   + ( 
     &   + fcoruh
#if !defined density
     &    - pres_gradu(k)/rho(k)
#endif
#ifdef density
     &   - g*pres_gradu(k)
     &   - gdash_advec_x
     &   + dens_advec_x
#endif
     &    +taux_1 
#ifdef nonlinear_advec
     &    - advec_uu
     &    - advec_vu
     &    + polar_u
#endif
#ifdef entrain
     &    - 0.5*(wd(i,j,k)+wd(im,j,k))*u(i,j,k,taun)
     &    - 0.5*(we(i,j,k)+we(im,j,k))*u(i,j,k,taun)
#endif
#ifdef constant_laplacian
     &    + diffuse_uu
#endif
#ifdef smagorinsky_laplacian
     &    + diff_uu(i,k,j)
     &    + diff_vu(i,k,j)
     &    + polar_diff_u(i,k,j)
#endif
#ifdef smagorinsky_biharmonic
     &    - biharmon_u
#endif
     &    )*dtts*rkmu(i,j)

	rH1 = sum_2ptp(3,i,j,im,j)/ sum_2ptp(31,i,j,im,j)
	rH1 = max(1.0, rH1)
	u(i,j,k,taup) = u(i,j,k,taup)/rH1
	
	if (u(i,j,k,taup) .ne. u(i,j,k,taup)) then !Prajeesh isnan check
	 u(i,j,k,taup) = 0.0
	endif


        call drag_across_layers (1)
	
	
#ifdef density
	 in_key = 2
	 call rgravity (in_key)
#endif
	
	rH1 = sum_2ptm(3,i,j,i,jm)/ sum_2ptm(31,i,j,i,jm)
	rH1 = max(1.0, rH1)

         v(i,j,k,taup) =
     &   v(i,j,k,taum)*rH1
     &    + ( 
     &    - fcorvh 
#if !defined density
     &    - pres_gradv(k)/rho(k)
#endif
#ifdef density
     &   - gdash*pres_gradv(k)
     &   - gdash_advec_y
     &   + dens_advec_y
#endif
     &    + tauy_1 
#ifdef nonlinear_advec
     &    - advec_uv
     &    - advec_vv
     &    - polar_v
#endif
#ifdef entrain
     &    - 0.5*(wd(i,j,k)+wd(i,jm,k))*v(i,j,k,taun)
     &    - 0.5*(we(i,j,k)+we(i,jm,k))*v(i,j,k,taun)
#endif
#ifdef constant_laplacian
     &    + diffuse_vv
#endif
#ifdef smagorinsky_laplacian
     &    + diff_vv(i,k,j)
     &    + diff_uv(i,k,j)
     &    + polar_diff_v(i,k,j)
#endif
#ifdef smagorinsky_biharmonic
     &     - biharmon_v
#endif
     &    )*dtts*rkmv(i,j)

	rH1 = sum_2ptp(3,i,jm,i,j)/ sum_2ptp(31,i,jm,i,j)
	rH1 = max(1.0, rH1)
	v(i,j,k,taup) = v(i,j,k,taup)/rH1

	if (v(i,j,k,taup) .ne. v(i,j,k,taup)) then
	 v(i,j,k,taup) = 0.0
	endif
	

        call drag_across_layers (2)

	rtemp1 = 0.5*rho(k)
	rtemp2 = 0.5*rho(k)*((rho(k+1)-rho(k))/rho(k))*g
	rKE=rtemp1*dz(k)*((u(i,j,k,taun)**2)+(v(i,j,k,taun)**2))
	rpert = h(i,j,k,taun) - dz(k)
	rPE = rtemp2*(rpert**2)
	rEnergy(i,j,k) = rEnergy(i,j,k) + (rKE + rPE)*dtts
	relaxterm_u1 = 0.0
	relaxterm_u2 = 0.0
	relaxterm_u3 = 0.0
	relaxterm_v1 = 0.0
	relaxterm_v2 = 0.0
	relaxterm_v3 = 0.0

	
	return
	end subroutine momentum



	function sum_2ptp(index, i1,j1, i2,j2)
        use size_mod, only : k, u, v, h, taup, rkmu, rkmv, rkmh
        implicit none
        
        real :: sum_2ptp 
        integer, intent(in) :: index, i1, j1, i2, j2
        
        if (index .eq. 1) then
        sum_2ptp = (u(i1,j1,k,taup) + u(i2,j2,k,taup))
        endif

        if (index .eq. 2) then
        sum_2ptp = (v(i1,j1,k,taup) + v(i2,j2,k,taup))
        endif

        if (index .eq. 3) then
        sum_2ptp = (h(i1,j1,k,taup) + h(i2,j2,k,taup))
        endif

        if (index .eq. 11) then
        sum_2ptp = (rkmu(i1,j1) + rkmu(i2,j2))
        sum_2ptp = max(1.0,sum_2ptp)
        endif

        if (index .eq. 21) then
        sum_2ptp = (rkmv(i1,j1) + rkmv(i2,j2))
        sum_2ptp = max(1.0,sum_2ptp)
        endif

        if (index .eq. 31) then
        sum_2ptp = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptp = max(1.0,sum_2ptp)
        endif
        end function sum_2ptp

	function sum_2ptm(index, i1,j1, i2,j2)

        use size_mod, only: taum, k, u, v, h, rkmv, rkmu, rkmh
        implicit none

        real :: sum_2ptm
        integer, intent(in) :: index, i1, j1, i2, j2

        
        if (index .eq. 1) then
        sum_2ptm = (u(i1,j1,k,taum) + u(i2,j2,k,taum))
        endif

        if (index .eq. 2) then
        sum_2ptm = (v(i1,j1,k,taum) + v(i2,j2,k,taum))
        endif

        if (index .eq. 3) then
        sum_2ptm = (h(i1,j1,k,taum) + h(i2,j2,k,taum))
        endif

        if (index .eq. 11) then
        sum_2ptm = (rkmu(i1,j1) + rkmu(i2,j2))
        sum_2ptm = max(1.0,sum_2ptm)
        endif

        if (index .eq. 21) then
        sum_2ptm = (rkmv(i1,j1) + rkmv(i2,j2))
        sum_2ptm = max(1.0,sum_2ptm)
        endif

        if (index .eq. 31) then
        sum_2ptm = (rkmh(i1,j1) + rkmh(i2,j2))
        sum_2ptm = max(1.0,sum_2ptm)
        endif

        end function sum_2ptm 

        subroutine drag_across_layers (ind_drag)
            use size_mod, only : i, j, k, taup, rkmu, rkmv, loop, u, v, km
            use param_mod, only : day2sec, diffuse_my, diffuse_tr, dtts

        implicit none
            integer, intent(in) :: ind_drag
            real :: relax_on_off, alpha_rgm, alpha_rgt, beta_mldm
            real :: beta_mldt, gama, rdrag, relax_South
            real :: widS, relax_East, widE

#include "cparam.h"

        if ( k .gt. 1 .and. k .lt. km-1 .and. ind_drag .eq. 1) then
        u(i,j,k,taup) = u(i,j,k,taup)   
     & +rdrag * (u(i,j,k-1,taup) - u(i,j,k,taup)) * dtts * rkmu(i,j) 
     & +rdrag * (u(i,j,k+1,taup) - u(i,j,k,taup)) * dtts * rkmu(i,j)
        endif
        
        if ( k .gt. 1 .and. k .lt. km-1 .and. ind_drag .eq. 2) then
        v(i,j,k,taup) = v(i,j,k,taup)   
     & +rdrag * (v(i,j,k-1,taup) - v(i,j,k,taup)) * dtts * rkmv(i,j)
     & +rdrag * (v(i,j,k+1,taup) - v(i,j,k,taup)) * dtts * rkmv(i,j)
        endif

        return
        end subroutine drag_across_layers

        end module momentum_mod 
