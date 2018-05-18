	subroutine momentum (loop_ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	subroutine to solve momentum equations
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
        real tx,ty, Lv
!#include "param.h"
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1

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
         sst_in = temp(i,1,j,1) ! temperature of first layer
         call heat_flux(sphm_in,airt_in,sst_in,uwndin,vwndin,cldin,qnet_in,
     &   qs, ql, qb, tx, ty, Lv)
	taux_1 = (tx / rho(1)  ) * (1 - mod(loop,2))
	tauy_1 = (ty / rho(1)  ) * (1 - mod(loop,2)) 
#endif


c        goto 183
	
	 h_ekman = 20.0

	 if (k .eq. 1) then
	 if (h(i,1,j,taun) .lt. h_ekman) then
	 taux_1 = taux_1 * h(i,1,j,taun)/h_ekman
	 tauy_1 = tauy_1 * h(i,1,j,taun)/h_ekman
	 endif
	 endif
	 if (k .eq. 2) then
	 if (h(i,1,j,taun) .lt. h_ekman) then
	 taux_1 = taux_1 * (h_ekman - h(i,1,j,taun))/h_ekman
	 tauy_1 = tauy_1 * (h_ekman - h(i,1,j,taun))/h_ekman
	 endif
	 endif
	
	 if (k .eq. 2 .and. h(i,1,j,taun) .ge. h_ekman) then
	 taux_1 = 0.0
	 tauy_1 = 0.0
	 endif

c 183    continue

	if (k .gt. 2) taux_1 = 0.0
	if (k .gt. 2) tauy_1 = 0.0
	
	
	
#ifdef restrt
	if ( loop .eq. 1) dtts = dtts
#else
	if ( loop .eq. 1) dtts = dtts/2.0
#endif
	
c	 if (rkmh(i,j) .ne. 0.0) then
c	  write(*,*) diffuse_uu, diffuse_vv, loop
c	 endif
	
#ifdef density
	 in_key = 1
	 call rgravity (in_key)
#endif
	rtemp_relax_u1 = dz(k)
	rtemp_relax_v1 = dz(k)
	if (u(i,k,j,taum) .gt. 1.0 .and. u(i,k,j,taum) .ne. 0.0) 
     &   relaxterm_u1 = rtemp_relax_u1*(1.0 - u(i,k,j,taum))/864000.0/3.0
	if (v(i,k,j,taum) .gt. 1.0 .and. v(i,k,j,taum) .ne. 0.0) 
     &   relaxterm_v1 = rtemp_relax_v1*(1.0 - v(i,k,j,taum))/864000.0/3.0
	if (u(i,k,j,taum) .lt. 1.0 .and. u(i,k,j,taum) .ne. 0.0) 
     &   relaxterm_u1 = rtemp_relax_u1*(1.0 - u(i,k,j,taum))/864000.0/3.0
	if (v(i,k,j,taum) .lt. 1.0 .and. v(i,k,j,taum) .ne. 0.0) 
     &   relaxterm_v1 = rtemp_relax_v1*(1.0 - v(i,k,j,taum))/864000.0/3.0
	
	relaxterm_u1 = 0.0
	relaxterm_u2 = 0.0
	relaxterm_u3 = 0.0
	relaxterm_v1 = 0.0
	relaxterm_v2 = 0.0
	relaxterm_v3 = 0.0
	
	rH1 = sum_2ptm(3,i,j,im,j)/ sum_2ptm(31,i,j,im,j)
	rH1 = max(1.0, rH1)
	
	
        u(i,k,j,taup) =
     &   u(i,k,j,taum)*rH1
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
c     &    + 0.5*(wd(i,k,j)+wd(i-1,k,j))*u(i,k,j,taun)
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
	u(i,k,j,taup) = u(i,k,j,taup)/rH1
	
c	u(i,k,j,taup) = u(i,k,j,taup) *
c     & rkmt(i,j)*rkmt(ip,j)*rkmt(ip,jp)*rkmt(i,jp)
	
	if (u(i,k,j,taup) .ne. u(i,k,j,taup)) then
	 u(i,k,j,taup) = 0.0
	endif


        call drag_across_layers (1)
	
	
c	write(*,*) u(i,2,j,taun)
#ifdef density
	 in_key = 2
	 call rgravity (in_key)
#endif
	
	rH1 = sum_2ptm(3,i,j,i,jm)/ sum_2ptm(31,i,j,i,jm)
	rH1 = max(1.0, rH1)

         v(i,k,j,taup) =
     &   v(i,k,j,taum)*rH1
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
c     &    + 0.5*(wd(i,k,j)+wd(i,k,j-1))*v(i,k,j,taun)
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
	v(i,k,j,taup) = v(i,k,j,taup)/rH1

!	v(i,k,j,taup) = v(i,k,j,taup) *
!     & rkmt(i,j)*rkmt(ip,j)*rkmt(ip,jp)*rkmt(i,jp)

	if (v(i,k,j,taup) .ne. v(i,k,j,taup)) then
	 v(i,k,j,taup) = 0.0
	endif
	

        call drag_across_layers (2)

	rtemp1 = 0.5*rho(k)
	rtemp2 = 0.5*rho(k)*((rho(k+1)-rho(k))/rho(k))*g
	rKE=rtemp1*dz(k)*((u(i,k,j,taun)**2)+(v(i,k,j,taun)**2))
	rpert = h(i,k,j,taun) - dz(k)
	rPE = rtemp2*(rpert**2)
	rEnergy(i,k,j) = rEnergy(i,k,j) + (rKE + rPE)*dtts
	relaxterm_u1 = 0.0
	relaxterm_u2 = 0.0
	relaxterm_u3 = 0.0
	relaxterm_v1 = 0.0
	relaxterm_v2 = 0.0
	relaxterm_v3 = 0.0

	
	return
	end



	function sum_2ptp(index, i1,j1, i2,j2)
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

        
        if (index .eq. 1) then
        sum_2ptp = (u(i1,k,j1,taup) + u(i2,k,j2,taup))
        endif

        if (index .eq. 2) then
        sum_2ptp = (v(i1,k,j1,taup) + v(i2,k,j2,taup))
        endif

        if (index .eq. 3) then
        sum_2ptp = (h(i1,k,j1,taup) + h(i2,k,j2,taup))
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
        end

	function sum_2ptm(index, i1,j1, i2,j2)
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

        
        if (index .eq. 1) then
        sum_2ptm = (u(i1,k,j1,taum) + u(i2,k,j2,taum))
        endif

        if (index .eq. 2) then
        sum_2ptm = (v(i1,k,j1,taum) + v(i2,k,j2,taum))
        endif

        if (index .eq. 3) then
        sum_2ptm = (h(i1,k,j1,taum) + h(i2,k,j2,taum))
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

        end

        subroutine drag_across_layers (ind_drag)
            use size_mod
            use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"
        

#include "cparam.h"


        if ( k .gt. 1 .and. k .lt. km-1 .and. ind_drag .eq. 1) then
        u(i,k,j,taup) = u(i,k,j,taup)   
     & +rdrag * (u(i,k-1,j,taup) - u(i,k,j,taup)) * dtts * rkmu(i,j) 
     & +rdrag * (u(i,k+1,j,taup) - u(i,k,j,taup)) * dtts * rkmu(i,j)
        endif
        
        if ( k .gt. 1 .and. k .lt. km-1 .and. ind_drag .eq. 2) then
        v(i,k,j,taup) = v(i,k,j,taup)   
     & +rdrag * (v(i,k-1,j,taup) - v(i,k,j,taup)) * dtts * rkmv(i,j)
     & +rdrag * (v(i,k+1,j,taup) - v(i,k,j,taup)) * dtts * rkmv(i,j)
        endif

        return
        end
