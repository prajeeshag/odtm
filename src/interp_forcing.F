	subroutine interp_forcing
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	linear interpolation for wind forcing
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"


	 loop_year = day2sec*365.0/dt
	 do i=1,imt
	 do j=1,jmt
#ifdef monthly_wind
	rinterp_index = (0.5*(dpm(month)+dpm(month+1)))*day2sec/dt
#endif
#ifdef daily_wind
	rinterp_index = 1.0*day2sec/dt
#endif
	 x1 = sphm_read(i,j,1)! specific humidity
	 x2 = sphm_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 sphm(i,j) = x1 !*rkmu(i,j)
	 else 
	 sphm(i,j) = (sphm(i,j) + fracx ) !*rkmu(i,j)
	 endif
	
	 x1 = uwnd_read(i,j,1)
	 x2 = uwnd_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 uwnd(i,j) = x1 !*rkmu(i,j)
	 else 
	 uwnd(i,j) = (uwnd(i,j) + fracx ) !*rkmu(i,j)
	 endif
	
	 x1 = vwnd_read(i,j,1)
	 x2 = vwnd_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 vwnd(i,j) = x1 !*rkmu(i,j)
	 else 
	 vwnd(i,j) = (vwnd(i,j) + fracx ) !*rkmu(i,j)
         endif
        
	 x1 = airt_read(i,j,1)
	 x2 = airt_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 airt(i,j) = x1 !*rkmu(i,j)
	 else 
	 airt(i,j) = (airt(i,j) + fracx ) !*rkmu(i,j)
         endif
        
	 x1 = ssw_read(i,j,1) ! short wave flux
	 x2 = ssw_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 ssw(i,j) = x1 !*rkmu(i,j)
	 else 
	 ssw(i,j) = (ssw(i,j) + fracx ) !*rkmu(i,j)
         endif
        
	 x1 = cld_read(i,j,1)  ! cloud fraction
	 x2 = cld_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 cld(i,j) = x1 !*rkmu(i,j)
	 else 
	 cld(i,j) = (cld(i,j) + fracx ) !*rkmu(i,j)
         endif
        
	 x1 = pme_read(i,j,1)  ! outward heat flux
	 x2 = pme_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 pme(i,j) = x1 !*rkmu(i,j)
	 else 
	 pme(i,j) = (pme(i,j) + fracx ) !*rkmu(i,j)
         endif
        
	 x1 = chl_read(i,j,1)  ! outward heat flux
	 x2 = chl_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 chl(i,j) = x1 !*rkmu(i,j)
	 else 
	 chl(i,j) = (chl(i,j) + fracx ) !*rkmu(i,j)
         endif

	 x1 = rvr_read(i,j,1)  ! outward heat flux
	 x2 = rvr_read(i,j,2)
	 fracx = (x2-x1)/rinterp_index
	 if ( wind_switch .eq. 1) then
	 rvr(i,j) = x1 !*rkmu(i,j)
	 else 
	 rvr(i,j) = (rvr(i,j) + fracx ) !*rkmu(i,j)
         endif

	enddo
	enddo
        
	wind_switch = 0  
	return
	stop
	end
