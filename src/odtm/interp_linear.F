	subroutine interp_linear
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	linear interpolation for wind forcing
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : i, j, month, taux, taux_force, taux_snap, tauy
        use size_mod, only : tauy_force, tauy_snap, wind_switch, imt, jmt
        use param_mod, only : day2sec, dpm, dt, fracx, fracy
        implicit none
        
        integer :: loop_year
        real :: rinterp_index, x1, x2, y1, y2

	 loop_year = day2sec*365.0/dt
	 do i=1,imt
	 do j=1,jmt
#ifdef monthly_wind
	rinterp_index = (0.5*(dpm(month)+dpm(month+1)))*day2sec/dt
#endif
#ifdef daily_wind
	rinterp_index = 1.0*day2sec/dt
#endif
	 x1 = taux(i,j,1) ! why do we have third index?
	 y1 = tauy(i,j,1)
	
	 x2 = taux(i,j,2)
	 y2 = tauy(i,j,2)

	 fracx = (x2-x1)/rinterp_index
	 fracy = (y2-y1)/rinterp_index

	 if ( wind_switch .eq. 1) then ! wind_switch =1 -> donot interpolate the data to model grid,  =0->interpolate the data to model grid, (grid time steps). Need not interpolate at the first time because, 1st timestep acn be feeded as same.
	 taux_force(i,j) = x1 !*rkmu(i,j) *rkmv(i,j)
	 tauy_force(i,j) = y1 !*rkmv(i,j) *rkmu(i,j)!SHIKHA next two lines are repeat of this, why not assign than calculate? saves cost
	 taux_snap(i,j) = x1 !*rkmu(i,j) *rkmv(i,j)
	 tauy_snap(i,j) = y1 !*rkmv(i,j) *rkmu(i,j)
	 else 
	 taux_force(i,j) = (taux_force(i,j) + fracx ) !*rkmu(i,j)*rkmv(i,j)
         tauy_force(i,j) = (tauy_force(i,j) + fracy ) !*rkmv(i,j)*rkmu(i,j)
	 taux_snap(i,j) = (taux_snap(i,j) + fracx ) !*rkmu(i,j)*rkmv(i,j)
         tauy_snap(i,j) = (tauy_snap(i,j) + fracy ) !*rkmv(i,j)*rkmu(i,j)
	 endif
	
	enddo
	enddo

#ifdef thermodynamic_forcing
	 do i=1,imt
	 do j=1,jmt
	 do k=1,kmdata
#ifdef monthly_wind
	rinterp_index = (0.5*(dpm(month)+dpm(month+1)))*day2sec/dt
#endif
#ifdef daily_wind
	rinterp_index = 1.0*day2sec/dt
#endif
	 x1 = dens_in(i,k,j,1)
	 x2 = dens_in(i,k,j,2)

	 fracx = (x2-x1)/rinterp_index

	 if ( wind_switch .eq. 1) then
	 dens_data(i,k,j) = x1
	 else 
	 dens_data(i,k,j) = dens_data(i,k,j) + fracx
	 endif
	
	enddo
	enddo
	enddo
#endif
!	wind_switch = 0  ! THIS HAS TO BE ENABLED IF NO THREMO-FORCING
	return
	stop
	end
