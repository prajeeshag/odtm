	program main
        use size_mod
        use param_mod
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c       main program for the 1&1/2 layer redu!ced gravity model
!c
!c
!c
!c
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!#include "size.h"
!#include "polar.h" ! polar coordinates for model
!#include "param.h"


        print *, 'myname is :', myname

         call initial_declaration

	 call initial_condition
	
	 call polar_coord


!c
!c Initial time-index values
!!c
	taum = 1
        taun = 2
        taup = 3
        taus = 4
	
!c
!c       do the integration
!c
	 
	loop_ind=0
	days = 3 !365*30
	month_start = 1
	loop_start = 1
	number_of_snap = 3 !12*3*30
#ifdef output_average
        call set_averages_to_zero
#endif
	month = month_start
	month_start_snap = month_start
	lpm = dpm(month)*day2sec/dt
	lpd = day2sec/dt
	month_wind = month_start
	iday_month = month_start
	iday_wind = 1
	iday_start = iday_month !c-1
	
	 
	 call stamp
	iday_start_snap = iday_start !- 7.0
	
	 call read_wind (month_wind)! why reading twice ? later on at line no 105
	 call read_forcing (month_wind)!may be - default wind, if no option
!cprescribed in makefile
#ifdef monthly_wind
	 call read_wind (month_wind)
	 call read_forcing (month_wind)
#endif
	
	 call check
	
#ifdef trace
	tracer_switch = 1
#else
	tracer_switch = 0
#endif
	loop_total = int(days*day2sec/dt)
	
!cccccccccccc timer.F cccccccc
	 call system_clock(itimer2,itimerrate,itimermax) 
	 time_switch = 1
!cccccccccccc timer.F cccccccc
	
	
	do loop = loop_start, (loop_total+loop_start)
	age_time = age_time + 1
!c	if (loop .eq. loop_start) call maph
	
	 loop_day = loop*dt/day2sec
#ifdef monthly_wind
	 if ( loop .gt. lpm) then 
	  month = month + 1
	  month_wind = month_wind + 1
	  lpm = lpm + dpm(month)*day2sec/dt
#ifdef monthly_climatology
	  if ( month .eq. 13) then 
	  month = 1
	  month_wind = 1
	  endif
#endif
	  if ( month .eq. 13) month = 1
	  call read_wind (month_wind)
	  call read_forcing (month_wind)
#ifdef  thermodynamic_forcing
	 call read_density (month_wind)
#endif
	
	 endif
#endif
	
#ifdef daily_wind
	 if ( loop .gt. lpd) then 
	  day = day + 1
	  iday_wind = iday_wind + 1
	  lpd = lpd + day2sec/dt
	  call read_wind (iday_wind)
	  call read_forcing (iday_wind)
	 endif
	 month_wind = iday_wind

#ifdef prescribe_S_boundary
        call read_boundary (iday_wind)
#endif
#endif
	
	 sum_adv=0.0
	
	 call interp_linear ! linear interpolation of wind forcing to model grid
        
	 call interp_forcing
	
#if defined smagorinsky_laplacian
	 call smagorinsky_coeff
	
	 call smagorinsky
#endif
	
#ifdef entrain
	 call entrain_detrain
#endif
#ifdef thermodynamic_forcing
	call average_density
#endif

	do i=1,imt-1
	do k=1,km-1
	do j=1,jmt-1
	 if (rkmh(i,j) .ne. 0.0) then
	 nmid = (jmt/2)+1
	 rnmid = (j-nmid)*dyd+0.25 + reflat
	 call clinic
	 endif
	enddo
	enddo
	enddo

        day_night = cos(loop*(2*3.14/(day2sec/dt))) + 1.0
        
#ifdef open_NS
	 call openb
#endif
#ifdef open_EW
	 call openb
#endif
	
#ifdef prescribeflow
	do k=1,km
	deltax(k) = 0.0
	rsumu(k) = 0.0
	do j=1,jmt
	rsumu(k) = rsumu(k) + u(1,k,j,taun)*h(1,k,j,taun)*dy
	enddo
	enddo
	do k=1,km
	rsumv(k) = 0.0
	rsumh(k) = 0.0
	rsumx = 0.0
	rsumy = 0.0
	do i=1,imt
	rsumv(k) = rsumv(k) + v(i,k,1,taun)*h(i,k,1,taun)*dx
	rsumh(k) = rsumh(k) + h(i,k,1,taun)*dx
	rsumx  = rsumx + rkmt(i,1)
	enddo
	enddo
	do j=1,jmt
	rsumy  = rsumy + rkmt(1,j)
	enddo
	do k=1,km
	deltax(k) = (rsumu(k) + rsumv(k))/rsumh(k)/rsumx
	enddo
	if ( mod (loop,loop_total/365./100.) .eq. 0) then
	write(*,*) rsumu(1) , rsumv(1) , deltax(1)
	endif
#endif

c        pvort(:,:,:) = 0.0
	do i=2,imt-1
	do k=1,km-1
	do j=2,jmt-1
	 nmid = (jmt/2)+1
	 rnmid = (j-nmid)*dyd+0.25 + reflat

	
#ifdef atmosphere
	 call atmos
#endif
	
	
        
         call stability_check (loop)

	 call pressure_integral (loop_ind)
	
	 rlct = rkmu(i,j) + rkmv(i,j)
	 if (rlct .ne. 0.0) then
	 call momentum (loop_ind)
	 endif
	

#ifdef trace
	 call tracer 
#endif
#ifdef age_tracer
	 call age
#endif
	
#ifdef density 
	 call layer_density
#endif

 1020	continue
	enddo
	enddo
	enddo


        call mixed_layer_physics

        call balance_pme
        if ( mod(loop,1) .eq. 0) then
         call couple_rgmld
        endif

#ifdef open_NS
	 call openb
#endif
#ifdef open_EW
	 call openb
#endif
	
#ifdef particle_trajectory
	call ptraj
#endif
	
	
        
	 call filter
	
	
	
#ifdef inversion
	 if ( mod(loop_day+1,1) .eq. 0) then
	 call inverse_model
	 endif
#endif
	

!c
!cc interchange time-index for leap-frog scheme.
!c 

        write(*,*) temp(imt/2, 1, jmt/2 + 20, 1), h(imt/2, 1, jmt/2 + 20, taun),
     & loop, SHCoeff(imt/2, 5, jmt/2 + 20)
        write(*,*) temp(imt/2, 1, jmt/2 + 20, 1), h(imt/2, 1, jmt/2 + 20, taun),
     & loop, SHCoeff(imt/2, 1, jmt/2 + 20)

#ifdef output_average
        call output_manager (number_of_snap)
#endif


	if ( mod (loop,loop_total/number_of_snap) .eq. 0) then
	 call restart
	endif
	 


!cccccccccccccccccccccccccccccccccccccccccccc
!c rotate the timestep once to achieve      c
!c leap-frog time difference.               c
!cccccccccccccccccccccccccccccccccccccccccccc
!c                 taum                     c
!c                v   ^                     c
!c               v o o ^                    c
!c              v ( | ) ^                   c
!c             v (  ~  ) ^                  c
!c            v   _ . _   ^                 c
!c         taup > > > > > taun              c
!c	      ktaum= taum                  c
!c	      taum = taun                  c
!c	      taun = taup                  !c
!cc	      taup = ktaum                 !c
!ccccccccccccccccccccccccccccccccccccccccccccc
	 do i=1,imt
	 do j=1,jmt
	 do k=1,km
	if ( u(i, k, j,taun) .ne. u(i, k, j,taun)
     &    .or. u(i, k, j,taun) .lt. -10.0 .or.
     &         u(i, k, j,taun) .gt. 10.0 ) then
	 loop_ind = loop_ind + 1
	 call write_output (loop_ind)
         stop 'stop=>blow-up'
	 endif
	 enddo
	 enddo
	 enddo
	
        
	if ( mod (loop,loop_total/number_of_snap) .eq. 0) then
	 loop_ind = loop_ind + 1
c	write(*,*)'Writing output for loop -> ', loop
	
	 call timer
	
	
	 call write_output (loop_ind)
#ifdef output_average
        call set_averages_to_zero
#endif
	endif
	
	
	enddo
	write (*,*)
	write (*,*)'Integration finished'

	
	end
