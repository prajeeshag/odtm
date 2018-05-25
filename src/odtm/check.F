	subroutine check
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	
c	subroutin to display model parameters and exp. settings
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : k, dz, rho, he, days, imt, jmt, km, theta
        use param_mod, only : deg2rad, diffuse, diffuse_th, dt, dx, dx2, dxd
        use param_mod, only : dy, dy2, g, re
        implicit none
        
        real :: amin, rsumk, rtemp1, rtemp2, rmode_speed, rtemp3
        integer :: kk
#if !defined restrt && !defined monthly_climatology && !defined daily_wind
	write (*,*)
	write (*,*)' Error: one of the following option must be ',
     &  'enabled =>   restrt or monthly_climatology'
c	 stop  '==>check.F'
#endif

#if defined restrt && monthly_climatology
	 write (*,*)
	 write (*,*)
	 write(*,*)' Model is starting from a restart file and',
     &  ' using monthly climatology forcing'
	 write (*,*)
#endif

#if defined monthly_climatology && !defined restrt
	 write(*,*) 
	 write(*,*) 
	 write(*,*)' Model is starting with monthly climatology forcing'
	 write(*,*) 
#endif

#if defined restrt && !defined monthly_climatology
	 write(*,*) 
	 write(*,*) 
	 write(*,*)' Model is starting from a restart file'
	 write(*,*) 
#endif
#if defined monthly_wind && defined daily_wind
	 write(*,*)
	 write(*,*)
	 write(*,*)' only one option is possible:'  
	 write(*,*)' Dmonthly_wind or Ddaily_wind'
	 write(*,*)
#endif

#if defined trace && !defined trestrt
	 write(*,*)
	 write(*,*)
	 write(*,*)'Tracer is Initializing in a pre-defined domain'
	 write(*,*)'see initial.F for the details'
	 write(*,*)
#endif

#if defined trestrt && !defined trace
	 write(*,*)
	 write(*,*)'option -> trace <- is not enabled'
	 write(*,*)'option -> trestrt <- cannot operate without,
     & option -> trace <-'
	 write(*,*)
	 stop '=> check.F'
#endif

	
#ifdef trestrt
	open (11, file='restart.swm',form='unformatted',
     &   access='sequential')
	read (11) lon, lat, dxd, dyd, dt, diffuse, diffuse_th,
     &  diffuse_ah, tracer_switch
	close (11)
	if (tracer_switch .eq. 1) then
	 write(*,*)'Tracer is restarting from the restart file'
	else
	 write(*,*)'restart.swm  does not contain tracer values.'
	 write(*,*)'This restart file is only for option -> restrt <-.'
	 write(*,*)'This restart file cannot use with a combination of' 
	 write(*,*)'option -> restrt <- and option -> trestrt <-'
	 write(*,*)
	 write(*,*) lon, lat, dxd, dyd, dt, diffuse, diffuse_th,
     &  diffuse_ah, tracer_switch
	 stop '==> check.F'
	endif
#endif

#if defined restrt && !defined trace && !defined trestrt
	 open (11, file='restart.swm',form='unformatted',
     &   access='sequential')
        read (11) lon, lat, dxd, dyd, dt, diffuse, diffuse_th,
     &  diffuse_ah, tracer_switch
        close (11)
	 if (tracer_switch .eq. 1) then
	 write(*,*)'restart.swm contains tracer values. This restart'
	 write(*,*)'file can use only with options -> trestrt <-.'
	 write(*,*)
	 stop '=> check.F'
	endif
#endif


c
c checking viscocity limits ===== A = (dx)^3.beta
c
	 Amin = dx*dx*dx*2.e-11
	 if ( Amin .gt. diffuse) then
	 write (*,*) 'Amin is = ', Amin,' Diffuse = ', diffuse 
	 write (*,*) 'Viscocity is less than the required minimum value'
c	 stop
	write(*,*)
	write(*,*)
	write (*,*) 'WARNING!!! check.F found Viscocity coefficient'
	write (*,*) '"diffuse" is less than required minimum value'
	write (*,*) '(press contiunue anyway)'
c	read (*,*)
	 else
	 write (*,*) 'Amin is = ', Amin,' Diffuse = ', diffuse 
	 write (*,*) 'Viscocity is sufficient compared to minimum value'
	 endif
c
c check C.F.L limit
c
	rsumk = 0.0
	do k=1,1
	rsumk = rsumk + dz(k)
	enddo
	 dx2 = dxd*Re*cos(theta(1))*deg2rad
	 dx2 = dx2*dx2
	 rtemp1 = (1.0/dx2 + 1.0/dy2)
	 rtemp1 = 1.0/sqrt(rtemp1)
	 rtemp2 = ((rho(2)-rho(1))/rho(2))*g*rsumk
c	 rtemp2 = 0.40*1.0e4
	 rtemp2 = sqrt(rtemp2)
	 rtemp3 = rtemp1/(2*(rtemp2+2.0))  ! 2.0 is the |U|max
	 write (*,*)
	 write (*,19) 'C.F.L limit = ',rtemp3, ' , dt = ',dt
	 if ( rtemp3 .lt. dt ) then
	 write (*,*)  'C.F.L limit vialation. Time step (dt) is higher than required minimum value'
c	 stop '=>check.F'
	 endif
	
	rmode_speed = rtemp2
	
	 write (*,*)
	
	write (*,*)' Model parameters'
	write (*,*)' ----------------'
	write (*,*) 
	write (*,20)' layer thickness (h) = ', (dz(kk),kk=1,km-1)
	write (*,20)' baroclinc mod speed (m/s) = ', rmode_speed
	write (*,20)' minimum layer thickness (he) = ', (he(kk),kk=1,km-1)
	write (*,20)' reduced gravity = ', rtemp2/rsumk
	write (*,20)' laplacian viscosity = ', diffuse
	write (*,20)' laplacian mixing  = ', diffuse_th
	write (*,21)' time step dt = ', dt,' sec'
	write (*,26)' grid dimension (lon x lat) = ', imt,' x ',jmt
	write (*,21)' grid spacing (lon x lat) = ', dx,' x ',dy
	write (*,22)' total rundays = ', days,' = ', days/365.0,' years '
	write (*,24)'Running==='
	write (*,*)'Time ellapsed after each snapshot (minutes)'
	 
 19	format (a15,f8.0,a10,f8.0)
 20	format (a40,15f8.2)
 21	format (a40,f9.2,a4,f9.2)
 22	format (a40,f8.2,a4,f8.2, a6)
 23	format (a40,a20)
 24	format (a10,a20)
 25	format (a60)
 26	format (a40,i4,a4,i4)

	return
	stop
	end
