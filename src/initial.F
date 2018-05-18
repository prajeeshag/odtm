	subroutine initial_condition
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      initialze model with predefined velocity and surface elevation.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
        real tempin(201), saltin(201)
!#include "param.h"

	open (9, file='topo.swm')
	read (9,*)((rkmt(ii,jj),ii=1,imt),jj=1,jmt)
	close(9)
#ifdef density
	open  (1,file='dens.swm',form='unformatted',access='sequential')
	read  (1) ((dens_data(i,1,j),i=1,imt),j=1,jmt)
	close (1)
#endif
c	open (11, file='upwel_mask.swm', form='unformatted',
c     &access='sequential')
c	read (11) ((tr_source_term_2D(i,j),i=1,imt),j=1,jmt)
c	close (11)

	taum = 1
        taun = 2
        taup = 3
        taus = 4
! since a C grid is being used, u,v, h are defined at 3 different points, hence
! each has a dfferent mask.    u--------h
!                                       |
!                                       |       
!                                       |
!                                       v	
	rkmu(:,:) = 1.0
	rkmv(:,:) = 1.0
	rkmh(:,:) = 1.0
	do ii=1,imt
	do jj=1,jmt
	 if (rkmt(ii,jj) .eq. 0.0 .and. rkmt(ii,jj+1) .eq. 0.0) then
	 rkmu(ii,jj) = 0.0
	 endif
	enddo
	enddo
	do ii=1,imt-1
	do jj=1,jmt
	 if (rkmt(ii,jj) .eq. 0.0 .and. rkmt(ii+1,jj) .eq. 0.0) then
	 rkmv(ii,jj) = 0.0
	 endif
	enddo
	enddo
	do ii=1,imt-1
	do jj=1,jmt-1
	 if (rkmu(ii,jj) .eq. 0.0 .and. rkmu(ii+1,jj) .eq. 0.0 .and.
     &	     rkmv(ii,jj) .eq. 0.0 .and. rkmv(ii,jj+1) .eq. 0.0) then
	 rkmh(ii,jj) = 0.0
	 endif
	enddo
	enddo
	
c	open (11, file='wd_mask.swm')
c	read (11,*) ((wd_mask(i,j), i=1,imt), j=1,jmt)
c	close (11)
	
	rvalid_grid=0.0
	rvalid_grid_mass=0.0
	do ii=1,imt
	do jj=1,jmt
c	wd_mask (ii,jj) = wd_mask (ii,jj)*rkmt(ii,jj)
	we_upwel(ii,jj,1) = 0.0
	we_upwel(ii,jj,2) = 0.0
	
	if (rkmh(ii,jj) .ne. 0.0 ) then
	 rvalid_grid_mass = rvalid_grid_mass + 1
	endif
	
	 if (wd_mask(ii,jj)*rkmh(ii,jj) .eq. 1.0) 
     &           rvalid_grid = rvalid_grid+1
	do kk=1,km
	 h(ii,kk,jj,taum)=dz(kk)
	 h(ii,kk,jj,taun)=h(ii,kk,jj,taum) 
	 h(ii,kk,jj,taup)=h(ii,kk,jj,taum) 
	enddo
	do kk=1,km
	we(ii,kk,jj) = 0.0
	wd(ii,kk,jj) = 0.0
	do ll=1,4
	 u(ii,kk,jj,ll) = 0.0
	 v(ii,kk,jj,ll) = 0.0
	 eta(ii,kk,jj,ll) = 0.0
	enddo
	enddo
	enddo
	enddo
	rvalid_grid_background = rvalid_grid
	rvalid_grid_background_mass = rvalid_grid_mass
	write(*,*) rvalid_grid

	 do i=1,imt
	 do j=1,jmt
	 do k=1,km
	 do nt=1,nn
	  if ( i .gt. 170 .and. i .lt. 190 .and. j.gt.40 .and. j.lt.90) then
	  t(i,k,j,nt,taum) = 10.0
	  t(i,k,j,nt,taun) = 10.0
	 else 
	  t(i,k,j,nt,taum)=0.0
	  t(i,k,j,nt,taun)=0.0
	 endif
	 enddo
	 enddo
	 enddo
	 enddo
	
	do ii=1,imt
	do jj=1,jmt
	 if (rkmt(ii,jj) .eq. 0.0) mask (ii,jj) = 1.0
	 if (rkmt(ii,jj) .eq. 1.0) mask (ii,jj) = 0.0
	rrkmt(ii,jj) = rkmt(ii,jj)
c	rkmt(ii,jj) = h(ii,1,jj,taun)/h(ii,1,jj,taun)
	
	enddo
	enddo
	
c	do i=1,imt
c        gdx(i) = reflon + (i-1)*dxd
c        enddo
c        gdx(0) = reflon - dxd
c        gdx(imt+1) = gdx(imt) + dxd

c        do j=1,jmt
c        gdy(j) = reflat_start + (j-1)*dyd
c        enddo
c        gdy(0) = reflat_start - dyd
c        gdy(jmt+1) = gdy(jmt) + dyd
	
	open (11,file='xaxis_dat.dat')
	do i=1,imt
	read (11,*)gdx(i)
	enddo
	close (11)

	open (11,file='yaxis_dat.dat')
	do j=1,jmt
	read (11,*)gdy(j)
	enddo
	close (11)

	do i=1,imt-1
         rdx(i) = (gdx(i+1) - gdx(i))*deg2rad
        enddo
         rdx(0) = gdx(1)*deg2rad
         rdx(imt) = rdx(imt-1)*deg2rad
         rdx(imt+1) = rdx(imt)*deg2rad

        do j=1,jmt-1
         rdy(j) = (gdy(j+1) - gdy(j))*deg2rad
        enddo
         rdy(0) = gdy(1)*deg2rad
         rdy(jmt) = gdy(jmt-1)*deg2rad
         rdy(jmt+1) = gdy(jmt)*deg2rad
	
#ifdef trace
	open  (1,file='temp.swm',form='unformatted',access='sequential')
        do l=1,12
        read  (1) (((temp_read(i,k,j,l),i=1,imt),k=1,201),j=1,jmt)
        enddo
        close (1)
	open  (1,file='salt.swm',form='unformatted',access='sequential')
        do l=1,12
        read  (1) (((salt_read(i,k,j,l),i=1,imt),k=1,201),j=1,jmt)
        enddo
        close (1)
        do i=1,imt
        do j=1,jmt
c        open (1,file='temp.txt')
c        do kk=1,201
c         read(1,*) rtemp
c         temp_read(i,kk,j) = rtemp
c        enddo
c        close (1)
        do k=1,201
c         if ( k .le. 30) temp_read(i,k,j) = 25.0
c         if ( k .ge. 31 .and. k .le. 87) temp_read(i,k,j) = 20.0
c         if ( k .ge. 88 .and. k .le. 169) temp_read(i,k,j) = 15.0
c         if ( k .gt. 169 ) temp_read(i,k,j) = 5.0
c         temp_read(i,k,j) = 27.5 - 0.3*k

c         salt_read(i,k,j) = 35.0
         tempin(k) = temp_read(i,k,j,1)  
         saltin(k) = salt_read(i,k,j,1) 
        enddo
        do k=1,51
         temp(i,k,j,1) = temp_read(i,k,j,1)
         salt(i,k,j,1) = salt_read(i,k,j,1)
         temp(i,k,j,2) = temp_read(i,k,j,1)
         salt(i,k,j,2) = salt_read(i,k,j,1)
         uvel(i,k,j,1) = 0.0
         vvel(i,k,j,1) = 0.0
         uvel(i,k,j,2) = 0.0
         vvel(i,k,j,2) = 0.0
        enddo
        kmax = 201
        do k=1,km-1   
         call interp_extrap_initial (i,j,k,kmax,tempin,
     &   saltin,tempout,saltout)
        t(i,k,j,1,taun) = tempout
        t(i,k,j,2,taun) = saltout
        t(i,k,j,1,taum) = tempout
        t(i,k,j,2,taum) = saltout
        enddo
        t(i,km,j,1,taun) = 8.0
        t(i,km,j,2,taun) = 35
        t(i,km,j,1,taum) = 8.0
        t(i,km,j,2,taum) = 35.0
        enddo
        enddo

#endif
#ifdef restart_odtm
        call read_restart
#endif
	 do i=1,imt
	 do j=1,jmt
	 do k=1,km
	 do nt=1,nn
#if defined trace && !defined seasonal_tracer
c	 if (j .gt. yy0) then
c
c specified profile of tracers
c
c	  rexp1 = ((j-yy0)/7.0)**2
c	  t(i,k,j,nt,taum)= 0.0
c	  t(i,k,j,nt,taun)= 0.0
	age_psi(i,k,j,nt) = 0.0
        age_tpsi(i,k,j,nt) = 0.0
        age_t(i,k,j,nt) = 0.0
#else
	age_psi(i,k,j,nt) = 0.0
        age_tpsi(i,k,j,nt) = 0.0
        age_t(i,k,j,nt) = 0.0
	age_time = 0.0
c	 t(i,k,j,nt,taum)=0.0
c         t(i,k,j,nt,taun)=0.0
#endif
	 enddo
	 enddo
	 enddo
	 enddo
	
	
#ifdef particle_trajectory
	 x(1) = reflon
	 do i=2,imt
	 x(i) = x(1) + i*dxd
	 enddo
	 y(1) = reflat_start
	 do j=2,jmt
	 y(j) = y(1) + j*dyd
	 enddo
	 call initpt
#endif

#ifdef atmosphere
	do i=1,imt
	do j=1,jmt
	 exp1 = ((i-imt/1.5)/10.)**2 + ((j-jmt/2.)/20.)**2
	 u_atmos(i,j,taum) = 00.0*exp(-1.0*exp1)
	 v_atmos(i,j,taum) = 0.0
	 h_atmos(i,j,taum) = 1.0*exp(-1.0*exp1)
	 u_atmos(i,j,taun) = 0.0*exp(-1.0*exp1)
	 v_atmos(i,j,taun) = 0.0
	 h_atmos(i,j,taun) = 1.0*exp(-1.0*exp1)
	 q_atmos_heating (i,j) = 100.0*exp(-1.0*exp1)
	enddo
	enddo
	Q_atmos = 0.0
#endif
	
	
	return
	end
