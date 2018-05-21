program main
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c       main program for the 1&1/2 layer redu!ced gravity model
    !c
    !c
    !c
    !c
    !c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    use size_mod, only : days, i, iday_start, iday_start_snap, itimer2
    use size_mod, only : itimermax, itimerrate, j, k, loop, loop_start, lpd
    use size_mod, only : lpm, month, month_start, month_start_snap, month_wind
    use size_mod, only : taum, taun, taup, taus, time_switch, tracer_switch 
    use size_mod, only : h, iday_wind, rkmh, rkmu, rkmv, shcoeff, temp, u
    use size_mod, only : imt, jmt, km, gdx, gdy, kmaxMYM, dz
    use size_mod, only : t

    use param_mod, only : day2sec, dpm, dt, dyd, loop_day, loop_ind, loop_total
    use param_mod, only : nmid, number_of_snap, reflat, rnmid, sum_adv, rdepth
    
    use momentum_mod, only : momentum
    use tracer_mod, only : tracer
    use couple_mod, only : couple_rgmld
    use presgrad_mod, only : pressure_integral
    
    use mpp_mod, only : mpp_npes, mpp_pe, mpp_error, stdout, FATAL, WARNING, NOTE, mpp_init
    use mpp_mod, only : mpp_exit, mpp_max, mpp_sum
    use mpp_io_mod, only : mpp_io_init, mpp_open, mpp_close, MPP_RDONLY, MPP_ASCII, MPP_MULTI
    use fms_mod,  only : field_exist, field_size, read_data, fms_init, fms_end
    use mpp_domains_mod, only : domain2d, domain1d, mpp_define_layout, mpp_define_domains
    use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_domain_components, mpp_update_domains
    use diag_manager_mod, only : diag_manager_init, register_diag_field, register_static_field
    use diag_manager_mod, only : diag_axis_init, send_data, diag_manager_end
    use diag_data_mod, only : FILL_VALUE
    use time_manager_mod, only : set_calendar_type, NO_CALENDAR, JULIAN, NOLEAP
    use time_manager_mod, only : time_type, set_time, set_date, operator(+), assignment(=)
    use time_manager_mod, only : print_time, set_ticks_per_second

    implicit none

    integer :: iday_month, ii
    real :: age_time, day_night, rlct, depth_mld
    type(time_type) :: time, time_step

    integer :: domain_layout(2), halo=1

    integer :: id_lon, id_lat, id_sst, id_depth_mld, id_depth, id_sss
    integer :: id_h, id_eta, id_u, id_v, id_tx, id_ty, id_temp, id_salt
    integer :: id_we, id_dens, id_pvort, id_mask, id_dxu, id_dyv
    integer :: id_temp_mld, id_salt_mld, id_uvel, id_vvel, id_diag, id_sh, id_sm
    integer :: id_mld, id_tke, id_rif, id_mlen, id_st_h, id_st_m

    type(domain2d) :: domain

    call init_odtm()

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
    days = 1 !365*30
    month_start = 1
    loop_start = 1
    number_of_snap = 6 !12*3*30
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

        time = time + time_step

        age_time = age_time + 1
        !c  if (loop .eq. loop_start) call maph
    
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
        write(*,*) temp(imt/2, 1, jmt/2 + 20, 1), h(imt/2, 1, jmt/2 + 20, taun), &
                   loop, SHCoeff(imt/2, 5, jmt/2 + 20)
        write(*,*) temp(imt/2, 1, jmt/2 + 20, 1), h(imt/2, 1, jmt/2 + 20, taun), &
                    loop, SHCoeff(imt/2, 1, jmt/2 + 20)

#ifdef output_average
        call output_manager (number_of_snap)
#endif
        
        call send_data_diag(time)


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
        !c        ktaum= taum                  c
        !c        taum = taun                  c
        !c        taun = taup                  !c
        !cc       taup = ktaum                 !c
        !ccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,imt
            do j=1,jmt
                do k=1,km
                    if ( u(i, k, j,taun) .ne. u(i, k, j,taun) &
                        .or. u(i, k, j,taun) .lt. -10.0 .or. &
                        u(i, k, j,taun) .gt. 10.0 ) then

                        loop_ind = loop_ind + 1
                        call write_output (loop_ind)
                        stop 'stop=>blow-up'
                    endif
                enddo
            enddo
        enddo
    
        
        if ( mod (loop,loop_total/number_of_snap) .eq. 0) then
            loop_ind = loop_ind + 1
            call timer
            call write_output (loop_ind)
#ifdef output_average
            call set_averages_to_zero
#endif
        endif
    enddo
    
    write (*,*)
    write (*,*)'Integration finished'
    
    time = time + time_step
    call diag_manager_end(time)
    call fms_end()
    
    contains


    subroutine init_odtm()

         integer :: ii
 
        call mpp_init()
        call fms_init()
        call set_calendar_type(NOLEAP)
        call diag_manager_init()


        call mpp_define_layout((/1,imt,1,jmt/),mpp_npes(),domain_layout)

        call mpp_define_domains((/1,imt,1,jmt/), domain_layout, domain, xhalo=halo, yhalo=halo )
        
        call initial_declaration

        call initial_condition

        call polar_coord
    
        time_step = set_time(seconds=int(dt))
        time = set_date(1995, 1, 1, 0, 0, 0)

        id_lon = diag_axis_init('lon', gdx(1:imt), 'degrees_east', cart_name='X', &
            long_name='longitude', domain2=domain)

        id_lat = diag_axis_init('lat', gdy(1:jmt), 'degrees_north', cart_name='Y', &
            long_name='latitude', domain2=domain) 

        id_depth_mld = diag_axis_init('depth_mld', (/(real(ii)*0.5,ii=1,kmaxMYM)/), 'meters', &
            cart_name='Z', long_name='depth')

        rdepth(1) = dz(1)

        do ii=2,km
            rdepth(ii)=rdepth(ii-1) + dz(ii)
        enddo

        id_depth = diag_axis_init('depth', rdepth, 'meters', &
            cart_name='Z', long_name='depth')

        id_sst = register_diag_field('odtm', 'sst', (/id_lon,id_lat/), init_time=Time, &
                 long_name='Sea Surface Temperature', units='deg-C',missing_value=FILL_VALUE)

        id_sss = register_diag_field('odtm', 'sss', (/id_lon,id_lat/), init_time=Time, &
                 long_name='Sea Surface Salinity', units='psu',missing_value=FILL_VALUE)

        id_temp = register_diag_field('odtm', 'temp', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='Temperature', units='deg-C',missing_value=FILL_VALUE)

        id_salt = register_diag_field('odtm', 'salt', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='Salinity', units='deg-C',missing_value=FILL_VALUE)

        id_h = register_diag_field('odtm', 'h', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='Hieght', units='meters',missing_value=FILL_VALUE)

        id_eta = register_diag_field('odtm', 'eta', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='eta', units='meters',missing_value=FILL_VALUE)

        id_u = register_diag_field('odtm', 'u', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='U-velocity', units='ms-1',missing_value=FILL_VALUE)

        id_v = register_diag_field('odtm', 'v', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='V-velocity', units='ms-1',missing_value=FILL_VALUE)

        id_tx = register_diag_field('odtm', 'tx', (/id_lon,id_lat/), init_time=Time, &
                 long_name='taux', units='?',missing_value=FILL_VALUE)

        id_ty = register_diag_field('odtm', 'ty', (/id_lon,id_lat/), init_time=Time, &
                 long_name='tauy', units='?',missing_value=FILL_VALUE)

        id_we = register_diag_field('odtm', 'we', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)
       
        id_dens = register_diag_field('odtm', 'dens', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='Density', units='?',missing_value=FILL_VALUE)

        id_pvort = register_diag_field('odtm', 'pvort', (/id_lon,id_lat,id_depth/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_mask = register_static_field('odtm', 'mask', (/id_lon,id_lat/), long_name='?', units='?', &
                    missing_value=FILL_VALUE )

        id_dxu = register_static_field('odtm', 'dxu', (/id_lon,id_lat/), long_name='?', units='?', &
                    missing_value=FILL_VALUE)

        id_dyv = register_static_field('odtm', 'dyv', (/id_lon,id_lat/), long_name='?', units='?', &
                    missing_value=FILL_VALUE)
         

        id_temp_mld = register_diag_field('odtm', 'temp_mld', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='Temperature', units='deg-C',missing_value=FILL_VALUE)
        
        id_salt_mld = register_diag_field('odtm', 'salt_mld', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='Salinity', units='psu',missing_value=FILL_VALUE)

        id_uvel = register_diag_field('odtm', 'uvel', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='U-velocity', units='ms-1',missing_value=FILL_VALUE)

        id_vvel = register_diag_field('odtm', 'vvel', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='V-velocity', units='ms-1',missing_value=FILL_VALUE)

        id_diag = register_diag_field('odtm', 'diag', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_sh = register_diag_field('odtm', 'sh', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_sm = register_diag_field('odtm', 'sm', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)
        
        id_mld = register_diag_field('odtm', 'mld', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_tke = register_diag_field('odtm', 'tke', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)
        
        id_rif = register_diag_field('odtm', 'rif', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_mlen = register_diag_field('odtm', 'mlen', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_st_h = register_diag_field('odtm', 'st_h', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_st_m = register_diag_field('odtm', 'st_m', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)


    end subroutine init_odtm


    subroutine send_data_diag(time)

        type(time_type) :: time
        integer :: used

        real :: tmp2(imt,jmt), tmp3(imt,jmt,1)

       
        if (id_sst > 0) then 
            tmp2 = FILL_VALUE
            where (rkmh>0.) tmp2 = t(1:imt,1,1:jmt,1,taun)
            used = send_data(id_sst, tmp2, time)
        endif

        if (id_sss > 0) then 
            tmp2 = FILL_VALUE
            where (rkmh>0.) tmp2 = t(1:imt,1,1:jmt,2,taun)
            used = send_data(id_sss, tmp2, time)
        endif

        if (id_temp>0) then
            do ii = 1, km
                tmp3 = FILL_VALUE
                where (rkmh>0.) tmp3(:,:,1) = t(:,ii,:,1,taun)
                used = send_data(id_temp,tmp3,time,ks_in=ii,ke_in=ii)
            enddo
        endif

        if (id_salt>0) then
            do ii = 1, km
                tmp3 = FILL_VALUE
                where (rkmh>0.) tmp3(:,:,1) = t(:,ii,:,2,taun)
                used = send_data(id_salt,tmp3,time,ks_in=ii,ke_in=ii)
            enddo
        endif

        if (id_h>0) then
            do ii = 1, km
                tmp3 = FILL_VALUE
                where (rkmh>0.) tmp3(:,:,1) = h(:,ii,:,taun)
                used = send_data(id_h,tmp3,time,ks_in=ii,ke_in=ii)
            enddo
        endif


    end subroutine send_data_diag

    
end program main
