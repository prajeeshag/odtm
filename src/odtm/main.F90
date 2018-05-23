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

    use size_mod, only : days, i, iday_start, itimer2
    use size_mod, only : itimermax, itimerrate, j, k, loop, loop_start, lpd
    use size_mod, only : lpm, month, month_start, month_wind
    use size_mod, only : taum, taun, taup, taus, time_switch, tracer_switch 
    use size_mod, only : iday_wind, rkmh, rkmu, rkmv, shcoeff, rkmt
    use size_mod, only : imt, jmt, km, gdx, gdy, kmaxMYM, dz, nn
    use size_mod, only : t, eta, u, v, temp, h, we, pvort, salt, dxu, dyv
    use size_mod, only : uvel, vvel, smcoeff, SHCoeff, diag_ext1, diag_ext2
    use size_mod, only : diag_ext3, diag_ext4, diag_ext5, diag_ext6
    use size_mod, only : sphm, uwnd, vwnd, airt, ssw, cld, pme, chl, rvr
    use size_mod, only : taux_force, tauy_force

    use param_mod, only : day2sec, dpm, dt, dyd, loop_day, loop_total
    use param_mod, only : nmid, reflat, rnmid, sum_adv
    use param_mod, only : denss, rmld_misc
    
    use momentum_mod, only : momentum
    use tracer_mod, only : tracer
    use couple_mod, only : couple_rgmld
    use presgrad_mod, only : pressure_integral
    use interp_extrap_initial_mod, only : interp_extrap_initial
    
    use mpp_mod, only : mpp_npes, mpp_pe, mpp_error, stdout, FATAL, WARNING, NOTE, mpp_init
    use mpp_mod, only : mpp_exit, mpp_max, mpp_sum
    use fms_mod,  only : field_exist, field_size, read_data, fms_init, fms_end
    use fms_io_mod, only : register_restart_field, restart_file_type, save_restart, restore_state
    use fms_io_mod, only : open_namelist_file, open_file, close_file, file_exist

    use mpp_domains_mod, only : domain2d, domain1d, mpp_define_layout, mpp_define_domains
    use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_domain_components, mpp_update_domains
    use diag_manager_mod, only : diag_manager_init, register_diag_field, register_static_field
    use diag_manager_mod, only : diag_axis_init, send_data, diag_manager_end
    use diag_data_mod, only : FILL_VALUE
    use data_override_mod, only : data_override_init, data_override
    use time_manager_mod, only : set_calendar_type, NO_CALENDAR, JULIAN, NOLEAP, date_to_string
    use time_manager_mod, only : time_type, set_time, set_date, operator(+), assignment(=)
    use time_manager_mod, only : print_time, set_ticks_per_second, increment_date, operator(>=)

    implicit none

    integer :: iday_month, ii
    real :: age_time, day_night, rlct, depth_mld
    type(time_type) :: time, time_step, time_restart

    integer :: domain_layout(2), halo=1, used

    integer :: id_lon, id_lat, id_sst, id_depth_mld, id_depth, id_sss, id_airt
    integer :: id_h, id_eta, id_u, id_v, id_tx, id_ty, id_temp, id_salt
    integer :: id_we, id_dens, id_pvort, id_mask, id_dxu, id_dyv
    integer :: id_temp_mld, id_salt_mld, id_u_mld, id_v_mld, id_diag, id_sh, id_sm
    integer :: id_mld, id_tke, id_rif, id_mlen, id_st_h, id_st_m, id_pme
    integer :: id_sphm, id_uwnd, id_vwnd, id_ssw, id_cld, id_chl, id_rvr

    type(domain2d) :: domain

    type(restart_file_type) :: restart_odtm

    real :: tmp2(imt,jmt), tmp3(imt,jmt,km), tmp3m(imt,jmt,kmaxMYM), rdepth(km)

    logical :: lmask(imt,jmt), lmask3(imt,jmt,km), lmask3m(imt,jmt,kmaxMYM)
    logical :: override
    character (len=32) :: timestamp
    integer :: restart_interval(6) = 0
    
    namelist /main_nml/ restart_interval

    call init_odtm()

    !c Initial time-index values
    taum = 1
    taun = 2
    taup = 3
    taus = 4

    !c       do the integration
    days = 365 !365*30
    month_start = 1
    loop_start = 1
    month = month_start
    lpm = dpm(month)*day2sec/dt
    lpd = day2sec/dt
    month_wind = month_start
    iday_month = month_start
    iday_wind = 1
    iday_start = iday_month !c-1
    
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
        endif
#endif
    
#ifdef daily_wind
        if ( loop .gt. lpd) then 
            day = day + 1
            iday_wind = iday_wind + 1
            lpd = lpd + day2sec/dt
        endif

        month_wind = iday_wind

#ifdef prescribe_S_boundary
        call read_boundary (iday_wind)
#endif
#endif
    
        sum_adv=0.0
    
        call data_override('OCN','sphm',sphm,time,override)
        if (.not.override) call mpp_error(WARNING, 'sphm not overriden')
        used = send_data(id_sphm, sphm, time)

        call data_override('OCN','uwnd',uwnd,time,override)
        if (.not.override) call mpp_error(WARNING, 'uwnd not overriden')
        used = send_data(id_uwnd, uwnd, time)

        call data_override('OCN','vwnd',vwnd,time,override)
        if (.not.override) call mpp_error(WARNING, 'vwnd not overriden')
        used = send_data(id_vwnd, vwnd, time)

        call data_override('OCN','airt',airt,time,override)
        if (.not.override) call mpp_error(WARNING, 'airt not overriden')
        used = send_data(id_airt, airt, time)

        call data_override('OCN','ssw',ssw,time,override)
        if (.not.override) call mpp_error(WARNING, 'ssw not overriden')
        used = send_data(id_ssw, ssw, time)

        call data_override('OCN','cld',cld,time,override)
        if (.not.override) call mpp_error(WARNING, 'cld not overriden')
        used = send_data(id_cld, cld, time)

        call data_override('OCN','pme',pme,time,override)
        if (.not.override) call mpp_error(WARNING, 'pme not overriden')
        used = send_data(id_pme, pme, time)

        call data_override('OCN','chl',chl,time,override)
        if (.not.override) call mpp_error(WARNING, 'chl not overriden')
        used = send_data(id_chl, chl, time)

        call data_override('OCN','rvr',rvr,time,override)
        if (.not.override) call mpp_error(WARNING, 'rvr not overriden')
        used = send_data(id_rvr, rvr, time)

        call data_override('OCN','taux_force',taux_force,time,override)

        call data_override('OCN','tauy_force',tauy_force,time,override)

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
                    call stability_check ()

                    call pressure_integral ()
    
                    rlct = rkmu(i,j) + rkmv(i,j)
                    if (rlct .ne. 0.0) then
                        call momentum ()
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

        call send_data_diag(time)


        time = time + time_step
    
        if ( time >= time_restart ) then
            
           timestamp = date_to_string(time)

           time_restart = increment_date(time, restart_interval(1), restart_interval(2), &
                restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )

           call save_restart(restart_odtm,timestamp) 

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
                        
                        call diag_manager_end(time)
                        stop 'stop=>blow-up'
        
                    endif
                enddo
            enddo
        enddo
    enddo
    
    write (*,*)
    write (*,*)'Integration finished'

    call save_restart(restart_odtm) 
    
    call diag_manager_end(time)
    call fms_end()
    
    contains


    subroutine init_odtm()

        integer :: ii, used, unit

        unit = open_namelist_file()
        rewind(unit)
        read(unit,nml=main_nml) 
         
        call mpp_init()
        call fms_init()
        call set_calendar_type(NOLEAP)

        time_step = set_time(seconds=int(dt))

        time = set_date(1995, 1, 1, 0, 0, 0)
      
        if (all(restart_interval==0)) restart_interval(1) = 10
 
        time_restart = increment_date(time, restart_interval(1), restart_interval(2), &
        restart_interval(3), restart_interval(4), restart_interval(5), restart_interval(6) )

        call diag_manager_init()

        call mpp_define_layout((/1,imt,1,jmt/),mpp_npes(),domain_layout)

        call mpp_define_domains((/1,imt,1,jmt/), domain_layout, domain, xhalo=halo, yhalo=halo )

        call data_override_init(Ocean_domain_in=domain)
        
        call initial_declaration

        call init_grid()

        !call initial_condition

        call polar_coord
    


        lmask=.false.; lmask3 = .false.; lmask3m = .false.

        lmask(2:imt-1,2:jmt-1)=rkmh(2:imt-1,2:jmt-1)>0

        do ii = 1,km
            lmask3(2:imt-1, 2:jmt-1, ii) = rkmh(2:imt-1,2:jmt-1)>0
        enddo

        do ii = 1,kmaxMYM
            lmask3m(2:imt-1, 2:jmt-1, ii) = rkmh(2:imt-1,2:jmt-1)>0
        enddo

        id_lon = diag_axis_init('lon', gdx(1:imt), 'degrees_east', cart_name='X', &
            long_name='longitude', domain2=domain)

        id_lat = diag_axis_init('lat', gdy(1:jmt), 'degrees_north', cart_name='Y', &
            long_name='latitude', domain2=domain) 

        id_depth_mld = diag_axis_init('depth_mld', (/(real(ii)*5.0,ii=1,kmaxMYM)/), 'meters', &
            cart_name='Z', long_name='depth')

        rdepth(1) = dz(1)

        do ii=2,km
            rdepth(ii)=rdepth(ii-1) + dz(ii)
        enddo

        id_depth = diag_axis_init('depth', rdepth, 'meters', &
            cart_name='Z', long_name='depth')

        id_airt = register_diag_field('odtm', 'airt', (/id_lon,id_lat/), init_time=Time, &
                 long_name='Air Temperature', units='deg-C',missing_value=FILL_VALUE)

        id_sphm = register_diag_field('odtm', 'sphm', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_uwnd = register_diag_field('odtm', 'uwnd', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_vwnd = register_diag_field('odtm', 'vwnd', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_ssw = register_diag_field('odtm', 'ssw', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_cld = register_diag_field('odtm', 'cld', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_pme = register_diag_field('odtm', 'pme', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_rvr = register_diag_field('odtm', 'rvr', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

        id_chl = register_diag_field('odtm', 'chl', (/id_lon,id_lat/), init_time=Time, &
                 long_name='?', units='?',missing_value=FILL_VALUE)

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

        id_u_mld = register_diag_field('odtm', 'u_mld', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
                 long_name='U-velocity', units='ms-1',missing_value=FILL_VALUE)

        id_v_mld = register_diag_field('odtm', 'v_mld', (/id_lon,id_lat,id_depth_mld/), init_time=Time, &
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

        tmp2 = 0.
        where(lmask) tmp2 = 1.

        used = send_data(id_mask, rkmt, time)
        used = send_data(id_dxu, dxu, time, mask=lmask)
        used = send_data(id_dyv, dyv, time, mask=lmask)

        unit = open_file(file='RESTART/._tmp_',action='write')
        call close_file(unit,'delete')

    end subroutine init_odtm


    subroutine send_data_diag(time)

        type(time_type) :: time
        integer :: used

       
        used = send_data(id_sst,t(1:imt,1,1:jmt,1,taun), time, mask=lmask)

        used = send_data(id_sss,t(1:imt,1,1:jmt,2,taun), time, mask=lmask)

        if (id_temp>0) then
            do ii = 1, km
                tmp3(:,:,ii) = t(:,ii,:,1,taun)
            enddo
            used = send_data(id_temp,tmp3,time, mask=lmask3)
        endif

        if (id_salt>0) then
            do ii = 1, km
                tmp3(:,:,ii) = t(:,ii,:,2,taun)
            enddo
            used = send_data(id_salt,tmp3,time, mask=lmask3)
        endif

        if (id_h>0) then
            do ii = 1, km
                tmp3(:,:,ii) = h(:,ii,:,taun)
            enddo
            used = send_data(id_h,tmp3,time, mask=lmask3)
        endif

        if (id_eta>0) then
            do ii = 1, km
                tmp3(:,:,ii) = eta(:,ii,:,1)
            enddo
            used = send_data(id_eta,tmp3,time, mask=lmask3)
        endif

        if (id_u>0) then
            do ii = 1, km
                tmp3(:,:,ii) = u(:,ii,:,taun)
            enddo
            used = send_data(id_u,tmp3,time, mask=lmask3)
        endif

        if (id_v>0) then
            do ii = 1, km
                tmp3(:,:,ii) = v(:,ii,:,taun)
            enddo
            used = send_data(id_v,tmp3,time, mask=lmask3)
        endif

        if (id_we>0) then
            do ii = 1, km
                tmp3(:,:,ii) = we(:,ii,:)
            enddo
            used = send_data(id_we,tmp3,time, mask=lmask3)
        endif

        if (id_dens>0) then
            do ii = 1, km
                tmp3(:,:,ii) = denss(:,ii,:)
            enddo
            used = send_data(id_dens,tmp3,time, mask=lmask3)
        endif

        if (id_pvort>0) then
            do ii = 1, km
                tmp3(:,:,ii) = pvort(:,ii,:)
            enddo
            used = send_data(id_pvort,tmp3,time, mask=lmask3)
        endif
        
        if (id_temp_mld>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = temp(:,ii,:,1)
            enddo
            used = send_data(id_temp_mld,tmp3m,time, mask=lmask3m)
        endif

        if (id_salt_mld>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = salt(:,ii,:,1)
            enddo
            used = send_data(id_salt_mld,tmp3m,time, mask=lmask3m)
        endif

        if (id_u_mld>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = uvel(:,ii,:,taun)
            enddo
            used = send_data(id_u_mld,tmp3m,time, mask=lmask3m)
        endif

        if (id_v_mld>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = vvel(:,ii,:,taun)
            enddo
            used = send_data(id_v_mld,tmp3m,time, mask=lmask3m)
        endif

        if (id_diag>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = rmld_misc(:,ii,:)
            enddo
            used = send_data(id_diag,tmp3m,time, mask=lmask3m)
        endif

        if (id_sh>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = shcoeff(:,ii,:)
            enddo
            used = send_data(id_sh,tmp3m,time, mask=lmask3m)
        endif

        if (id_sm>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = smcoeff(:,ii,:)
            enddo
            used = send_data(id_sm,tmp3m,time, mask=lmask3m)
        endif

        if (id_tke>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = diag_ext1(:,ii,:)
            enddo
            used = send_data(id_tke,tmp3m,time, mask=lmask3m)
        endif

        if (id_rif>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = diag_ext2(:,ii,:)
            enddo
            used = send_data(id_rif,tmp3m,time, mask=lmask3m)
        endif

        if (id_mlen>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = diag_ext3(:,ii,:)
            enddo
            used = send_data(id_mlen,tmp3m,time, mask=lmask3m)
        endif

        if (id_st_h>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = diag_ext4(:,ii,:)
            enddo
            used = send_data(id_st_h,tmp3m,time, mask=lmask3m)
        endif

        if (id_st_m>0) then
            do ii = 1, kmaxMYM
                tmp3m(:,:,ii) = diag_ext5(:,ii,:)
            enddo
            used = send_data(id_st_m,tmp3m,time, mask=lmask3m)
        endif

    end subroutine send_data_diag
    


    subroutine init_grid

        use size_mod, only : rdx, rdy, rkmt, we_upwel, wd, rrkmt
        use size_mod, only : temp_read, salt_read

        use param_mod, only : deg2rad, mask
    
        !c initialze model grid.
    
        implicit none
    
        real :: tempin(201), saltin(201)
        integer :: ii, jj, kk, kmax, l, ll, nt
        integer :: dimz(4), nlon, nlat
        real :: saltout, tempout
        real, allocatable :: tmp2(:,:)
        integer :: id_restart
        character (len=32) :: grid_file='INPUT/grid_spec.nc'
        character (len=32) :: temp_clim_file='INPUT/temp_clim.nc'
        character (len=32) :: salt_clim_file='INPUT/salt_clim.nc'
        character (len=32) :: restart_file='odtm_restart.nc'
  
        taum = 1
        taun = 2
        taup = 3
        taus = 4

        if (.not.field_exist(grid_file, "geolon_t"))  &
            call mpp_error(FATAL,'geolon_t not present in '//trim(grid_file))

        if (.not.field_exist(grid_file, "geolat_t"))  &
            call mpp_error(FATAL,'geolat_t not present in '//trim(grid_file))

        if (.not.field_exist(grid_file, "rkmt"))  &
            call mpp_error(FATAL,'rkmt not present in '//trim(grid_file))

       
        call field_size(grid_file, "geolon_t", dimz)
        
        nlon = dimz(1); nlat = dimz(2)
        
        allocate(tmp2(nlon,nlat))
        
        call read_data(grid_file, 'geolon_t', tmp2, no_domain=.true.)
   
        gdx = 0. ; gdy = 0.
        gdx(1:nlon) = tmp2(:,1)

        call read_data(grid_file, 'geolat_t', tmp2, no_domain=.true.)
    
        gdy(1:nlat) = tmp2(1,:)

        do ii=1,imt-1
           rdx(ii) = (gdx(ii+1) - gdx(ii))*deg2rad
        enddo

        rdx(0) = gdx(1)*deg2rad
        rdx(imt) = rdx(imt-1)*deg2rad
        rdx(imt+1) = rdx(imt)*deg2rad
    
        do ii=1,jmt-1
         rdy(ii) = (gdy(ii+1) - gdy(ii))*deg2rad
        enddo

        rdy(0) = gdy(1)*deg2rad
        rdy(jmt) = gdy(jmt-1)*deg2rad
        rdy(jmt+1) = gdy(jmt)*deg2rad

        call read_data(grid_file, 'rkmt', rkmt, no_domain=.true.)


        deallocate(tmp2)


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
                if (rkmt(ii,jj) .eq. 0.0 .and. rkmt(ii,jj+1) .eq. 0.0) &
                        rkmu(ii,jj) = 0.0
            enddo
        enddo
    
        do ii=1,imt-1
            do jj=1,jmt
                if (rkmt(ii,jj) .eq. 0.0 .and. rkmt(ii+1,jj) .eq. 0.0) &
                        rkmv(ii,jj) = 0.0
            enddo
        enddo
    
        do ii=1,imt-1
            do jj=1,jmt-1
                if (rkmu(ii,jj) .eq. 0.0 .and. rkmu(ii+1,jj) .eq. 0.0 .and. &
                    rkmv(ii,jj) .eq. 0.0 .and. rkmv(ii,jj+1) .eq. 0.0) &
                         rkmh(ii,jj) = 0.0
            enddo
        enddo

        do ii=1,imt
            do jj=1,jmt
                if (rkmt(ii,jj) .eq. 0.0) mask (ii,jj) = 1.0
                if (rkmt(ii,jj) .eq. 1.0) mask (ii,jj) = 0.0
            enddo
        enddo

        rrkmt(:,:) = rkmt(:,:)

        if (.not. field_exist(temp_clim_file, 'temp')) &
            call mpp_error(FATAL, 'field temp not found in '//trim(temp_clim_file))

        if (.not. field_exist(salt_clim_file, 'salt')) &
            call mpp_error(FATAL, 'field salt not found in '//trim(salt_clim_file))

        call field_size(temp_clim_file, 'temp', dimz)

        allocate ( u(imt,km,jmt,4), v(imt,km,jmt,4) ) 

        allocate ( t(imt,km,jmt,nn,4) ) 

        allocate ( h(imt,km,jmt,4) ) 
        
        allocate ( temp(imt,kmaxMYM,jmt,2), salt(imt,kmaxMYM,jmt,2) )

        allocate ( uvel(imt,kmaxMYM,jmt,2), vvel(imt,kmaxMYM,jmt,2) )

        id_restart = register_restart_field(restart_odtm, restart_file, 'u', u(:,:,:,1), u(:,:,:,2))
        id_restart = register_restart_field(restart_odtm, restart_file, 'v', v(:,:,:,1), v(:,:,:,2))
        id_restart = register_restart_field(restart_odtm, restart_file, 't', t(:,:,:,1,1), t(:,:,:,1,2), mandatory=.true.)
        id_restart = register_restart_field(restart_odtm, restart_file, 's', t(:,:,:,2,1), t(:,:,:,2,2), mandatory=.true.)
        id_restart = register_restart_field(restart_odtm, restart_file, 'h', h(:,:,:,1), h(:,:,:,2))
        id_restart = register_restart_field(restart_odtm, restart_file, 'temp', temp(:,:,:,1), temp(:,:,:,2), mandatory=.true.)
        id_restart = register_restart_field(restart_odtm, restart_file, 'salt', salt(:,:,:,1), salt(:,:,:,2), mandatory=.true.)
        id_restart = register_restart_field(restart_odtm, restart_file, 'uvel', uvel(:,:,:,1), uvel(:,:,:,2))
        id_restart = register_restart_field(restart_odtm, restart_file, 'vvel', vvel(:,:,:,1), vvel(:,:,:,2))

        we_upwel(:,:,:) = 0.0
        
        do kk=1,km
            h(:,kk,:,taum)=dz(kk)
            h(:,kk,:,taun)=h(:,kk,:,taum) 
            h(:,kk,:,taup)=h(:,kk,:,taum) 
        enddo

        we(:,:,:) = 0.0
        wd(:,:,:) = 0.0
        u(:,:,:,:) = 0.0
        v(:,:,:,:) = 0.0
        eta(:,:,:,:) = 0.0

        t(:,:,:,:,taum) = 10.0
        t(:,:,:,:,taun) = 10.0
       
 
        do nt = 1, 12
            call read_data(temp_clim_file, 'temp', temp_read(:,:,:,nt), timelevel=nt)
            call read_data(salt_clim_file, 'salt', salt_read(:,:,:,nt), timelevel=nt)
        enddo

        uvel(:,:,:,:) = 0.0
        vvel(:,:,:,:) = 0.0

        if (file_exist('INPUT/'//trim(restart_file))) then
            call restore_state(restart_odtm)
        else
            call mpp_error(NOTE,'Model starting from initial state')
            do i=1,imt
                do j=1,jmt
                    do k=1,201
                        tempin(k) = temp_read(i,k,j,1)  
                        saltin(k) = salt_read(i,k,j,1) 
                    enddo
                    do k=1,51
                        temp(i,k,j,1) = temp_read(i,k,j,1)
                        salt(i,k,j,1) = salt_read(i,k,j,1)
                        temp(i,k,j,2) = temp_read(i,k,j,1)
                        salt(i,k,j,2) = salt_read(i,k,j,1)
                    enddo

                    kmax = 201
                    do k=1,km-1   
                        call interp_extrap_initial (i,j,k,kmax,tempin, &
                            saltin,tempout,saltout)
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
        endif
  
    end subroutine init_grid

end program main
