module param_mod
    !   parameter initialization for model
    use size_mod, only : imt, jmt, km, NPP, nn, lm

    implicit none

    double precision :: diffuse_uu, diffuse_vv, diffuse_h, diffuse_eta = 2.5e3
    double precision :: diffuse_uv, diffuse_vu, diffuse_u, diffuse_v

    real :: adv_u, adv_v
    real :: adv_t, adv_s
    real :: visc_hv, visc_hu
    real :: gdash
    real :: vert_mixu, vert_mixv
    real :: D, d2z, dz2, curl
    real :: advec_u, advec_v
    real :: advec_uu, advec_vu, advec_uv, advec_vv
    real :: gdash_advec_x, gdash_advec_y
    real :: dens_advec_x, dens_advec_y
    real :: advec_t, diffuse_t
    real :: tr_source_term
    real :: Eintrmdf, Eintrmdb
    real :: biharmon_u, biharmon_v, biharmon_h
    integer :: loop_day, loop_month, loop_total
    real :: sum_adv
    real :: ht
    real :: tropdubdx, tropdvbdy
    integer :: nmid
    real :: rvisc_4u, rvisc_4v
    real :: bfr
    real :: rnmid, rvalid_grid_background, rvalid_grid
    real :: rvalid_grid_mass, rvalid_grid_background_mass
    real :: rvolume, rvolume_background
    real :: fracx, fracy
    real :: xold, yold

    integer :: xposition, yposition
    integer :: xmone, xpone, ymone, ypone
    real :: sum_accr, sum_accl, sum_accu, sum_accd, dens_ave
    real :: pot_temp, pot_den
    real :: sum_acc, dens350P
    real :: ave
    logical :: NaN_check

!c===================== snaps.h ===========================================

    integer, parameter :: itime=1
    double precision :: taux_1, tauy_1
    integer :: ncid,vlonid,vlatid,vdepthid,vtimeid,vtimeid1,ncid1
    integer :: vlonid1, vlatid1, vdepthidmld, ncid2
    integer :: start4(4)
    integer :: startp(2)
    integer :: startinv(2)
    integer :: loop_ind
    integer :: vardim(4), vtid, vsid, vuid, vvid, vardim3(3)
    integer :: vetaid
    integer :: vardim33(3), vardim2(2), vardimmld(4), vardimmld3(3)
    integer :: varpdim(2), vardimmldext(4)
    integer :: varinvdim(2)
    integer :: vqnetid, vqsaltid, vtauxid, vtauyid
    integer :: vatmos_h_id, vatmos_u_id, vatmos_v_id
    integer :: vsaltid, vdensid, vvortid
    integer :: vtempid(nn)
    integer :: vmldtempid, vmldsaltid, vmlduvelid, vmldvvelid, vmldmldid
    integer :: vmldshid, vmldsmid, vmldqnetid
    integer :: vpxid, vpyid
    integer :: vdxuid, vdyvid
    integer :: vthetaid
    integer :: vmaskid
    integer :: vmld4did1, vmld4did2, vmld4did3
    integer :: vmld4did4, vmld4did5, vmld4did6
    integer :: vtimeid2, lonid2, latid2, idepthidmld2
    integer :: itimeid2, vlonid2, vlatid2, vdepthidmld2
    real :: rlon(imt),rlat(jmt),rdepth(km),rtime(itime)
    real :: uu(imt,km,jmt), vv(imt,km,jmt), rdepth_mld(51)
    real :: hh(imt,km,jmt) 
    integer :: mask(imt,jmt)
    real :: tauxx(imt,jmt), tauyy(imt,jmt)
    real :: u_atmos_snap (imt,jmt)
    real :: v_atmos_snap (imt,jmt)
    real :: h_atmos_snap (imt,jmt)

    real :: temp_snap(imt,km,jmt), salt_snap(imt,km,jmt)
    real :: temp_mld(imt,51,jmt), salt_mld(imt,51,jmt)
    real :: uvel_mld(imt,51,jmt), vvel_mld(imt,51,jmt)
    double precision :: SMCoeff_mld(imt,51,jmt), SHCoeff_mld(imt,51,jmt)
    real :: rmld_misc(imt,51,jmt), denss(imt,km,jmt)
    real :: diag_ext1_mld(imt, 51, jmt), diag_ext2_mld(imt, 51, jmt)
    real :: diag_ext3_mld(imt, 51, jmt), diag_ext4_mld(imt, 51, jmt)
    real :: diag_ext5_mld(imt, 51, jmt), diag_ext6_mld(imt, 51, jmt)
    real :: mld_mld(imt,jmt)

    real :: pvorticity(imt,km,jmt)
    real :: undef

    character(len=21) :: stamp_day(12) = (/ &
                    'days since 1995-01-01', &
                    'days since 1995-02-01', &
                    'days since 1995-03-01', &
                    'days since 1995-04-01', &
                    'days since 1995-05-01', &
                    'days since 1995-06-01', &
                    'days since 1995-07-01', &
                    'days since 1995-08-01', &
                    'days since 1995-09-01', &
                    'days since 1995-10-01', &
                    'days since 1995-11-01', &
                    'days since 1995-12-01'/)
         
    character(len=20) :: stamp_time(12) = (/ &
                    '01-JAN-1995 00:00:00', &
                    '01-FEB-1995 00:00:00', &
                    '01-MAR-1995 00:00:00', &
                    '01-APR-1995 00:00:00', &
                    '01-MAY-1995 00:00:00', &
                    '01-JUN-1995 00:00:00', &
                    '01-JUL-1995 00:00:00', &
                    '01-AUG-1995 00:00:00', &
                    '01-SEP-1995 00:00:00', &
                    '01-OCT-1995 00:00:00', &
                    '01-NOV-1995 00:00:00', &
                    '01-DEC-1995 00:00:00'/)

    character(len=4) :: Tr(nn)

    character(len=21) :: stamp_day_daily(365)
    character(len=21) :: stamp_time_daily(365)

    integer :: number_of_snap

    integer :: start3(3)=(/1,1,1/)
    integer :: start2(2)=(/1,1/)
    integer :: count2(2)=(/imt,jmt/)
    integer :: count4(4)=(/imt,km,jmt,1/)
    integer :: count3(3)=(/imt,jmt,itime/)
    integer :: count33(3)=(/imt,km,jmt/)
    integer :: count4mld(4)=(/imt,51,jmt,1/)
    integer :: countp(2)=(/NPP, 1/)
    integer :: countinv(2)=(/nn, 1/)

    !c  time stamps for netcdf snapshots.
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real :: dpm(lm) = (/ 31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, &
                         31.0, 30.0, 31.0, 30.0, 31.0, 31.0/)


    real, parameter :: g= 9.8 !0.0098 !0.0196 !0.090 !0.03
    real, parameter :: dens_ref = 1025.0
    real, parameter :: pi =22.0/7.0
    real, parameter :: pi2 = pi*pi
    real, parameter :: deg2rad = pi/180.0
    real, parameter :: omega = 0.7272205e-4
    real, parameter :: radearth = 6400.0e3
    real, parameter :: Re = 6368.0e3
    real, parameter :: deg2meter = deg2rad*Re
    real, parameter :: reflon = 30.0 !9.50 ! 29.50  !for indiand ocen
    real, parameter :: reflat_start = -15.0 !-40.0 !-35.25 !-40.25

    real, parameter :: dxd= 1.0/2.0
    real, parameter :: dyd= 1.0/2.0
    real, parameter :: reflat = -0.25 - (-reflat_start - (jmt*dyd + reflat_start))/2.0 

    real, parameter :: dx=dxd*deg2meter
    real, parameter :: dy=dyd*deg2meter
    real, parameter :: d2x=0.7*dx
    real, parameter :: d2y=0.7*dy
    real :: dx2=dx*dx
    real, parameter :: dy2=dy*dy
    real, parameter :: dx4=dx*dx*dx*dx
    real, parameter :: dy4=dy2*dy2
    !real, parameter :: dphi=dxd*deg2rad
    real :: dphi=dxd*deg2rad
    !real, parameter :: dthe=dyd*deg2rad
    real :: dthe=dyd*deg2rad
    real, parameter :: d2phi=2.0*dxd*deg2rad
    real, parameter :: d2the=2.0*dyd*deg2rad

    real, parameter :: dt= 1800 !600
    real :: dtts=2.0*dt
    real, parameter :: day2sec=86400.0
    real, parameter :: sec2day = 1/day2sec
    real, parameter :: alpha = 0.1
    
    real, parameter :: eddy_visc=4.0e5
    real, parameter :: eddy_mix=1.0e-1
    real, parameter :: beta = 2.2891e-11  ! beta at equator = 2.2891e-11, at 30= 2.0e-11
    real, parameter :: f0 = 0.0
    real, parameter :: dissip = 2.0e-7
    real, parameter :: diffuse = 1000
    real, parameter :: diffuse_th = 0.0e3
    real :: diffuse_tr = 1.0e3
    real :: diffuse_MY = 1.0e3
    real, parameter :: diffuse_ah = 2.5e3
    real, parameter :: Asmag_back = diffuse
    real, parameter :: smag_coeff = 0.3
    real, parameter :: biharmonic = 1.0e13
    real, parameter :: dens_relax = 1/(2.0*dt*4)
    real, parameter :: trelax = 1.0/(6.0*2.0*dt)

    real, parameter :: te = day2sec*1
    real, parameter :: td = day2sec*30

	real :: phi(imt), theta(jmt)
	real :: dphi2, dthe2
	real :: polar_u, polar_v
	real :: polar_diff_u(imt,km,jmt) 
	real :: polar_diff_v(imt,km,jmt)

end module param_mod




