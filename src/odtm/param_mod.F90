module param_mod
    !   parameter initialization for model
    use size_mod, only : lm

    implicit none

    real :: diffuse_uu, diffuse_vv, diffuse_h

    real :: advec_uu, advec_vu, advec_uv, advec_vv
    integer :: loop_day, loop_total
    real :: sum_adv
    real :: tropdubdx, tropdvbdy
    integer :: nmid
    real :: rnmid, rvalid_grid_background, rvalid_grid
    real :: rvalid_grid_mass, rvalid_grid_background_mass
    real :: fracx, fracy

    integer :: loop_ind

    !c  time stamps for netcdf snapshots.
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real :: dpm(lm) = (/ 31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, &
                         31.0, 30.0, 31.0, 30.0, 31.0, 31.0/)


    real, parameter :: g= 9.8 !0.0098 !0.0196 !0.090 !0.03
    real, parameter :: pi =22.0/7.0
    real, parameter :: deg2rad = pi/180.0
    real, parameter :: rad2deg = 1./deg2rad
    real, parameter :: omega = 0.7272205e-4
    real, parameter :: Re = 6368.0e3
    real, parameter :: deg2meter = deg2rad*Re
    real, parameter :: dxd= 1.0/2.0
    real, parameter :: dyd= 1.0/2.0
    real, parameter :: reflat_start = -15.0

    real, parameter :: dx=dxd*deg2meter
    real, parameter :: dy=dyd*deg2meter
    real :: dx2=dx*dx
    real, parameter :: dy2=dy*dy
    real :: dphi=dxd*deg2rad
    real :: dthe=dyd*deg2rad

    real, parameter :: dt= 1800 
    real, parameter :: dtts=2.0*dt
    real, parameter :: day2sec=86400.0
    real, parameter :: sec2day = 1/day2sec
    real, parameter :: alpha = 0.1
    real :: diffuse_MY = 1.0e3
    
    real, parameter :: diffuse = 1000
    real, parameter :: diffuse_th = 0.0e3
    real :: diffuse_tr = 1.0e3

    real, parameter :: te = day2sec*1
    real, parameter :: td = day2sec*30

	real :: polar_u, polar_v

end module param_mod




