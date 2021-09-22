
module size_mod
    
  !	grid specification for the model

  implicit none

  integer, parameter :: nLayer=6, km=nLayer, lm=13
  integer, parameter :: NPP=225, nn=2, kmaxMYM=51, kclim=201
  integer, parameter :: taum = 1, taun = 2, taup = 3, taus = 4
  
  integer :: imt, jmt, halo=1
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: i,j,k
  integer :: iday_wind, iday_start, iday_start_snap
  integer :: itimer2, itimer3, itimerrate, itimermax
  integer :: loop, month, month_start, loop_start, lpm, lpd
  integer :: month_wind, month_start_snap
  integer :: taum1,taun1,taup1
  real :: Q_atmos
  real :: fcoruh, fcorvh
  real :: time_switch
  real :: tracer_switch
  real :: we_baro
  real :: wind_switch

  real, pointer, dimension(:,:,:,:) :: u => null(), v => null()

  real, pointer, dimension(:,:,:,:,:) :: t => null()

  real, pointer, dimension(:,:,:,:) :: temp => null(), salt => null()
  real, pointer, dimension(:,:,:,:) :: uvel => null(), vvel => null()
  real, pointer, dimension(:,:,:,:) :: h => null()

  real, pointer, dimension(:,:,:) :: diag_ext1 => null(), diag_ext2 => null()
  real, pointer, dimension(:,:,:) :: diag_ext3 => null(), diag_ext4 => null()
  real, pointer, dimension(:,:,:) :: diag_ext5 => null(), diag_ext6 => null()

  real, pointer, dimension(:,:,:,:) :: temp_read => null(), salt_read => null()

  real, pointer, dimension(:,:,:) :: we => null(), wd => null(), we_mld => null()

  real, pointer, dimension(:,:,:,:) :: eta => null()

  real, pointer, dimension(:,:,:) :: pvort => null()

  real, pointer, dimension(:,:,:) :: rEnergy => null()

  real, pointer, dimension(:,:,:) :: SMCoeff => null(), SHCoeff => null()

  real, pointer, dimension(:,:,:) :: rmld_misc => null(), denss => null()

  real, pointer, dimension(:,:) :: mld_mld => null()
  integer, pointer, dimension(:,:) :: mask => null()
  real, pointer, dimension(:) :: we_prof => null()

  real, pointer, dimension(:) :: gdx => null(), gdy => null()
  real, pointer, dimension(:) :: gdxb => null(), gdyb => null()
  real, pointer, dimension(:) :: rdx => null(), rdy => null()

  real, pointer, dimension(:,:) :: tracedubdxmld => null(), tracedvbdymld => null()
  real, pointer, dimension(:,:) :: dxu => null(), dyu => null()
  real, pointer, dimension(:,:) :: dxv => null(), dyv => null()
  real, pointer, dimension(:,:) :: dxh => null(), dyh => null()
  real, pointer, dimension(:,:) :: dah => null()
  real, pointer, dimension(:,:) :: dau => null()
  real, pointer, dimension(:,:) :: dav => null()
  real, pointer, dimension(:,:) :: rdxu => null(), rdyu => null()
  real, pointer, dimension(:,:) :: rdxv => null(), rdyv => null()
  real, pointer, dimension(:,:) :: rdxh => null(), rdyh => null()

  real, pointer, dimension(:,:) :: rkmt => null(), rkmh => null()
  real, pointer, dimension(:,:) :: rkmu => null(), rkmv => null()
  real, pointer, dimension(:,:) :: rrkmt => null()
  logical, pointer, dimension(:,:) :: omask => null()

  real, pointer, dimension(:) :: pres_gradu => null(), pres_gradv => null()

  real, pointer, dimension(:,:) :: wd_mask => null()
  real, pointer, dimension(:,:,:) :: we_upwel => null()

  real, pointer, dimension(:,:,:) :: taux => null(), tauy => null()
  real, pointer, dimension(:,:) :: taux_force => null(), tauy_force => null()
  real, pointer, dimension(:,:) :: ssw => null()
  real, pointer, dimension(:,:) :: cld => null()
  real, pointer, dimension(:,:) :: pme => null(), pme_corr => null()
  real, pointer, dimension(:,:) :: chl => null()
  real, pointer, dimension(:,:) :: rvr => null()

  real, pointer, dimension(:,:) :: sphm => null(), airt => null()
  real, pointer, dimension(:,:) :: uwnd => null(), vwnd => null(), fcor => null()
  real, pointer, dimension(:,:) :: fcoru => null(), fcorv => null()

  real, pointer, dimension(:) :: he => null(), hd => null()
  real, pointer, dimension(:) :: dz => null(), dz_max => null(), dz_min => null(), zdz => null()
  real, pointer, dimension(:) :: dzmym => null(), zdepth => null()
  real, pointer, dimension(:) :: tr01_max => null(), tr02_max => null(), tr01_min => null(), tr02_min => null()
  real, pointer, dimension(:) :: rho => null(), rho_min => null(), rho_max => null()
  real, pointer, dimension(:) :: phi => null(), theta => null(), tan_theta => null()

    contains

    
    subroutine init_size()

        real :: rsum
        integer :: isg, ieg, jsg, jeg

        isg = 1 - halo
        ieg = imt + halo
        jsg = 1 - halo
        jeg = jmt + halo
 
        ! on data domain 
        allocate ( u(isd:ied,jsd:jed,km,4), v(isd:ied,jsd:jed,km,4) )
        allocate ( t(isd:ied,jsd:jed,km,nn,4) )
        allocate ( h(isd:ied,jsd:jed,km,4) )
        allocate ( temp(isd:ied,jsd:jed,kmaxMYM,2), salt(isd:ied,jsd:jed,kmaxMYM,2) )
        allocate ( uvel(isd:ied,jsd:jed,kmaxMYM,2), vvel(isd:ied,jsd:jed,kmaxMYM,2) ) 
        allocate ( rkmh(isd:ied,jsd:jed))
        allocate ( rkmu(isd:ied,jsd:jed), rkmv(isd:ied,jsd:jed) )


        ! on compute domain
        allocate ( fcor(isc:iec,jsc:jec) )
        allocate ( fcoru(isc:iec,jsc:jec) )
        allocate ( fcorv(isc:iec,jsc:jec) )
        allocate ( diag_ext1(isc:iec,jsc:jec,kmaxMYM), diag_ext2(isc:iec,jsc:jec,kmaxMYM) )
        allocate ( diag_ext3(isc:iec,jsc:jec,kmaxMYM), diag_ext4(isc:iec,jsc:jec,kmaxMYM) )
        allocate ( diag_ext5(isc:iec,jsc:jec,kmaxMYM), diag_ext6(isc:iec,jsc:jec,kmaxMYM) )

        allocate ( temp_read(isc:iec,jsc:jec,kclim,lm) )
        allocate ( salt_read(isc:iec,jsc:jec,kclim,lm) )

        allocate ( we(isd:ied,jsd:jed,km), wd(isd:ied,jsd:jed,km), we_mld(isc:iec,jsc:jec,0:km) )
        allocate ( eta(isc:iec,jsc:jec,km,4), pvort(isc:iec,jsc:jec,km), rEnergy(isc:iec,jsc:jec,km) )
        allocate ( SMCoeff(isc:iec,jsc:jec,kmaxMYM), SHCoeff(isc:iec,jsc:jec,kmaxMYM) )
        allocate ( rmld_misc(isc:iec,jsc:jec,kmaxMYM), denss(isc:iec,jsc:jec,km) )

        allocate ( mld_mld(isc:iec,jsc:jec) )
        allocate ( wd_mask(isc:iec, jsc:jec) )
        allocate ( we_upwel(isc:iec, jsc:jec, 2) )
        allocate ( taux(isc:iec,jsc:jec,2), tauy(isc:iec,jsc:jec,2) )
        allocate ( taux_force(isc:iec,jsc:jec), tauy_force(isc:iec,jsc:jec) )
        allocate ( ssw(isc:iec,jsc:jec) )
        allocate ( cld(isc:iec,jsc:jec) )
        allocate ( pme(isc:iec,jsc:jec), pme_corr(isc:iec,jsc:jec) )
        allocate ( chl(isc:iec,jsc:jec) )
        allocate ( rvr(isc:iec,jsc:jec) )
        allocate ( sphm(isc:iec,jsc:jec), airt(isc:iec,jsc:jec) )
        allocate ( uwnd(isc:iec,jsc:jec), vwnd(isc:iec,jsc:jec) )


        ! only on vertical
        allocate ( we_prof(0:kmaxMYM+1) )
        allocate ( tracedubdxmld(kmaxMYM,4), tracedvbdymld(kmaxMYM,4) )
        allocate ( pres_gradu(km), pres_gradv(km) )
        allocate ( he(km), hd(km) )
        allocate ( dz(km), dz_max(km), dz_min(km), zdz(km) )
        allocate ( dzmym(kmaxMYM), zdepth(kmaxMYM) )
        allocate ( tr01_max(km), tr02_max(km), tr01_min(km), tr02_min(km) )
        allocate ( rho(km), rho_min(km), rho_max(km) )

        ! On global Domain
        allocate ( phi(imt), theta(jmt) )
        allocate ( tan_theta(jmt) )
        allocate ( gdx(imt), gdy(jmt) )
        allocate ( gdxb(imt+1), gdyb(jmt+1) )
        allocate ( rdx(imt), rdy(jmt) )

!        allocate ( dxu(imt,jmt), dyu(imt,jmt) )
!        allocate ( dxv(imt,jmt), dyv(imt,jmt) )
!        allocate ( dxh(imt,jmt), dyh(imt,jmt) )
!
!        allocate ( dah(imt,jmt) )
!        allocate ( dau(imt,jmt) )
!        allocate ( dav(imt,jmt) )

        allocate ( dxu(isd:ied,jsd:jed), dyu(isd:ied,jsd:jed) )
        allocate ( dxv(isd:ied,jsd:jed), dyv(isd:ied,jsd:jed) )
        allocate ( dxh(isd:ied,jsd:jed), dyh(isd:ied,jsd:jed) )

        allocate ( dah(isd:ied,jsd:jed) )
        allocate ( dau(isd:ied,jsd:jed) )
        allocate ( dav(isd:ied,jsd:jed) )

        allocate ( rdxu(isg:ieg,jsg:jeg), rdyu(isg:ieg,jsg:jeg) )
        allocate ( rdxv(isg:ieg,jsg:jeg), rdyv(isg:ieg,jsg:jeg) )
        allocate ( rdxh(isg:ieg,jsg:jeg), rdyh(isg:ieg,jsg:jeg) )
        allocate ( omask(isg:ieg,jsg:jeg) )
        allocate ( mask(imt,jmt) )
        allocate ( rkmt(imt,jmt) ) 

        dz(1) = 70.0     ! layer-1
        dz(2) = 70.0     ! layer-2
        dz(3) = 70.0     ! layer-3
        dz(4) = 190.0
        dz(5) = 400.0
        dz(6) = 3200.0


        zdz(1) = dz(1) - dz(1)/2.0
        rsum = dz(1)

        do k=2,km
            rsum = rsum + dz(k)
            zdz(k) = rsum - dz(k)/2.0
        enddo


        rho(1) = 1023.66
        rho(2) = 1025.11
        rho(3) = 1026.40
        rho(4) = 1028.00
        rho(5) = 1030.00
        rho(6) = 1032.00

        rho_max(1) = 1028.31
        rho_max(2) = 1028.69
        rho_max(3) = 1029.19
        rho_max(4) = 1030.00
        rho_max(5) = 1033.00
        rho_max(6) = 1035.00

        rho_min(1) = 1020.02
        rho_min(2) = 1022.53
        rho_min(3) = 1024.43
        rho_min(4) = 1026.00
        rho_min(5) = 1028.00
        rho_min(6) = 1030.00

        tr01_max(1) = 31.5
        tr01_max(2) = 29.0
        tr01_max(3) = 25.5
        tr01_max(4) = 22.5
        tr01_max(5) = 20.5

        tr01_min(1) = 10.5
        tr01_min(2) = 09.8
        tr01_min(3) = 08.4
        tr01_min(4) = 06.4
        tr01_min(5) = 05.0

        tr02_max(1) = 40.2
        tr02_max(2) = 40.3
        tr02_max(3) = 40.5
        tr02_max(4) = 40.6
        tr02_max(5) = 40.7

        tr02_min(1) = 31.7
        tr02_min(2) = 33.5
        tr02_min(3) = 33.6
        tr02_min(4) = 33.8
        tr02_min(5) = 34.0

! he is the threshold height of layer at and below which entrainment starts
        he(1) = 30
        he(2) = 30
        he(3) = 30
        he(4) = 50
        he(5) = 80
        he(6) = 1000

! hd is the ht of layer, h>=hd- detrainment present
        hd(1) = 300.0 ! what is hd?
        hd(2) = 300.0
        hd(3) = 300.0
        hd(4) = 800.0
        hd(5) =1000.0


        dz_max(1) = 3*dz(1) ! what is dz_max?
        dz_max(2) = 3*dz(2)
        dz_max(3) = 3*dz(3)
        dz_max(4) = 3*dz(4)
        dz_max(5) = 3*dz(5)
        dz_max(6) = 3*dz(6)

        dz_min(1) = 10
        dz_min(2) = 10
        dz_min(3) = 10
        dz_min(4) = 30
        dz_min(5) = 30
        dz_min(6) = 1000

    end subroutine init_size

end module size_mod
