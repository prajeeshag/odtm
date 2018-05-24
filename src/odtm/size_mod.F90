
module size_mod
    
  !	grid specification for the model

  implicit none

  integer, parameter :: imt=180, jmt=120, nLayer=12, km=nLayer, lm=13
  integer, parameter :: NPP=225, nn=2, kmaxMYM=51, kclim=201

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

  real :: mld_mld(imt,jmt)
  integer :: mask(imt,jmt)
  real :: we_prof(0:kmaxMYM+1), we_baro

 ! real :: uN(imt,km,2), vN(imt,km,2), hhN(imt,km,2), etaN(imt,km,2)
 ! real :: uS(imt,km,2), vS(imt,km,2), hhS(imt,km,2), etaS(imt,km,2)
 ! real :: uE(2,km,jmt), vE(2,km,jmt), hhE(2,km,jmt), etaE(2,km,jmt)

  real :: gdx(0:imt+1), gdy(0:jmt+1)
  real :: rdx(0:imt+1), rdy(0:jmt+1)

  real :: tracedubdxmld(kmaxMYM,4), tracedvbdymld(kmaxMYM,4)
  real :: dxr(imt,jmt), dyr(imt,jmt)
  real :: dxu(imt,jmt), dyu(imt,jmt)
  real :: dxv(imt,jmt), dyv(imt,jmt)
  real :: dxh(imt,jmt), dyh(imt,jmt)
  real :: dah(imt,jmt)

  integer :: i,j,k

  real :: rkmt(imt,jmt), rkmh(imt,jmt), rkmhH(imt,jmt)
  real :: rkmu(imt,jmt), rkmv(imt,jmt)
  real :: rrkmt(imt,jmt), tmask(imt,jmt)

  real :: smask(imt,jmt), umask(imt,jmt)
  real :: vmask(imt,jmt)

  real :: pres_gradu(km), pres_gradv(km)

  real :: wd_mask(imt, jmt)
  real :: we_upwel(imt, jmt, 2)
  real :: tr_source_term_2D(imt, jmt)
  real :: taux(imt,jmt,2), tauy(imt,jmt,2)
  real :: taux_force(imt,jmt), tauy_force(imt,jmt)
  real :: taux_snap(imt,jmt), tauy_snap(imt,jmt)
  real :: sphm_read(imt,jmt,2), airt_read(imt,jmt,2)
  real :: uwnd_read(imt,jmt,2), vwnd_read(imt,jmt,2)
  real :: ssw_read(imt,jmt,2), ssw(imt,jmt)
  real :: cld_read(imt,jmt,2), cld(imt,jmt)
  real :: pme_read(imt,jmt,2), pme(imt,jmt), pme_corr(imt,jmt)
  real :: chl_read(imt,jmt,2), chl(imt,jmt)
  real :: rvr_read(imt,jmt,2), rvr(imt,jmt)

  real :: sphm(imt,jmt), airt(imt,jmt)
  real :: uwnd(imt,jmt), vwnd(imt,jmt)
  real :: Thetainv(nn)
  integer :: taum,taun,taup,taus
  integer :: taum1,taun1,taup1
  real :: Q_atmos
  real :: he(km), hd(km)
  real :: days
  integer :: loop, month, month_start, loop_start, lpm, lpd
  integer :: month_wind, month_start_snap
  real :: wind_switch
  integer :: iday_wind, iday_start, iday_start_snap
  real :: tracer_switch
  integer :: itimer2, itimer3, itimerrate, itimermax
  real :: time_switch

  real :: dz(km), dz_max(km), dz_min(km), zdz(km)
  real :: dzmym(kmaxMYM), zdepth(kmaxMYM)
  real :: tr01_max(km), tr02_max(km), tr01_min(km), tr02_min(km)
  real :: fcor(imt,jmt), fcoruh, fcorvh
  real :: rho(km), rho_min(km), rho_max(km)


    contains

    
    subroutine init_size()
        
        allocate ( u(imt,jmt,km,4), v(imt,jmt,km,4) )
        allocate ( t(imt,jmt,km,nn,4) )
        allocate ( h(imt,jmt,km,4) )
        allocate ( temp(imt,jmt,kmaxMYM,2), salt(imt,jmt,kmaxMYM,2) )
        allocate ( uvel(imt,jmt,kmaxMYM,2), vvel(imt,jmt,kmaxMYM,2) ) 

        allocate ( diag_ext1(imt,jmt,kmaxMYM), diag_ext2(imt,jmt,kmaxMYM) )
        allocate ( diag_ext3(imt,jmt,kmaxMYM), diag_ext4(imt,jmt,kmaxMYM) )
        allocate ( diag_ext5(imt,jmt,kmaxMYM), diag_ext6(imt,jmt,kmaxMYM) )

        allocate ( temp_read(imt,jmt,kclim,lm) )
        allocate ( salt_read(imt,jmt,kclim,lm) )

        allocate ( we(imt,jmt,km), wd(imt,jmt,km), we_mld(imt,jmt,0:km) )
        allocate ( eta(imt,jmt,km,4), pvort(imt,jmt,km), rEnergy(imt,jmt,km) )
        allocate ( SMCoeff(imt,jmt,kmaxMYM), SHCoeff(imt,jmt,kmaxMYM) )

        allocate ( rmld_misc(imt,jmt,kmaxMYM), denss(imt,jmt,km) )

    end subroutine init_size

end module size_mod
