
module size_mod
    
  !	grid specification for the model

  implicit none

  integer, parameter :: imt=180, jmt=120, nLayer=12, km=nLayer, lm=13
  integer, parameter :: NPP=225, nn=2, kmaxMYM=51

  double precision :: u(imt,km,jmt,4), v(imt,km,jmt,4)
  double precision :: t(imt,km,jmt,nn,4)
  double precision :: we(imt,km,jmt), wd(imt,km,jmt)
  double precision :: we_mld(imt,0:km,jmt)

  double precision :: we_prof(0:52), we_baro
  double precision :: uN(imt,km,2), vN(imt,km,2), hhN(imt,km,2), etaN(imt,km,2)
  double precision :: uS(imt,km,2),vS(imt,km,2),hhS(imt,km,2), etaS(imt,km,2)
  double precision :: uE(2,km,jmt), vE(2,km,jmt), hhE(2,km,jmt), etaE(2,km,jmt)
  double precision :: temp(imt,51,jmt,2), salt(imt,51,jmt,2)
  double precision :: uvel(imt,51,jmt,2), vvel(imt,51,jmt,2)

  real :: temp_read(imt,201,jmt,lm)
  real :: salt_read(imt,201,jmt,lm)
  real :: rK13(imt,kmaxMYM,jmt)
  real :: rK23(imt,kmaxMYM,jmt)
  real :: rK31(imt,kmaxMYM,jmt)
  real :: rK32(imt,kmaxMYM,jmt)
  real :: rK33(imt,kmaxMYM,jmt)
  real :: gdx(0:imt+1), gdy(0:jmt+1)
  real :: rdx(0:imt+1), rdy(0:jmt+1)
  real :: tracedubdxmld(51,4), tracedvbdymld(51,4)
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
  real :: h(imt,km,jmt,4), vmask(imt,jmt)
  real :: eta(imt,km,jmt,4)
  real :: pvort(imt,km,jmt)
  real :: pres_gradu(km), pres_gradv(km)
  real :: Asmag(imt,km,jmt)
  real :: DIFFuu(imt,km,jmt), DIFFvu(imt,km,jmt)
  real :: DIFFvv(imt,km,jmt), DIFFuv(imt,km,jmt)
  real :: age_psi(imt,km,jmt,nn)
  real :: age_tpsi(imt,km,jmt,nn)
  real :: age_t(imt,km,jmt,nn)
  real :: rEnergy(imt,km,jmt)
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
  real :: SMCoeff(imt,51,jmt), SHCoeff(imt,51,jmt)
  real :: diag_ext1(imt,51,jmt), diag_ext2(imt,51,jmt)
  real :: diag_ext3(imt,51,jmt), diag_ext4(imt,51,jmt)
  real :: diag_ext5(imt,51,jmt), diag_ext6(imt,51,jmt)
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
  real :: Asmag_back
  integer :: itimer2, itimer3, itimerrate, itimermax
  real :: time_switch
  real :: dz(km), dz_max(km), dz_min(km), zdz(km)
  real :: dzmym(kmaxMYM), zdepth(kmaxMYM)
  real :: tr01_max(km), tr02_max(km), tr01_min(km), tr02_min(km)
  real :: fcor(imt,jmt), fcoruh, fcorvh
  real :: rho(km), rho_min(km), rho_max(km)

#ifdef output_average
	real :: SHCoeff_output(imt,51,jmt), SMCoeff_output(imt,51,jmt)
#ifdef snap_mld_extended
#endif
#endif
end module size_mod
