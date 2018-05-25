         subroutine mixed_layer_physics
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       subroutine to parameterise mixing
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : airt, chl, cld, dzmym, gdy, i, j, pme, pme_corr
            use size_mod, only : rkmh, rvr, salt, sphm, ssw, temp, uvel, uwnd
            use size_mod, only : vvel, vwnd, we_baro, we_prof, zdepth, loop
            use size_mod, only : salt_read, SHCoeff, SMCoeff, temp_read, month
            use size_mod, only : diag_ext1, diag_ext2, diag_ext3
            use size_mod, only : diag_ext4, diag_ext5, diag_ext6
            use size_mod, only : rmld_misc, mld_mld, imt, jmt
            use param_mod, only : day2sec, diffuse_my, diffuse_tr, dt, dtts

            use interp_extrap_mod, only : interp_extrap_we_profile
            use heat_flux_mod, only : heat_flux
            use moto_mod, only : MT1

        implicit none
        real tempin(51), saltin(51)
        real rmldout
        real uvelin(51), vvelin(51), we_profin(0:52)
        real sswin(51)
        real SM(51), SH(51), YQ(51)
        real tauxin, tauyin, entr_flux(51)
        real diag1(51), diag2(51), diag3(51)
        real diag4(51), diag5(51), diag6(51)
        real tx, ty, Lv, pme_eps
        integer tau_index
      
        integer :: kmax, kkk, ii, jj, kk, icval
        real :: albedo, sfx, rAlpha2, rL1, rchla, rPAR, rinfraRed, rvisible
        real :: rvisible_extinct, cldin, uwndin, vwndin, speed, airt_in, sphm_in
        real :: sst_in, qs, ql, qb, rAlpha1, rL2, qnet_in, pme_in, rLATin
        real :: rLATMax, rvr_runoff, relax_South, relax_East, relax_on_off
        real :: alpha_rgm, alpha_rgt, beta_mldm, beta_mldt, gama, zero_out_tradv
        real :: rdrag, widS, widE
        !! Here is the start of Mixed Layer Model


         do i=1,imt
         do j=1,jmt

         kmax = 51
         dzmym(1) = 5.0
         zdepth(1) = 5.0
         do kkk=2,kmax
         dzmym(kkk) = 5.0
         zdepth(kkk) = zdepth(kkk-1) + 5.0
         enddo

         if (rkmh(i,j) .eq. 1.0 ) then

        ii = i
        jj = j
        call interp_extrap_we_profile (ii, jj, kmax)
c        pvort(i,km,j) = we_baro
         do kk = 1,kmax
         ! print *, 'test=', temp(i,kk,j,1), rkmh(i,j)
         tempin(kk) = temp(i,j,kk,1) *rkmh(i,j)
         saltin(kk) = salt(i,j,kk,1) *rkmh(i,j)
         uvelin(kk) = uvel(i,j,kk,1) *rkmh(i,j)
         vvelin(kk) = vvel(i,j,kk,1) *rkmh(i,j)
         we_profin(kk) = we_prof(kk)
         enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! treatement for surface light
        sswin(1)  = ssw(i,j)
        albedo     = 0.3
        sfx =  sswin(1)  *(1.0-albedo)     !rsum/ncount  !* (1 -mod(loop,2))
        do kk = 2,kmax
        rAlpha1    =       0.805
        rAlpha2    =       0.195
        rL1        = 0.043
        rL2        = 0.42
        rchla      = max(0.02,chl(i,j)) !0.02  ! mg/m3
        rPAR = 0.45
!        rsum = 0.0
!        ncount = 0
!         do kkk=2,kk
c         rsum = rsum + sswin(kkk-1)
c         ncount = ncount + 1
c         enddo
        rinfraRed =
     & sfx * rPAR * ( rAlpha1 * exp(-1.0d0 * (kk-1) * 5.0 * rL1)+
     &                 rAlpha2 * exp(-1.0d0 * (kk-1) * 5.0 * rL2) )
        rvisible  = sfx * (1-rPAR) * rvisible_extinct (rchla, kk-1)
        sswin(kk) = rinfraRed * 0.0 + rvisible
        enddo


! treatement for surface light
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         cldin= cld(i,j) * 0.0
         uwndin = uwnd(i,j)
         vwndin = vwnd(i,j)
         speed = (uwndin**2 + vwndin**2)**0.5
         airt_in = airt(i,j)
         sphm_in = sphm(i,j)
         sst_in = tempin(1) ! temperature of first layer
         call heat_flux(sphm_in,airt_in,sst_in,uwndin,vwndin,cldin,qnet_in,
     &   qs, ql, qb, tx, ty, Lv)
         pme_in  = (pme(i,j) + rvr(i,j)*0.0 + (ql/(Lv*1024.5))) * (-1.0)
         pme_eps = pme_corr(i,j)
         rLATin = gdy(j)
         rLATMax = 0.75
         if (abs(rLATin) .eq. rLATMax)
     &   rLATin = rLATMax*sign(rLATin , rLATMax)

         rvr_runoff = rvr(i,j)
         tauxin = tx
         tauyin = ty
         qnet_in = qnet_in !* (1 - mod(loop,2))

         rmld_misc(i,j,1)  = qnet_in
         rmld_misc(i,j,2)  = ssw(i,j)
         rmld_misc(i,j,3)  = pme_in
         rmld_misc(i,j,4)  = uwndin
         rmld_misc(i,j,5)  = vwndin
         rmld_misc(i,j,6)  = tauxin
         rmld_misc(i,j,7)  = tauyin
         rmld_misc(i,j,8)  = sphm_in
         rmld_misc(i,j,9)  = airt_in
         rmld_misc(i,j,10)  = cldin
         rmld_misc(i,j,11)  = chl(i,j)
         rmld_misc(i,j,12)  = ql
         rmld_misc(i,j,13)  = qs
         rmld_misc(i,j,14)  = qb
         rmld_misc(i,j,15)  = rvr(i,j)


         call MT1 (tempin, saltin, uvelin, vvelin, tauxin,
     &    tauyin, qnet_in, sswin, pme_in,
     &    rLATin, dt, SM, SH, we_profin, we_baro, icval, YQ, entr_flux,
     &    speed, pme_eps, rvr_runoff,
     &    diag1, diag2, diag3, diag4, diag5, diag6,ii,jj)
          
       
!        do kk=1,kmax
!        if (i .eq. 44 .and. j .eq. 31) write(*,*) tempin(1), saltin(1), 
!     &    uvelin(1), vvelin(1), tauxin, tauyin, qnet_in, sswin(1), pme_in,
!     &    rLATin, dt, SM(1), SH(1), we_profin(1), we_baro, icval, YQ(1), 
!     &    entr_flux(1),
!     &    speed, pme_eps, rvr_runoff

!        enddo

         do kk = 1,kmax
!         rmld_misc(i,kk,j) = entr_flux(kk)
         temp(i,j,kk,2) = tempin(kk) * rkmh(i,j)
         salt(i,j,kk,2) = saltin(kk) * rkmh(i,j)
         uvel(i,j,kk,2) = uvelin(kk) * rkmh(i,j)
         vvel(i,j,kk,2) = vvelin(kk) * rkmh(i,j)
         relax_South = 0.0
         relax_East = 0.0
#include "cparam.h"
         if (relax_South .eq. 1) then
         if (j .le. widS) temp(i,j,kk,2) = temp(i,j,kk,2) -
     & (tempin(kk) - temp_read(i,j,kk,month))*dtts*
     & rkmh(i,j)/(day2sec*100.0*((j-1)/widS))
         if (j .le. widS) salt(i,j,kk,2) = salt(i,j,kk,2) -
     & (saltin(kk) - salt_read(i,j,kk,month))*dtts*
     & rkmh(i,j)/(day2sec*100.0*((j-1)/widS))
c         if (j .le. widS) pvort(i,1,j) = (j-1)/widS
c         if (j .le. widS) pvort(i,2,j) = (widS+1-j)/widS
         endif
         if (relax_East .eq. 1) then
         if (i .ge. imt-widE+1) temp(i,j,kk,2) = temp(i,j,kk,2) -
     & (tempin(kk) - temp_read(i,j,kk,month))*dtts*
     & rkmh(i,j)/(day2sec*100.0*((widE+i-imt)/widE))
         if (i .ge. imt-widE+1) salt(i,j,kk,2) = salt(i,j,kk,2) -
     & (saltin(kk) - salt_read(i,j,kk,month))*dtts*
     & rkmh(i,j)/(day2sec*100.0*((widE+i-imt)/widE))
c         if (i .ge. imt-widE+1) pvort(i,1,j) = (imt-i)/widE
c         if (i .ge. imt-widE+1) pvort(i,2,j) = (widE+i-imt)/widE
         endif
         enddo

        if ( mod(loop,100000) .eq. 0) then
         do kk = 3,kmax-2
         temp(i,j,kk,2) = (1.0/16.0)*
     &   (-tempin(kk-2) +4.0*tempin(kk-1) + 10.0*tempin(kk)
     &    +4.0*tempin(kk+1) - tempin(kk+2) ) * rkmh(i,j)
         salt(i,j,kk,2) = (1.0/16.0)*
     &   (-saltin(kk-2) +4.0*saltin(kk-1) + 10.0*saltin(kk)
     &    +4.0*saltin(kk+1) - saltin(kk+2) ) * rkmh(i,j)
         uvel(i,j,kk,2) = (1.0/16.0)*
     &   (-uvelin(kk-2) +4.0*uvelin(kk-1) + 10.0*uvelin(kk)
     &    +4.0*uvelin(kk+1) - uvelin(kk+2) ) * rkmh(i,j)
         vvel(i,j,kk,2) = (1.0/16.0)*
     &   (-vvelin(kk-2) +4.0*vvelin(kk-1) + 10.0*vvelin(kk)
     &    +4.0*vvelin(kk+1) - vvelin(kk+2) ) * rkmh(i,j)
         enddo
         endif
         do kk = 1,kmax
         SMCoeff(i,j,kk) = entr_flux(kk) * rkmh(i,j)
         SHCoeff(i,j,kk) = SH(kk) * rkmh(i,j)
         enddo
#ifdef snap_mld_extended
         do kk = 1,kmax
         diag_ext1(i,j,kk) = diag1(kk) * rkmh(i,j)
         diag_ext2(i,j,kk) = diag2(kk) * rkmh(i,j)
         diag_ext3(i,j,kk) = diag3(kk) * rkmh(i,j)
         diag_ext4(i,j,kk) = diag4(kk) * rkmh(i,j)
         diag_ext5(i,j,kk) = diag5(kk) * rkmh(i,j)
         diag_ext6(i,j,kk) = diag6(kk) * rkmh(i,j)
         enddo
#endif

         mld_mld(i,j) = icval

c          call mld_consistency (kmax,ii,jj)
         endif

        enddo
        enddo
c            write(*,*) SHCoeff(imt/2-10:imt/2+10,1,jmt/2+20)
c            write(*,*) rmld_misc(imt/2:imt/2,12,jmt/2+20)

        return
        end
!! Here is the end of Mixed Layer Model
