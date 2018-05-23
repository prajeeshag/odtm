        subroutine interp_extrap (ii,jj,kk,kmax,tau_index,
     &  tempout,saltout,uvelout,vvelout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c       linear interpolation for wind forcing
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : eta, h, loop, salt, t, temp, u, uvel, v
            use size_mod, only : vvel, we_baro, we_mld, we_prof, taun, km
            use size_mod, only : imt, jmt
            use param_mod, only : day2sec, diffuse_my, diffuse_tr, dt
        implicit none
        integer, intent(in) :: ii, jj, kk, kmax, tau_index 

        real, intent(out) :: tempout, saltout
        real, intent(out) ::  uvelout, vvelout
        real rsum1, rsum2, rsum3, rsum4
        integer hstart, hend
        
        real :: alpha_rgm, alpha_rgt, beta_mldm, beta_mldt
        integer :: icount, irepeat, kkk
        real :: rdiv, rdrag, relax_East, relax_South, relax_on_off, rsum, widE, widS
        real :: zero_out_tradv

        hstart = 0
        do irepeat = 1,kk-1
        hstart = hstart + h(ii,irepeat,jj,taun)
        enddo

        hend = h(ii,kk,jj,taun) + hstart
        rsum1 = 0.0
        rsum2 = 0.0
        rsum3 = 0.0
        rsum4 = 0.0
        icount = 0
        
        tempout = 0.0
        saltout = 0.0
        uvelout = 0.0
        vvelout = 0.0

           ! if (ii==imt/2.and.jj==jmt/2.and.tau_index==2) then
           !     print *, 'interp_extrap.F :', hstart, hend, kk, tau_index, kmax, 1+hstart/5, hend/5 
           ! endif

        do kkk = 1+hstart/5, hend/5
        icount = icount + 1

        if (kkk .gt. kmax) goto 101
        rsum1 = rsum1 + temp(ii,kkk,jj,tau_index) 
        rsum2 = rsum2 + salt(ii,kkk,jj,tau_index)
        rsum3 = rsum3 + uvel(ii,kkk,jj,tau_index)
        rsum4 = rsum4 + vvel(ii,kkk,jj,tau_index)
        enddo




        tempout = rsum1/max(1,icount)
        saltout = rsum2/max(1,icount)
        uvelout = rsum3/max(1,icount)
        vvelout = rsum4/max(1,icount)

 101    continue

! if the mixed layer crosses layer 1 of reduced gravity model,
! to handle a portion of layer having reduced gravity, we assume a dummy
! rg layer with constant values same as RG, and this take the average,
! thereby assuming a mixed layer in the complete RG layer 
        if (hstart .gt. kmax*5.0) goto 102
        if (kkk .gt. kmax) then
        tempout = (t(ii,kk,jj,1,tau_index)*
     &    (h(ii,kk,jj,taun)-icount*5.0) + 
     &   rsum1*5.0 )/max(1.0, (h(ii,kk,jj,taun)) )

        saltout = (t(ii,kk,jj,2,tau_index)*
     &    (h(ii,kk,jj,taun)-icount*5.0) + 
     &   rsum2*5.0 )/max(1.0, (h(ii,kk,jj,taun)) )

        uvelout = (u(ii,kk,jj,tau_index)*
     &    (h(ii,kk,jj,taun)-icount*5.0) + 
     &   rsum3*5.0 )/max(1.0, (h(ii,kk,jj,taun)) )

        vvelout = (v(ii,kk,jj,tau_index)*
     &    (h(ii,kk,jj,taun)-icount*5.0) + 
     &   rsum4*5.0 )/max(1.0, (h(ii,kk,jj,taun)) )

        
        endif

        

 102    continue
           ! if (ii==imt/2.and.jj==jmt/2) then
           !     print *, 'interp_extrap.F  2 :', tempout
           ! endif
        return
        end subroutine




        subroutine interp_extrap_we_profile (ii, jj, kmax)
        use size_mod, only : eta, h, loop, we_baro, we_mld, we_prof, taum, taun, km
            use param_mod, only : day2sec, diffuse_my, diffuse_tr, dt
        implicit none
        integer hstart, hend, ii, jj, kk, kmax
        real :: gama, alpha_rgm, alpha_rgt, beta_mldm, beta_mldt
        integer :: icount, irepeat, kkk
        real :: rdiv, rdrag, relax_East, relax_South, relax_on_off, rsum, widE, widS
        real :: zero_out_tradv

        gama = 1.0

#include "cparam.h"  ! <== This file overrides above parameter setup

        we_mld(ii,0,jj) = (eta(ii,1,jj,taun)-eta(ii,1,jj,taum))/dt  ! surface boundary condition
        we_prof(:) = 0.0

        do kk=1,km-1
        hstart = 0
        do irepeat = 1,kk-1
        hstart = hstart + h(ii,irepeat,jj,taun)
        enddo
        hend = h(ii,kk,jj,taun) + hstart

c        write(*,*) hstart, hend

        we_prof(1+hstart/5-1) = we_mld(ii,kk-1,jj)

        rdiv = (1 + hend/5 - hstart/5)
        if (hend/5 .gt. kmax) then
        irepeat = 0
        do kkk = 1+hstart/5, hend/5
        irepeat = irepeat + 1
        if (kkk .ge. kmax) then
        rdiv = irepeat * 1.0
        goto 102
        endif
        enddo
        endif
 102    continue

        do kkk = 1+hstart/5, hend/5
        if (kkk .ge. kmax-1) goto 101
        we_prof(kkk) = we_prof(kkk-1) + 
     &          ((we_mld(ii,kk,jj) - we_mld(ii,kk-1,jj))/
     &          (rdiv)) 
        enddo
        enddo
 101    continue
        rsum = 0.0
        icount = 0
        do kkk=1,kmax
        rsum = rsum + we_prof(kkk)  * gama
        if ( we_prof(kkk) .ne. 0.0) icount = icount + 1
        we_prof(kkk) = we_prof(kkk) * gama

c        if (ii .eq. 28 .and. jj .eq. 52) then
c        write(*,10)  we_prof(kkk), we_mld(ii,0,jj), we_mld(ii,1,jj),
c     &   we_mld(ii,2,jj), we_mld(ii,3,jj), rsum/max(1,icount)
c        endif
c 10     format(e12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5)

        enddo
        we_baro = rsum/max(1,icount)

        return
        end


