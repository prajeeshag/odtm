       
        module interp_extrap_initial_mod

        contains
 
        subroutine interp_extrap_initial (ii,jj,kk,kmax,tempin,saltin,
     &  tempout,saltout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c       linear interpolation for wind forcing
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : h, taun, t, u, v
        implicit none

        real :: tempin(kmax), saltin(kmax), tempout, saltout
        real :: rsum1, rsum2
        integer :: hstart, hend, ii, jj, kk, tau_index, kmax
        integer :: irepeat, rsum3, rsum4, icount, kkk

        hstart = 0
        do irepeat = 1,kk-1
        hstart = hstart + h(ii,jj,irepeat,taun)
        enddo

        hend = h(ii,jj,kk,taun) + hstart
        rsum1 = 0.0
        rsum2 = 0.0
        rsum3 = 0.0
        rsum4 = 0.0
        icount = 0

        tempout = 0.0
        saltout = 0.0

        do kkk = 1+hstart/5, hend/5
        icount = icount + 1
        if (kkk .gt. kmax-1) goto 101
        rsum1 = rsum1 + tempin(kkk) 
        rsum2 = rsum2 + saltin(kkk)
c        if (kkk .gt. kmax) goto 101
        enddo

        tempout = rsum1/max(1,icount)
        saltout = rsum2/max(1,icount)
 101     continue  


      if (hstart .gt. kmax*5.0) goto 102
      if (kkk .gt. kmax) then

        tempout = (t(ii,jj,kk,1,taun) * (h(ii,jj,kk,taun)-icount*5.0) +
     &   rsum1*5.0 )/max(1.0,h(ii,jj,kk,taun))
        saltout = (t(ii,jj,kk,2,taun) * (h(ii,jj,kk,taun)-icount*5.0) +
     &   rsum2*5.0 )/max(1.0,h(ii,jj,kk,taun))

        endif


 102    continue
        return
        end subroutine interp_extrap_initial

        end module interp_extrap_initial_mod 
