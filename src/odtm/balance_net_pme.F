         subroutine balance_pme
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       subroutine to solve baroclinic pressure gradient
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : pme_corr, i, imt, j, jmt, rmld_misc
        implicit none

        integer :: ievap, iprecip
        real :: revap

        pme_corr(:,:) = 0.0

        ievap = 0
        iprecip = 0
        revap = 0.0
        do i=1,imt
        do j=1,jmt
         revap = revap + rmld_misc(i,j,3)
         if (rmld_misc(i,j,3) .gt. 0.0) ievap = ievap + 1
         if (rmld_misc(i,j,3) .lt. 0.0) iprecip = iprecip + 1
        enddo
        enddo

        do i=1,imt
        do j=1,jmt
         if (revap .gt. 0.0) then
         if (rmld_misc(i,j,3) .gt. 0.0) pme_corr(i,j) =  1.0*revap/ievap
         endif
         if (revap .lt. 0.0) then
         if (rmld_misc(i,j,3) .lt. 0.0) pme_corr(i,j) =  1.0*revap/iprecip
         endif
        enddo
        enddo

        return
        end subroutine balance_pme
