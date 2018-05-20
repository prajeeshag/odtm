	subroutine entrain_detrain
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	attempt to impliment artificial entrainment term
c	McCreary and Lu (1994), J.P.O, Vol.24, pp. 466-497
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : we, wd, h, hd, he, i, j, k, rkmh, taun
            use size_mod, only : imt, jmt, km
            use param_mod, only : td, te
        implicit none
        real :: rwesum
	rwesum = 0.0
	we(:,:,:) = 0.0
	wd(:,:,:) = 0.0
c
c entrain from layer below to layer above if mass receds he(k)
c
	
	do k=1,km-1
	do i=2,imt-1
	do j=2,jmt-1
	
	if (he(k) .gt. h(i,k,j,taun) )  then
	 we(i,k,j) =  rkmh(i,j) *( he(k) - h(i,k,j,taun) )**2
     &           /(he(k)*te)
	 wd(i,k+1,j) =  -1.0*rkmh(i,j) *( he(k) - h(i,k,j,taun) )**2
     &           /(he(k)*te)
	
         wd(i,km,j) = 0.0
	endif

	enddo
	enddo
	enddo
	
c
c detrain from layer above to below if mass exceeds hd(k)
c
	do k=1,km-1
        do i=2,imt-1
        do j=2,jmt-1

        if (hd(k) .lt. h(i,k,j,taun) )  then
         wd(i,k,j) =  wd(i,k,j) - rkmh(i,j)*(hd(k)-h(i,k,j,taun))**2
     &           /(hd(k)*td)
         we(i,k+1,j) =  we(i,k+1,j) + rkmh(i,j)*(hd(k)-h(i,k,j,taun))**2
     &           /(hd(k)*td)

        endif

        enddo
        enddo
        enddo

	return
	stop
	end
