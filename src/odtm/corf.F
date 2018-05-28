	subroutine corf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	find coriolis component
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : i, j, fcor, fcoruh, fcorvh, fcoru, fcorv
        use advec_mod, only : sum_4pt, sum_2pt
        implicit none
        integer :: ip, im, jp, jm
        real :: rtemp1, rtemp2
 
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1
	
	rtemp1 = sum_4pt (2,im,jp,i,jp,i,j,im,j) /
     &           sum_4pt (21,im,jp,i,jp,i,j,im,j)
	rtemp2 = sum_2pt (3,i,j,im,j) /sum_2pt (31,i,j,im,j)
	fcoruh = fcoru(i,j) * rtemp1 * rtemp2  


	rtemp1 = sum_4pt (1,i,j,ip,j,ip,jm,i,jm) /
     &           sum_4pt (11,i,j,ip,j,ip,jm,i,jm)
	rtemp2 = sum_2pt (3,i,j,i,jm) /sum_2pt (31,i,j,i,jm)
	
	fcorvh = fcorv(i,j) * rtemp1 * rtemp2  
	

	return
	stop
	end
