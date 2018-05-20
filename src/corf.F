	subroutine corf
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	find coriolis component
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
        use advec_mod, only : sum_4pt, sum_2pt
        implicit none
        integer :: ip, im, jp, jm
        real :: fcoru, rtemp1, rtemp2, fcorv 
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1
	
#ifdef area_weighted123
	do ii = im,ip
        do jj = jm,jp
         u(ii,k,jj,taun) = u(ii,k,jj,taun)*dxu(ii,jj)*dyv(ii,jj)
         v(ii,k,jj,taun) = v(ii,k,jj,taun)*dyv(ii,jj)*dxu(ii,jj)
         rkmu(ii,jj) = rkmu(ii,jj) * dxu(ii,jj)*dyv(ii,jj)
         rkmv(ii,jj) = rkmv(ii,jj) * dyv(ii,jj)*dxu(ii,jj)
        enddo
        enddo
#endif
	
	fcoru =  (fcor(im,j)+fcor(i,j))*0.5
	rtemp1 = sum_4pt (2,im,jp,i,jp,i,j,im,j) /
     &           sum_4pt (21,im,jp,i,jp,i,j,im,j)
	rtemp2 = sum_2pt (3,i,j,im,j) /sum_2pt (31,i,j,im,j)
	fcoruh = fcoru * rtemp1 * rtemp2  


	fcorv = (fcor(i,j) + fcor(i,jm))*0.5
	rtemp1 = sum_4pt (1,i,j,ip,j,ip,jm,i,jm) /
     &           sum_4pt (11,i,j,ip,j,ip,jm,i,jm)
	rtemp2 = sum_2pt (3,i,j,i,jm) /sum_2pt (31,i,j,i,jm)
	
	
	fcorvh = fcorv * rtemp1 * rtemp2  
	

#ifdef area_weighted123
	do ii = im,ip
        do jj = jm,jp
         u(ii,k,jj,taun) = u(ii,k,jj,taun)/dxu(ii,jj)/dyv(ii,jj)
         v(ii,k,jj,taun) = v(ii,k,jj,taun)/dyv(ii,jj)/dxu(ii,jj)
         rkmu(ii,jj) = rkmu(ii,jj) / dxu(ii,jj)/dyv(ii,jj)
         rkmv(ii,jj) = rkmv(ii,jj) / dyv(ii,jj)/dxu(ii,jj)
        enddo
        enddo
#endif

	return
	stop
	end
