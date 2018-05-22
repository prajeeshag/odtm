        module advec_mod
        implicit none
        private
        
        public :: advec, sum_2pt, sum_4pt

        contains 

        subroutine advec
        ! implimenting fully non-linear advection

        use size_mod, only : i, j, k, taun, u, v, h, rkmu, rkmv, rdx, rdy, dxu, dxv
        use size_mod, only : dyu, dyv
        use param_mod, only : advec_uu, advec_uv, advec_vu, advec_vv, dphi
        use param_mod, only : dthe, polar_u, polar_v, re, theta

        implicit none

        real :: rtemp1, rtemp2, rtemp3, rtemp4, rtemp5, rtemp6, rtemp7, rtemp8
        real :: rUp, rDw, rtheta_temp
        integer :: ip, im, jp, jm, jj, ii
   
	rtemp1=0.0
	rtemp2=0.0
	rtemp3=0.0
	rtemp4=0.0
	rtemp5=0.0
	rtemp6=0.0
	rtemp7=0.0
	rtemp8=0.0
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1

#ifdef area_weighted123
	do ii = im,ip
	do jj = jm,jp
	 u(ii,k,jj,taun) = u(ii,k,jj,taun)*dxu(ii,jj)*dyu(ii,jj)
	 v(ii,k,jj,taun) = v(ii,k,jj,taun)*dyv(ii,jj)*dxv(ii,jj)
	 rkmu(ii,jj) = rkmu(ii,jj) * dxu(ii,jj)*dyu(ii,jj)
	 rkmv(ii,jj) = rkmv(ii,jj) * dyv(ii,jj)*dxv(ii,jj)
	enddo
	enddo
#endif
	
!cccccccccccccccc u.du/dx + v.du/dy cccccccccccccccccccccccccccccccccccccccc
	dphi = (rdx(i-1) + rdx(i))/2.0
	dthe = rdy(j)


	rtemp1 = sum_2pt (1,i,j,ip,j) / sum_2pt(11,i,j,ip,j)
	rtemp2 = rtemp1 * rtemp1 * h(i,k,j,taun)
	rtemp3 = sum_2pt (1,i,j,im,j) / sum_2pt(11,i,j,im,j)
	rtemp4 = rtemp3 * rtemp3 * h(im,k,j,taun)
	
	advec_uu = (rtemp2 - rtemp4)/ dxu(i,j)
        
	
	rtemp1 = sum_2pt (1,i,jp,i,j) /sum_2pt (11,i,jp,i,j)
	rtemp2 = sum_2pt (2,im,jp,i,jp) /sum_2pt (21,im,jp,i,jp)
	rtemp3 = sum_4pt(3,im,jp,i,jp,i,j,im,j)/
     &           sum_4pt(31,im,jp,i,jp,i,j,im,j)
	rtemp4 = rtemp1 * rtemp2 * rtemp3

	rtemp5 = sum_2pt (1,i,j,i,jm) /sum_2pt (11,i,j,i,jm)
	rtemp6 = sum_2pt (2,i,j,im,j) /sum_2pt (21,i,j,im,j)
	rtemp7 = sum_4pt (3,im,j,i,j,i,jm,im,jm)/
     &           sum_4pt (31,im,j,i,j,i,jm,im,jm)
	rtemp8 = rtemp5 * rtemp6 * rtemp7

	rUp = (dxv(im,jp) + dxv(i,jp))*0.5
	rDw = (dxv(im,j) + dxv(i,j))*0.5
	advec_vu = (rtemp4*rUp - rtemp8*rDw)/ (dyu(i,j)*dxu(i,j))

	
	rtemp1 = u(i,k,j,taun) / dxu(i,j) / dyu(i,j)
	rtemp2 = sum_4pt (2,i,j,i,jp,im,jp,im,j)/
     &           sum_4pt (21,i,j,i,jp,im,jp,im,j)
	rtemp3 = sum_2pt (3,i,j,im,j) /sum_2pt (31,i,j,im,j)

	polar_u = 2.0*(rtemp1*rtemp2*rtemp3)/(Re*(1/tan(theta(j))))


	rtemp1=0.0
	rtemp2=0.0
	rtemp3=0.0
	rtemp4=0.0
	rtemp5=0.0
	rtemp6=0.0
	rtemp7=0.0
	rtemp8=0.0
!cccccccccccccccc u.dv/dx + v.dv/dy cccccccccccccccccccccccccccccccccccccccc
	dphi = rdx(i)
	dthe = (rdy(j-1) + rdy(j))/2.0
	 rtheta_temp = (theta(j) + theta(j-1))*0.5


	rtemp1 = sum_2pt (2,i,jp,i,j) /sum_2pt (21,i,jp,i,j)
	rtemp2 = rtemp1 * rtemp1 * h(i,k,j,taun)
	rtemp3 = sum_2pt (2,i,j,i,jm) /sum_2pt (21,i,j,i,jm)
	rtemp4 = rtemp3 * rtemp3 * h(i,k,jm,taun)

	rUp = (dxv(i,jp) + dxv(i,j))*0.5
	rDw = (dxv(i,j) + dxv(i,jm))*0.5
	advec_vv = (rtemp2*rUp - rtemp4*rDw)/ (dyv(i,j)*dxv(i,j))


	rtemp1 = sum_2pt (1,ip,jm,ip,j) /sum_2pt (11,ip,jm,ip,j)
	rtemp2 = sum_2pt (2,i,j,ip,j) /sum_2pt (21,i,j,ip,j)
	rtemp3 = sum_4pt (3,i,j,ip,j,ip,jm,i,jm)/
     &           sum_4pt (31,i,j,ip,j,ip,jm,i,jm)    
	rtemp4 = rtemp1 * rtemp2 * rtemp3
	rtemp5 = sum_2pt (1,i,j,i,jm) /sum_2pt (11,i,j,i,jm)
	rtemp6 = sum_2pt (2,i,j,im,j) /sum_2pt (21,i,j,im,j)
	rtemp7 = sum_4pt (3,im,j,i,j,i,jm,im,jm)/
     &           sum_4pt (31,im,j,i,j,i,jm,im,jm)
	rtemp8 = rtemp5 * rtemp6 * rtemp7
	
	advec_uv = (rtemp4 - rtemp8)/ dxv(i,j)


	rtemp1 = sum_4pt (1,i,j,ip,j,ip,jm,i,jm)/
     &           sum_4pt (11,i,j,ip,j,ip,jm,i,jm)
	rtemp2 = v(i,k,j,taun) / dyv(i,j) / dxv(i,j)
	rtemp3 = sum_2pt (3,i,j,i,jm) /sum_2pt (31,i,j,i,jm)
	rtemp4 = (rtemp1*rtemp1 - rtemp2*rtemp2) * rtemp3

	polar_v = rtemp4/(Re*(1/tan(theta(j))))
	
#ifdef area_weighted123
	do ii = im,ip
	do jj = jm,jp
	 u(ii,k,jj,taun) = u(ii,k,jj,taun)/dxu(ii,jj)/dyu(ii,jj)
	 v(ii,k,jj,taun) = v(ii,k,jj,taun)/dyv(ii,jj)/dxv(ii,jj)
	 rkmu(ii,jj) = rkmu(ii,jj) / dxu(ii,jj)/dyu(ii,jj)
	 rkmv(ii,jj) = rkmv(ii,jj) / dyv(ii,jj)/dxv(ii,jj)
	enddo
	enddo
#endif
	
	
	return
	stop
	end subroutine advec


	function sum_4pt (index, i1,j1, i2,j2, i3,j3, i4,j4)
        use size_mod, only : u, v, h, taun, k, rkmu, rkmv, rkmh
        implicit none
        
        real :: sum_4pt
        integer, intent(in) :: index, i1,j1, i2,j2, i3,j3, i4,j4
	
	if (index .eq. 1) then
	sum_4pt = (u(i1,k,j1,taun) + u(i2,k,j2,taun) + 
     &             u(i3,k,j3,taun) + u(i4,k,j4,taun))
	endif

	if (index .eq. 2) then
	sum_4pt = (v(i1,k,j1,taun) + v(i2,k,j2,taun) + 
     &             v(i3,k,j3,taun) + v(i4,k,j4,taun))
	endif

	if (index .eq. 3) then
	sum_4pt = (h(i1,k,j1,taun) + h(i2,k,j2,taun) + 
     &             h(i3,k,j3,taun) + h(i4,k,j4,taun))
	endif

	if (index .eq. 11) then
	sum_4pt = (rkmu(i1,j1) + rkmu(i2,j2) + rkmu(i3,j3) + rkmu(i4,j4))
	sum_4pt = max(1.0,sum_4pt)
	endif

	if (index .eq. 21) then
	sum_4pt = (rkmv(i1,j1) + rkmv(i2,j2) + rkmv(i3,j3) + rkmv(i4,j4))
	sum_4pt = max(1.0,sum_4pt)
	endif

	if (index .eq. 31) then
	sum_4pt = (rkmh(i1,j1) + rkmh(i2,j2) + rkmh(i3,j3) + rkmh(i4,j4))
	sum_4pt = max(1.0,sum_4pt)
	endif

	end function sum_4pt


	function sum_2pt (index, i1,j1, i2,j2)
        use size_mod, only: k, taun, u, v, h, rkmu, rkmv, rkmh
        implicit none
        
        real :: sum_2pt
        integer, intent(in) :: index, i1,j1, i2,j2

	
	if (index .eq. 1) then
	sum_2pt = (u(i1,k,j1,taun) + u(i2,k,j2,taun))
	endif

	if (index .eq. 2) then
	sum_2pt = (v(i1,k,j1,taun) + v(i2,k,j2,taun))
	endif

	if (index .eq. 3) then
	sum_2pt = (h(i1,k,j1,taun) + h(i2,k,j2,taun))
	endif

	if (index .eq. 11) then
	sum_2pt = (rkmu(i1,j1) + rkmu(i2,j2))
	sum_2pt = max(1.0,sum_2pt)
	endif

	if (index .eq. 21) then
	sum_2pt = (rkmv(i1,j1) + rkmv(i2,j2))
	sum_2pt = max(1.0,sum_2pt)
	endif

	if (index .eq. 31) then
	sum_2pt = (rkmh(i1,j1) + rkmh(i2,j2))
	sum_2pt = max(1.0,sum_2pt)
	endif

	end function sum_2pt

        end module advec_mod
