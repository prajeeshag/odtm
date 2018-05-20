	subroutine smooth (idim1, idim2, idim3, idim4, idim5,
     &  imts, jmts, kmts, imte, jmte, kmte)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	subroutine to smooth out un-wanted noise in the solutions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : dxu, dxv, dyu, dyv, dz, h, i, j, k, rkmh, rkmhh
        use size_mod, only : rkmu, rkmv, t, taun1, u, v
        implicit none
        
        integer :: idim1, idim2, idim3, idim4, idim5, imts, jmts, kmts, imte
        integer :: jmte, kmte, ii, im, ip, jj, jm, jp, kk, kmone, kp, nt
        real :: rdiv1, rdiv2, ru, rv, rh, rt
c=======================================================================
        ip = imte
        im = imts
	i = (ip+im)/2
        jp = jmte
        jm = jmts
	j = (jp+jm)/2
        kmone = kmts
        kp = kmte
        k = (kmone + kp)/2

#ifdef area_weighted
        do ii = im,ip
        do jj = jm,jp
	do kk=kmts, kmte
         u(ii,kk,jj,4) = u(ii,kk,jj,4)*dxu(ii,jj)*dyu(ii,jj)
         v(ii,kk,jj,4) = v(ii,kk,jj,4)*dyv(ii,jj)*dxv(ii,jj)
         u(ii,kk,jj,3) = u(ii,kk,jj,3)*dxu(ii,jj)*dyu(ii,jj)
         v(ii,kk,jj,3) = v(ii,kk,jj,3)*dyv(ii,jj)*dxv(ii,jj)
         if (h(ii,kk,jj,2) .eq. 0.0) h(ii,kk,jj,4) = dz(kk)
         if (h(ii,kk,jj,2) .eq. 0.0) h(ii,kk,jj,3) = dz(kk)
        enddo
         rkmu(ii,jj) = rkmu(ii,jj) * dxu(ii,jj)*dyu(ii,jj)
         rkmv(ii,jj) = rkmv(ii,jj) * dyv(ii,jj)*dxv(ii,jj)
        enddo
        enddo
#endif


	do k=kmts, kmte
	do taun1 = 4,3,-1

	if ( idim1 .eq. 1 .and. idim5 .eq. 5) then

        rdiv1 =   (rkmu(im,j) + rkmu(ip,j) +
     &             rkmu(i,jm) + rkmu(i,jp))
        rdiv2 =   (rkmu(im,jm) + rkmu(ip,jm) +
     &             rkmu(im,jp) + rkmu(ip,jp))
        rU = 0.25 * u(i,k,j,taun1)  + 
     &       0.5* (u(im,k,j,taun1) + u(ip,k,j,taun1) +
     &               u(i,k,jm,taun1) + u(i,k,jp,taun1))/max(1.0,rdiv1)+
     &       0.25* (u(im,k,jm,taun1) + u(ip,k,jm,taun1) +
     &                u(im,k,jp,taun1) + u(ip,k,jp,taun1))/max(1.0,rdiv2)
        u(i,k,j,taun1-2) = rU / dxu(i,j)/dyu(i,j)

        endif

	
	if ( idim2 .eq. 2 .and. idim5 .eq. 5) then

        rdiv1 =   (rkmv(im,j) + rkmv(ip,j) +
     &             rkmv(i,jm) + rkmv(i,jp))
        rdiv2 =   (rkmv(im,jm) + rkmv(ip,jm) +
     &             rkmv(im,jp) + rkmv(ip,jp))
        rV = 0.25 * v(i,k,j,taun1)  + 
     &       0.5* (v(im,k,j,taun1) + v(ip,k,j,taun1) +
     &               v(i,k,jm,taun1) + v(i,k,jp,taun1))/max(1.0,rdiv1)+
     &       0.25* (v(im,k,jm,taun1) + v(ip,k,jm,taun1) +
     &                v(im,k,jp,taun1) + v(ip,k,jp,taun1))/max(1.0,rdiv2)
        v(i,k,j,taun1-2) = rV / dyv(i,j)/dxv(i,j) 

	endif
	
	if ( idim3 .eq. 3 .and. idim5 .eq. 5) then

        rdiv1 =   (rkmhH(im,j) + rkmhH(ip,j) +
     &             rkmhH(i,jm) + rkmhH(i,jp))
        rdiv2 =   (rkmhH(im,jm) + rkmhH(ip,jm) +
     &             rkmhH(im,jp) + rkmhH(ip,jp))
        rH = 0.25 * h(i,k,j,taun1)  + 
     &       0.5* (h(im,k,j,taun1) + h(ip,k,j,taun1) +
     &               h(i,k,jm,taun1) + h(i,k,jp,taun1))/max(1.0,rdiv1)+
     &       0.25* (h(im,k,jm,taun1) + h(ip,k,jm,taun1) +
     &                h(im,k,jp,taun1) + h(ip,k,jp,taun1))/max(1.0,rdiv2)
        h(i,k,j,taun1-2) = rH * rkmhH(i,j)

        

	endif
        
	if ( idim3 .eq. 4 .and. idim5 .eq. 5) then
        nt = 1
        rdiv1 =   (rkmh(im,j) + rkmh(ip,j) +
     &             rkmh(i,jm) + rkmh(i,jp))
        rdiv2 =   (rkmh(im,jm) + rkmh(ip,jm) +
     &             rkmh(im,jp) + rkmh(ip,jp))
        rT = 0.25 * t(i,k,j,nt,taun1)  + 
     &       0.5* (t(im,k,j,nt,taun1) + t(ip,k,j,nt,taun1) +
     &        t(i,k,jm,nt,taun1) + t(i,k,jp,nt,taun1))/max(1.0,rdiv1)+
     &       0.25* (t(im,k,jm,nt,taun1) + t(ip,k,jm,nt,taun1) +
     &        t(im,k,jp,nt,taun1) + t(ip,k,jp,nt,taun1))/max(1.0,rdiv2)
        if (rdiv1 .eq. 4 .and. rdiv2 .eq. 4) t(i,k,j,nt,taun1-2) = rT

	endif


	enddo
	enddo
	
#ifdef area_weighted
        do ii = im,ip
        do jj = jm,jp
	do kk=kmts, kmte
         u(ii,kk,jj,4) = u(ii,kk,jj,4)/dxu(ii,jj)/dyu(ii,jj)
         v(ii,kk,jj,4) = v(ii,kk,jj,4)/dyv(ii,jj)/dxv(ii,jj)
         u(ii,kk,jj,3) = u(ii,kk,jj,3)/dxu(ii,jj)/dyu(ii,jj)
         v(ii,kk,jj,3) = v(ii,kk,jj,3)/dyv(ii,jj)/dxv(ii,jj)
         if (h(ii,kk,jj,2) .eq. 0.0) h(ii,kk,jj,4) = 0.0
         if (h(ii,kk,jj,2) .eq. 0.0) h(ii,kk,jj,3) = 0.0
        enddo
         rkmu(ii,jj) = rkmu(ii,jj) / dxu(ii,jj)/dyu(ii,jj)
         rkmv(ii,jj) = rkmv(ii,jj) / dyv(ii,jj)/dxv(ii,jj)
        enddo
        enddo
#endif
	return
	stop
	end
