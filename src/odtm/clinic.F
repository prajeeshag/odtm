	subroutine clinic
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	subroutine to solve baroclinic pressure gradient
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : dxu, dyu, dz, dz_max, gdx, gdy, h, i, j, loop
        use size_mod, only : rdx, rkmh, rdy, taum1, taun1, taup1, u, v, wd
        use size_mod, only : we, we_mld, dyv, k, taum, taun, taup
        use param_mod, only : day2sec, diffuse_h, diffuse_th, dphi, dt, dthe
        use param_mod, only : dtts, pi, re, theta, tropdubdx, tropdvbdy

        use advec_mod, only : sum_2pt
        implicit none
       real tx,ty, Lv
        integer :: ip, im, jp, jm,ii, jj
        real :: rtemp1, rtemp2, rtheta_vu, rtheta_vd, rtemp3, rtemp4, rmask_right
        real :: rmask_left, rmask_top,rmask_bot, rthetapone, rthetamone, rtemp5
        real :: rtemp6, x0, y0, A0, rexp, rsig, rsig1, dtts_back 
        logical :: NaN_check
	taum1 = 1
	taun1 = 2
	taup1 = 3
	
	ip = i + 1
	im = i - 1
	jp = j + 1
	jm = j - 1

	tropdubdx = 0.0
	tropdvbdy = 0.0
	
	rtemp1 = u(ip,k,j,taun)* sum_2pt (3,ip,j,i,j) /
     &                           sum_2pt (31,ip,j,i,j)
	rtemp1 = rtemp1 
	rtemp2 = u(i,k,j,taun)*  sum_2pt (3,im,j,i,j) /
     &                           sum_2pt (31,im,j,i,j)
	rtemp2 = rtemp2 
	dphi = (rdx(i+1)+rdx(i))/2.0

	tropdubdx = (rtemp1 - rtemp2)/
     &              ( (dxu(i,j)+dxu(ip,j))/2  )


	rtheta_vu = cos((theta(j)+theta(j+1))*0.5)
        rtheta_vd = cos((theta(j)+theta(j-1))*0.5)
	dthe = (rdy(j) + rdy(j+1))/2.0

	rtemp3 = v(i,k,jp,taun) * sum_2pt (3,i,j,i,jp) *rtheta_vu /
     &                            sum_2pt (31,i,j,i,jp)
	rtemp3 = rtemp3 
	rtemp4 = v(i,k,j,taun) * sum_2pt (3,i,j,i,jm)  *rtheta_vd/
     &                           sum_2pt (31,i,j,i,jm)
	rtemp4 = rtemp4

	tropdvbdy = ( rtemp3 - rtemp4 )/
     &              (  cos(theta(j)) *
     &                (dyv(i,jp)+dyv(i,j))/2 )


	
#ifdef restrt
        if ( loop .eq. 1) dtts = dtts
#else
        if ( loop .eq. 1) dtts = dt
#endif

	rmask_right = 1.0
	rmask_left = 1.0
	rmask_top = 1.0
	rmask_bot = 1.0
	
	rthetapone = (theta(j+1)+theta(j) )*0.5
	rthetamone = (theta(j-1)+theta(j) )*0.5
	rtemp1 = Re*Re*cos(theta(j))*cos(theta(j)) 
	rtemp1 = 1.0/rtemp1
	 if (rkmh(i+1,j) .eq. 0.0) rmask_right = 0.0
	 if (rkmh(i-1,j) .eq. 0.0) rmask_left = 0.0
	 if (rkmh(i,j+1) .eq. 0.0) rmask_top = 0.0
	 if (rkmh(i,j-1) .eq. 0.0) rmask_bot = 0.0
 
	dphi = (rdx(i+1)+rdx(i))/2.0
	rtemp2 = rmask_right*(h(i+1,k,j,taum) - h(i,k,j,taum) )/dphi
	dphi = (rdx(i)+rdx(i-1))/2.0
	rtemp3 = rmask_left*(h(i,k,j,taum) - h(i-1,k,j,taum) )/dphi
	dphi = rdx(i)
	rtemp4 = (rtemp2 - rtemp3)/dphi
	rtemp5 = rtemp1*rtemp4  ! reserve
	rtemp1 = Re*Re*cos(theta(j))
	rtemp1 = 1.0/rtemp1
	dthe = (rdy(j+1) + rdy(j))/2.0
	rtemp2 = rmask_top*(h(i,k,j+1,taun) - h(i,k,j,taun))/dthe
	dthe = (rdy(j-1) + rdy(j))/2.0
	rtemp3 = rmask_bot*(h(i,k,j,taun) - h(i,k,j-1,taun))/dthe
	dthe = rdy(j)
	rtemp4 = (cos(rthetapone)*rtemp2 - cos(rthetamone)*rtemp3)/dthe
	rtemp6 = rtemp1*rtemp4  ! reserve
	
	diffuse_h = rtemp5 + rtemp6
	
	 x0 = 135.0
	 y0 = 45.0
	 A0 = 25.0e0
	 rexp = ((i-x0)/5.0)**2 + ((j-y0)/5.0)**2
	 rsig = 2.0*pi*loop/((day2sec*365/dt))
	 rsig1 = sin(rsig)
c	 rexp1 = dz + A0*exp(-rexp)*rsig1
	
	 dtts_back = dtts
c         if (j .lt. 20) then
c	 dtts = dtts*(j/20.0)
c	 diffuse_th = diffuse_th*20/(j*1.0)
c         endif
	
	rtemp1 = (gdx(i+1) - gdx(i))
        rtemp2 = (gdy(j+1) - gdy(j))
        
c        if (rkmt(i,j) * rkmh(i,j) * rkmu(i,j) * rkmv(i,j) .eq. 0.0) then
c        rtemp = 0.0
c        else
c        rtemp = 0.0
c        endif

        we_mld(i,k,j) = (tropdubdx + tropdvbdy) + we(i,k,j) 
c        pvort(i,k,j) =  we_mld(i,k,j)

        h(i,k,j,taup) = h(i,k,j,taum)*rkmh(i,j) + (
     &  -(
     &     tropdubdx 
     &     +  tropdvbdy  
     &      )
#ifdef entrain
     &     +  we(i,k,j)*rkmh(i,j)*dtts_back/dtts
     &     +  wd(i,k,j)*rkmh(i,j)*dtts_back/dtts
#endif
     &     + diffuse_th * diffuse_h
     &  ) *dtts*rkmh(i,j)

	
#ifdef atmosphere
	 call couple
#endif
        
        if ( h(i,k,j,taup) .gt. dz_max(k)) then
        do ii=im,ip
        do jj=jm,jp
         h(ii,k,jj,taup) =  h(ii,k,jj,taup) - (10**(-1.5))*  h(ii,k,jj,taun)*
     &   rkmh(ii,jj)
c        pvort(ii,k,jj) = 1.0
        enddo
        enddo
        endif

        NaN_check  = isnan(h(i,k,j,taup))
        
        if (NaN_check .eqv. .true.) h(i,k,j,taup) = dz(k)*rkmh(i,j)
	return
	end
