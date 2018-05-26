	subroutine polar_coord
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c  attempt to impliment polar coordinate
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : imt, jmt, i, j, k, gdx, gdy, fcor, gdxb, gdyb
        use size_mod, only : dxu, dyu, dxv, dyv, dxh, dyh, dah
        use size_mod, only : rdxu,rdyu, rdxv, rdyv, rdxh, rdyh
        use size_mod, only : phi, theta, tan_theta, fcoru, fcorv
        use size_mod, only : isc, iec, jsc, jec
        use param_mod, only : deg2rad, omega, pi, re, deg2meter
        implicit none

        integer :: ip, im, jp, jm
        real :: deg2meter_x

	do i=1,imt
	 phi(i) = gdx(i)*deg2rad  !converting latitude to radian
	enddo
	
	
	do j=1,jmt
	 theta(j) = gdy(j)*deg2rad      !converting longitude to radian
	enddo
	
c
c define coriolis parameter f=2*omega*sin(theta)
c
	do i=isc, iec
	    do j=jsc, jec
	        fcor(i,j) = 2.0*omega*sin(theta(j))
            fcoru(i,j) = 2.0*omega*sin(gdy(j)*deg2rad)
            fcorv(i,j) = 2.0*omega*sin(gdyb(j)*deg2rad) 
	    enddo
	enddo
	
      dxu(:,:) = 0.
      dyu(:,:) = 0.
      dxv(:,:) = 0.
      dyv(:,:) = 0.
      dxh(:,:) = 0.
      dyh(:,:) = 0.

      rdxu(:,:) = 0.
      rdyu(:,:) = 0.
      rdxv(:,:) = 0.
      rdyv(:,:) = 0.
      rdxh(:,:) = 0.
      rdyh(:,:) = 0.
    
	do i=1,imt
	    do j=1,jmt
	        ip = i + 1
	        im = i - 1 
	        jp = j + 1
	        jm = j - 1
            deg2meter_x = deg2meter * cos(theta(j))
            if(i>1) dxu(i,j) = deg2meter_x * (gdx(i) - gdx(im))
            if (dxu(i,j)>0.) rdxu(i,j) = 1./dxu(i,j) 

            dyu(i,j) = deg2meter * (gdyb(jp) - gdyb(j))
            if (dyu(i,j)>0.) rdyu(i,j) = 1./dyu(i,j) 

            dxv(i,j) = deg2meter_x * (gdxb(ip) - gdxb(i))
            if (dxv(i,j)>0.) rdxv(i,j) = 1./dxv(i,j) 

            if(j>1) dyv(i,j) = deg2meter * (gdy(j) - gdy(jm))
            if (dyv(i,j)>0.) rdyv(i,j) = 1./dyv(i,j) 

            dxh(i,j) = deg2meter_x * (gdxb(ip) - gdxb(i))
            if (dxh(i,j)>0.) rdxh(i,j) = 1./dxh(i,j) 

            dyh(i,j) = deg2meter * (gdyb(jp) - gdyb(j))
            if (dyh(i,j)>0.) rdyh(i,j) = 1./dyh(i,j) 

	        dah(i,j) = (dxh(i,j)*dyh(i,j))
	    enddo
	enddo

        tan_theta(:) = 0.
        tan_theta(1:jmt) = tan(theta(1:jmt))
	
	return
	end subroutine polar_coord
