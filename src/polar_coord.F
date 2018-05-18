	subroutine polar_coord
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c  attempt to impliment polar coordinate
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"
	do i=1,imt
	 phi(i) = gdx(i)*deg2rad  !converting latitude to radian
	enddo
	
	
	do j=1,jmt
	 theta(j) = gdy(j)*deg2rad      !converting longitude to radian
	enddo
	
c
c define coriolis parameter f=2*omega*sin(theta)
c
	do i=1,imt
	do j=1,jmt
	 fcor(i,j) = 2.0*omega*sin(theta(j))
	enddo
	enddo
	

	gdx(0) = gdx(1)
	gdx(imt+1) = gdx(imt)
	gdy(0) = gdy(1)
	gdy(jmt+1) = gdy(jmt)
	do i=1,imt
	do j=1,jmt
	ip = i + 1
	im = i - 1 
	jp = j + 1
	jm = j - 1 
	dxr(i,j) =((gdx(i)+gdx(ip))/2.0 - 
     &             (gdx(i)+gdx(im))/2.0 )*Re*cos(theta(j))*pi/180.0
	dyr(i,j) =((gdy(jp)+gdy(j))/2.0 -
     &             (gdy(j)+gdy(jm))/2.0 )*Re*pi/180.0
	enddo
	enddo
	
	do i=2,imt-1
	do j=2,jmt-1
	ip = i + 1
	im = i - 1 
	jp = j + 1
	jm = j - 1 
	dxu(i,j) = (dxr(i,j)+dxr(i,jp))/2.0
	dyu(i,j) = (dyr(i,jp)+dyr(i,j))/2.0
	dxv(i,j) = (dxr(ip,j)+dxr(i,j))/2.0
	dyv(i,j) = (dyr(i,j)+dyr(ip,j))/2.0
	dxh(i,j) = (dxr(i,j)+dxr(i,jp))/2.0
	dyh(i,j) = (dyr(i,j)+dyr(ip,j))/2.0
	dah(i,j) = (dxh(i,j)*dyh(i,j))
	enddo
	enddo
	
	do i=2,imt-1
	 dxu(i,1) = dxu(i,2)
	 dyu(i,1) = dyu(i,2)
	 dxv(i,1) = dxv(i,2)
	 dyv(i,1) = dyv(i,2)
	 dxh(i,1) = dxh(i,2)
	 dyh(i,1) = dyh(i,2)
	 dah(i,1) = dah(i,2)

	 dxu(i,jmt) = dxu(i,jmt-1)
	 dyu(i,jmt) = dyu(i,jmt-1)
	 dxv(i,jmt) = dxv(i,jmt-1)
	 dyv(i,jmt) = dyv(i,jmt-1)
	 dxh(i,jmt) = dxh(i,jmt-1)
	 dyh(i,jmt) = dyh(i,jmt-1)
	 dah(i,jmt) = dah(i,jmt-1)
	enddo

	do j=2,jmt-1
	 dxu(1,j) = dxu(2,j)
	 dyu(1,j) = dyu(2,j)
	 dxv(1,j) = dxv(2,j)
	 dyv(1,j) = dyv(2,j)
	 dxh(1,j) = dxh(2,j)
	 dyh(1,j) = dyh(2,j)
	 dah(1,j) = dah(2,j)

	 dxu(imt,j) = dxu(imt-1,j)
	 dyu(imt,j) = dyu(imt-1,j)
	 dxv(imt,j) = dxv(imt-1,j)
	 dyv(imt,j) = dyv(imt-1,j)
	 dxh(imt,j) = dxh(imt-1,j)
	 dyh(imt,j) = dyh(imt-1,j)
	 dah(imt,j) = dah(imt-1,j)
	enddo
	

	 dxu(1,1) = dxu(2,2)
	 dyu(1,1) = dyu(2,2)
	 dxv(1,1) = dxv(2,2)
	 dyv(1,1) = dyv(2,2)
	 dxh(1,1) = dxh(2,2)
	 dyh(1,1) = dyh(2,2)
	 dah(1,1) = dah(2,2)

	 dxu(imt,jmt) = dxu(imt-1,jmt-1)
	 dyu(imt,jmt) = dyu(imt-1,jmt-1)
	 dxv(imt,jmt) = dxv(imt-1,jmt-1)
	 dyv(imt,jmt) = dyv(imt-1,jmt-1)
	 dxh(imt,jmt) = dxh(imt-1,jmt-1)
	 dyh(imt,jmt) = dyh(imt-1,jmt-1)
	 dah(imt,jmt) = dah(imt-1,jmt-1)

	 dxu(imt,1) = dxu(imt-1,1)
	 dyu(imt,1) = dyu(imt-1,1)
	 dxv(imt,1) = dxv(imt-1,1)
	 dyv(imt,1) = dyv(imt-1,1)
	 dxh(imt,1) = dxh(imt-1,1)
	 dyh(imt,1) = dyh(imt-1,1)
	 dah(imt,1) = dah(imt-1,1)

	 dxu(1,jmt) = dxu(1,jmt-1)
	 dyu(1,jmt) = dyu(1,jmt-1)
	 dxv(1,jmt) = dxv(1,jmt-1)
	 dyv(1,jmt) = dyv(1,jmt-1)
	 dxh(1,jmt) = dxh(1,jmt-1)
	 dyh(1,jmt) = dyh(1,jmt-1)
	 dah(1,jmt) = dah(1,jmt-1)

	return
	stop
	end
