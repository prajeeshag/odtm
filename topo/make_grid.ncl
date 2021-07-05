; USER INPUT STARTS

nz_x = 3	; nmber of zones in the topography in zonal direction
nz_y = 3	; nmber of zones in the topography in meridionaldirection

x = (/30,80,90,120/)	; start and end longitudes of each zone
dx = (/0.5,0.25,0.5/)	; resolution within each zone

y = (/-30,-15,8,30/)	; start and end latitudes of each zone
dy = (/0.5,0.25,0.5/)	; resolution within each zone

nxsm = (/5,5,5,5/)	; thickness (in degrees) of blending zone 
nysm = (/5,5,5,5/)	; thickness (in degrees) of blending zone 

; USER INPUT ENDS
;----------------------------------------------------------------------


 zdx = x(1:nz_x)-x(0:nz_x-1)
 zdy = y(1:nz_y)-y(0:nz_y-1)

 zdxs = new(nz_x+nz_x-1,float) 
 zdys = new(nz_y+nz_y-1,float) 

 xs = new(nz_x+nz_x,float) 
 ys = new(nz_y+nz_y,float) 

 k = -1
 do i = 0, nz_x+nz_x-2, 2
    k = k + 1

    zdxs(i) = zdx(k)
    ;right
    if (i.lt.nz_x+nz_x-2) then
       zdxs(i) = zdxs(i) - dx(k) * nxsm(i)
    end if

    ;left
    if (i.gt.0) then
       zdxs(i) = zdxs(i) - dx(k) * nxsm(i-1)
    end if
 end do

 k = -1
 do i = 1, nz_x+nz_x-3, 2
    k = k + 1
    zdxs(i) = dx(k) * nxsm(i-1)
    zdxs(i) = zdxs(i) + dx(k+1) * nxsm(i)
 end do

 xs(0) = x(0) 
 do i = 1, nz_x+nz_x-1
    xs(i) = xs(i-1) + zdxs(i-1)
 end do



 k = -1
 do i = 0, nz_y+nz_y-2, 2
    k = k + 1

    zdys(i) = zdy(k)
    ;right
    if (i.lt.nz_y+nz_y-2) then
       zdys(i) = zdys(i) - dy(k) * nysm(i)
    end if

    ;left
    if (i.gt.0) then
       zdys(i) = zdys(i) - dy(k) * nysm(i-1)
    end if
 end do

 k = -1
 do i = 1, nz_y+nz_y-3, 2
    k = k + 1
    zdys(i) = dy(k) * nysm(i-1)
    zdys(i) = zdys(i) + dy(k+1) * nysm(i)
 end do

 ys(0) = y(0) 
 do i = 1, nz_y+nz_y-1
    ys(i) = ys(i-1) + zdys(i-1)
 end do

 tmp = xs(0)
 zone=0
 k = 0
 nele=0
 do while (zone.lt.dimsizes(zdxs))
    if (zone%2.eq.0) then
       ddx = dx(k)
       tmp = tmp + ddx
       if (tmp.le.xs(zone+1)) then
          nele = nele + 1
       else
		  print(zone+"")
		  zone = zone + 1
		  k = k + 1
		  tmp = tmp - ddx
	       end if
	    else
	       inc = (dx(k) - dx(k-1))/(nxsm(k)+nxsm(k-1))
	       ddx = ddx+inc
	       if (inc.lt.0 .and. ddx.lt.dx(k)) then
		  ddx = dx(k)
	       end if
	       if (inc.gt.0 .and. ddx.gt.dx(k)) then
		  ddx = dx(k)
	       end if
	       tmp = tmp + ddx
	       if (tmp.le.xs(zone+1)) then
		  nele = nele + 1
	       else
		  print(zone+"")
          zone = zone + 1
          tmp = tmp - ddx
       end if 
    end if
 end do
 print(nele+"")

 xcord = new(nele,float)


 tmp = xs(0)
 zone=0
 k = 0
 nele=-1
 do while (zone.lt.dimsizes(zdxs))
    if (zone%2.eq.0) then
       ddx = dx(k)
       tmp = tmp + ddx
       if (tmp.le.xs(zone+1)) then
          nele = nele + 1
          xcord(nele) = tmp
       else
          print(zone+"")
          zone = zone + 1
          k = k + 1
          tmp = tmp - ddx
       end if
    else
       inc = (dx(k) - dx(k-1))/(nxsm(k)+nxsm(k-1))
       ddx = ddx+inc
       if (inc.lt.0 .and. ddx.lt.dx(k)) then
          ddx = dx(k)
       end if
       if (inc.gt.0 .and. ddx.gt.dx(k)) then
          ddx = dx(k)
       end if
       if(zone.eq.1) then
       end if
       tmp = tmp + ddx
       if(zone.eq.1) then
       end if
       if (tmp.le.xs(zone+1)) then
          nele = nele + 1
          xcord(nele) = tmp
       else
          print(zone+"")
          zone = zone + 1
          tmp = tmp - ddx
       end if 
    end if
 end do


 tmp = ys(0)
 zone=0
 k = 0
 nele=0
 do while (zone.lt.dimsizes(zdys))
    if (zone%2.eq.0) then
       ddy = dy(k)
       tmp = tmp + ddy
       if (tmp.le.ys(zone+1)) then
          nele = nele + 1
       else
          print(zone+"")
          zone = zone + 1
          k = k + 1
          tmp = tmp - ddy
       end if
    else
       inc = (dy(k) - dy(k-1))/(nysm(k)+nysm(k-1))
       ddy = ddy+inc
       if (inc.lt.0 .and. ddy.lt.dy(k)) then
          ddy = dy(k)
       end if
       if (inc.gt.0 .and. ddy.gt.dy(k)) then
          ddy = dy(k)
       end if
       tmp = tmp + ddy
       if (tmp.le.ys(zone+1)) then
          nele = nele + 1
       else
          print(zone+"")
          zone = zone + 1
          tmp = tmp - ddy
       end if 
    end if
 end do
 print(nele+"")

 ycord = new(nele,float)


 tmp = ys(0)
 zone=0
 k = 0
 nele=-1
 do while (zone.lt.dimsizes(zdys))
    if (zone%2.eq.0) then
       ddy = dy(k)
       tmp = tmp + ddy
       if (tmp.le.ys(zone+1)) then
          nele = nele + 1
          ycord(nele) = tmp
       else
          print(zone+"")
          zone = zone + 1
          k = k + 1
          tmp = tmp - ddy
       end if
    else
       inc = (dy(k) - dy(k-1))/(nysm(k)+nysm(k-1))
       ddy = ddy+inc
       if (inc.lt.0 .and. ddy.lt.dy(k)) then
          ddy = dy(k)
       end if
       if (inc.gt.0 .and. ddy.gt.dy(k)) then
          ddy = dy(k)
       end if
       if(zone.eq.1) then
       end if
       tmp = tmp + ddy
       if(zone.eq.1) then
       end if
       if (tmp.le.ys(zone+1)) then
          nele = nele + 1
          ycord(nele) = tmp
       else
          print(zone+"")
          zone = zone + 1
          tmp = tmp - ddy
       end if 
    end if
 end do


 xcord@units="degrees_east"
 ycord@units="degrees_north"
 nx = dimsizes(xcord)
 xcordt = (xcord(1:nx-1)+xcord(0:nx-2))*0.5
 ny = dimsizes(ycord)
 ycordt = (ycord(1:ny-1)+ycord(0:ny-2))*0.5
 xcordt@units="degrees_east"
 ycordt@units="degrees_north"


 f = addfile("etopo5.cdf","r")
 var = f->ROSE
 xi = f->ETOPO05_X
 yi = f->ETOPO05_Y

 topo = linint2_Wrap(xi, yi, var, True, xcordt, ycordt, 0)

 topo1 = where(topo.lt.-100,1.,0.)
 copy_VarMeta(topo,topo1)
 
 geolon_t = new((/nx-1,ny-1/),float)
 geolat_t = new((/nx-1,ny-1/),float)

 geolon_b = new((/nx,ny/),float)
 geolat_b = new((/nx,ny/),float)

 do i = 0, ny-2
    geolon_t(:,i) = xcordt
 end do
 do i = 0, nx-2
    geolat_t(i,:) = ycordt
 end do
 do i = 0, ny-1
    geolon_b(:,i) = xcord
 end do
 do i = 0, nx-1
    geolat_b(i,:) = ycord
 end do

 geolat_t!0="xt"
 geolat_t!1="yt"
 geolon_t!0="xt"
 geolon_t!1="yt"
 topo1!0 = "yt"
 topo1!1 = "xt"
 geolat_b!0="xb"
 geolat_b!1="yb"
 geolon_b!0="xb"
 geolon_b!1="yb"
 geolon_b!1="yb"

 printVarSummary(transpose(topo1))
 printVarSummary(geolon_t)
 
 system("rm -f grid_spec.nc")
 fo = addfile("grid_spec.nc","c")
 fo->geolon_t = transpose(geolon_t)
 fo->geolat_t = transpose(geolat_t)
 fo->geolon_b = transpose(geolon_b)
 fo->geolat_b = transpose(geolat_b)
 fo->rkmt = topo1

