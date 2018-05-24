#!/bin/bash

files1=$1
files2=$2

if [ -z "$file1" ]; then

	files=$(ls -rt output*.nc)

	file1=$(echo $files | awk '{print $1}')

	file2=$(echo $files | awk '{print $(NF-1)}')
fi


echo
echo
echo
echo ------------------------------------
echo $file1 $file2
echo ------------------------------------
echo
echo
echo

cat <<EOF > check_dev.ncl

begin

	file1 = "$file1"

	file2 = "$file2"

	f1 = addfile(file1,"r")
	f2 = addfile(file2,"r")

	temp1 = f1->temp
	temp_mld1 = f1->temp_mld

	salt1 = f1->salt
	salt_mld1 = f1->salt_mld

	h1 = f1->h

	u1 = f1->u

	v1 = f1->v

	eta1 = f1->eta

    	
	temp2 = f2->temp
	temp_mld2 = f2->temp_mld

	salt2 = f2->salt
	salt_mld2 = f2->salt_mld

	h2 = f2->h

	u2 = f2->u

	v2 = f2->v

	eta2 = f2->eta

	lon = f2->lon
	nlon = dimsizes(lon)
	lat = f2->lat
	nlat = dimsizes(lat)

	printVarSummary(temp1)

    wks   = gsn_open_wks ("pdf", "check_dev_plots")

  	resC                       = True                 ; plot mods desired
  	resC@gsnDraw = True
  	resC@gsnFrame = True
  	resC@mpFillOn = True
  	resC@cnFillOn = True
  	resC@cnLinesOn = False
  	resC@mpProjection  = "CylindricalEquidistant"
  	resC@mpMinLonF = lon(0)
  	resC@mpMaxLonF = lon(nlon-1)
  	resC@mpMinLatF = lat(0)
  	resC@mpMaxLatF = lat(nlat-1)
	resC@gsnMaximize = True
  	resC@gsnAddCyclic          = False                 ; add cyclic point
  	;resC@cnFillPalette        = "BlWhRe"
  	;symMinMaxPlt(eof,20,False,resC)

	;----------------------(Surface)---------------------------
	resC@gsnLeftString = "SST"
  	plot = gsn_csm_contour_map(wks,temp_mld2(0,0,:,:),resC)

	resC@gsnLeftString = "SSS"
  	plot = gsn_csm_contour_map(wks,salt_mld2(0,0,:,:),resC)

	resC@gsnLeftString = "H"
  	plot = gsn_csm_contour_map(wks,h2(0,0,:,:),resC)
    	
	resC@gsnLeftString = "eta"
  	plot = gsn_csm_contour_map(wks,eta2(0,:,:),resC)

	resC@gsnLeftString = "U"
  	plot = gsn_csm_contour_map(wks,u2(0,0,:,:),resC)

	resC@gsnLeftString = "V"
  	plot = gsn_csm_contour_map(wks,v2(0,0,:,:),resC)


    ;-------------Tendencies (Surface) -----------------------	

	diff = temp_mld2(0,0,:,:) - temp_mld1(0,0,:,:)
	copy_VarCoords(temp_mld2(0,0,:,:),diff)
  	resC@cnFillPalette        = "BlWhRe"
  	symMinMaxPlt(diff,20,False,resC)
	resC@gsnLeftString = "SST differece"
  	plot = gsn_csm_contour_map(wks,diff,resC)

	diff = salt_mld2(0,0,:,:) - salt_mld1(0,0,:,:)
	copy_VarCoords(temp_mld2(0,0,:,:),diff)
  	resC@cnFillPalette        = "BlWhRe"
  	symMinMaxPlt(diff,20,False,resC)
	resC@gsnLeftString = "SSS differece"
  	plot = gsn_csm_contour_map(wks,diff,resC)

	diff = h2(0,0,:,:) - h1(0,0,:,:)
	copy_VarCoords(temp_mld2(0,0,:,:),diff)
  	resC@cnFillPalette        = "BlWhRe"
  	symMinMaxPlt(diff,20,False,resC)
	resC@gsnLeftString = "H differece"
  	plot = gsn_csm_contour_map(wks,diff,resC)

	diff = eta2(0,:,:) - eta1(0,:,:)
	copy_VarCoords(temp_mld2(0,0,:,:),diff)
  	resC@cnFillPalette        = "BlWhRe"
  	symMinMaxPlt(diff,20,False,resC)
	resC@gsnLeftString = "eta differece"
  	plot = gsn_csm_contour_map(wks,diff,resC)

	diff = u2(0,0,:,:) - u1(0,0,:,:)
	copy_VarCoords(temp_mld2(0,0,:,:),diff)
  	resC@cnFillPalette        = "BlWhRe"
  	symMinMaxPlt(diff,20,False,resC)
	resC@gsnLeftString = "U differece"
  	plot = gsn_csm_contour_map(wks,diff,resC)

	diff = v2(0,0,:,:) - v1(0,0,:,:)
	copy_VarCoords(temp_mld2(0,0,:,:),diff)
  	resC@cnFillPalette        = "BlWhRe"
  	symMinMaxPlt(diff,20,False,resC)
	resC@gsnLeftString = "V differece"
  	plot = gsn_csm_contour_map(wks,diff,resC)

end
EOF

ncl check_dev.ncl 
