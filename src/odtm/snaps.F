	subroutine write_output (loop_index)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       subroutine to write snapshots of model output in netcdf format.
c       output can be restircted to selected variables in Makefile by
c       enabling curresponding options in the option list.
c 
c       author: vinu k. v.            email: vinu@ees.hokudai.ac.jp
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        use size_mod, only : dxu, dz, eta, gdx, gdy, h, h_atmos, h_output, i, j
        use size_mod, only : loop, month_start_snap, pvort, pvort_output, renergy
        use size_mod, only :  rkmu, rkmv, taun1, taux_force, tauy_force, tx_output
        use size_mod, only : ty_output, u, u_atmos, u_output, v, v_atmos, v_output
        use size_mod, only : we, we_output, diag_ext1, diag_ext1_output
        use size_mod, only : diag_ext2, diag_ext2_output, diag_ext3, diag_ext4
        use size_mod, only : diag_ext5, diag_ext6, dyv, k, rkmh, salt, salt_output 
        use size_mod, only : shcoeff, shcoeff_output, smcoeff, smcoeff_output, t
        use size_mod, only : taup, temp, temp_output, tr01_output, tr02_output, uvel
        use size_mod, only : uvel_output, vvel, vvel_output, dens_output
        use size_mod, only : diag_ext3_output, diag_ext4_output, diag_ext5_output
        use size_mod, only : diag_ext6_output, rmld_misc_output, taun, nn, imt, jmt, km
        use param_mod, only : h_atmos_snap, hh, rdepth, rdepth_mld
        use param_mod, only : salt_snap, start3, start4, startinv, startp, tauxx
        use param_mod, only : tauyy, tr, u_atmos_snap, undef, uu, v_atmos_snap, vv
        use param_mod, only : day2sec, dt, ncid, ncid1, ncid2, pvorticity, rlat
        use param_mod, only : rlon, rtime, vetaid, vmld4did1, vmld4did2, vmld4did3
        use param_mod, only : vmld4did4, vmldmldid, vmldqnetid, vmldsaltid, vmldshid
        use param_mod, only : vmldsmid, vmldtempid, vmlduvelid, vmldvvelid, vtempid
        use param_mod, only : vtid, vtimeid, vtimeid1, vtimeid2, vuid, vvid
        use param_mod, only : idepthidmld2, itimeid2, latid2, lonid2, npp
        use param_mod, only : vardim, vardim2, vardim3, vardim33, vardimmld, varinvdim
        use param_mod, only : varpdim, vdensid, vmld4did5, vmld4did6, vsaltid, vtauxid
        use param_mod, only : vtauyid, vvortid, count2, count3, count4, mask
        use param_mod, only : stamp_day, stamp_time, start2, temp_mld, temp_snap
        use param_mod, only : vardimmld3, vardimmldext, vdepthid, vdepthidmld
        use param_mod, only : vdepthidmld2, vdxuid, vdyvid, vlatid, vlatid1, vlatid2
        use param_mod, only : vlonid, vlonid1, vlonid2, vmaskid, count4mld
        use param_mod, only : diag_ext1_mld, diag_ext2_mld, diag_ext3_mld 
        use param_mod, only : diag_ext4_mld, diag_ext5_mld, diag_ext6_mld, mld_mld
        use param_mod, only : salt_mld, shcoeff_mld, smcoeff_mld, uvel_mld, vvel_mld


        implicit none
	include 'netcdf.inc'
    
	    real reta(imt,km,jmt,1)
        integer :: icount, idepthid, idepthidmld, ii, invid, ir, itimeid, itimeid1
        integer :: jr, kk, latid, latid1, lonid, lonid1, loop_index
        integer :: nt
        integer :: pid, rcode 

	taun1 = 2
	
	do nt=1,nn
	write(Tr(nt),'(a,i2.2)') 'Tr',nt
	enddo
	
c
c      storing starting dimension values
c
	start4(1) = 1
	start4(2) = 1
	start4(3) = 1
	start4(4) = loop_index
	start3(1) = 1
	start3(2) = 1
	start3(3) = loop_index
	startp(1) = 1
	startp(2) = loop_index
	startinv(1) = 1
	startinv(2) = loop_index

	rdepth(1) = dz(1)
	do ii=2,km
	rdepth(ii)=rdepth(ii-1) + dz(ii)
	enddo
        do ii=1,51
        rdepth_mld(ii) = ii*5.0
        enddo
	

c define undef value
	undef = 0.0

c NOTE*** the array ave(km) is used to find inter-facial depth as
c anomaly from the average depth of the respective layers. i.e. the
c at each level the suface elevation is written as h-h(average.of.domain).
	do kk=1,km
	enddo
c NOTE***
	
	i=0
	do ir=1,imt
	 i=i+1
	 j=0
	do jr=1,jmt
	 j=j+1
	do kk = 1,km
	 tauxx(i,j)=taux_force(ir,jr) !taux_snap(ir,jr)
	 tauyy(i,j)=tauy_force(ir,jr) !
	 u_atmos_snap(i,j) = u_atmos(i,j,taup)
	 v_atmos_snap(i,j) = v_atmos(i,j,taup)
	 h_atmos_snap(i,j) = h_atmos(i,j,taup)
	 uu(i,kk,j) = u(ir,kk,jr,2) !*rkmu(ir,jr)*rkmv(ir,jr)
	 vv(i,kk,j) = v(ir,kk,jr,2) !*rkmu(ir,jr)*rkmv(ir,jr)
	 hh(i,kk,j) = h(ir,kk,jr,2) !*rkmu(ir,jr)*rkmv(ir,jr)
	 salt_snap(i,kk,j) = we(ir,kk,jr)
	
#ifdef output_average
        tauxx(i,j)=tx_output(i,j) 
        tauyy(i,j)=ty_output(i,j) 
        uu(i,kk,j)=u_output(i,kk,j)
        vv(i,kk,j)=v_output(i,kk,j)
        hh(i,kk,j)=h_output(i,kk,j)
	salt_snap(i,kk,j) = we_output(ir,kk,jr)
#endif
	
	 if (eta(i,kk,j,1) .ne. 0.0) then
	 reta(i,kk,j,1) = eta(i,kk,j,1) * rkmu(i,j)*rkmv(i,j)
c	 reta(i,kk,j,1) = (eta(i,kk,j,1) - 0.8052951694)*rrkmt(i,j) ! CASE-2
	 else
	 reta(i,kk,j,1) = 0.0
	 endif
	
	 rEnergy(ir,kk,jr) = 0.0
	 pvorticity (i,kk,j) = pvort(i,kk,j) !*rkmh(i,j)
#ifdef output_average
         pvorticity (i,kk,j) = pvort_output(i,kk,j)
#endif
	
	enddo
	enddo
	enddo

c assign lon, lat, depth and time arrays
	icount=0
	do ii=1,imt
	icount=icount+1
	rlon(icount)= gdx(icount)
	enddo
	icount=0
	do ii=1,jmt
	icount=icount+1
	rlat(icount)=gdy(icount)
	enddo
	rtime(1)=loop*1.0*dt/day2sec



c open netcdf file
	if (loop_index .eq. 1) ncid = nccre ('snapshots.nc',ncclob,rcode)
	if (loop_index .eq. 1) ncid1 = nccre ('snaps_mld.nc',ncclob,rcode)
#ifdef snap_mld_extended
	if (loop_index .eq. 1) ncid2 = nccre ('snaps_mld_ext.nc',ncclob,rcode)
#endif
	if (loop_index .ne. 1) then 
	ncid = ncopn ('snapshots.nc',ncwrite,rcode)
	ncid1 = ncopn ('snaps_mld.nc',ncwrite,rcode)
#ifdef snap_mld_extended
	ncid2 = ncopn ('snaps_mld_ext.nc',ncwrite,rcode)
#endif
	call ncredf (ncid, rcode)
	call ncredf (ncid1, rcode)
#ifdef snap_mld_extended
	call ncredf (ncid2, rcode)
#endif
	vtimeid = ncvid (ncid,'time',rcode)
	vtimeid1 = ncvid (ncid1,'time',rcode)
	vtimeid2 = ncvid (ncid2,'time',rcode)
#ifdef snap_ssh
	vtid = ncvid (ncid, 'h', rcode)
	vetaid = ncvid (ncid, 'eta', rcode)
#endif
#ifdef snap_velocity
	vuid = ncvid (ncid, 'u', rcode)
	vvid = ncvid (ncid, 'v', rcode)
#endif
#ifdef snap_tracer
	do nt=1,nn
	vtempid(nt) = ncvid (ncid, Tr(nt), rcode)
	enddo
#endif
#ifdef snap_mld
	vmldtempid = ncvid (ncid1, 'temp', rcode)
	vmldsaltid = ncvid (ncid1, 'salt', rcode)
	vmlduvelid = ncvid (ncid1, 'uvel', rcode)
	vmldvvelid = ncvid (ncid1, 'vvel', rcode)
	vmldmldid = ncvid (ncid1, 'mld', rcode)
        vmldshid = ncvid (ncid1, 'SH', rcode)
        vmldsmid = ncvid (ncid1, 'SM', rcode)
        vmldqnetid = ncvid (ncid1, 'diag', rcode)
#endif
#ifdef snap_mld_extended
        vmld4did1 = ncvid (ncid2, 'tke', rcode)
        vmld4did2 = ncvid (ncid2, 'Rif', rcode)
        vmld4did3 = ncvid (ncid2, 'mlen', rcode)
        vmld4did4 = ncvid (ncid2, 'St_h', rcode)
        vmld4did5 = ncvid (ncid2, 'St_m', rcode)
        vmld4did6 = ncvid (ncid2, 'diag_ext', rcode)
#endif
#ifdef snap_we
	vsaltid = ncvid (ncid, 'we', rcode)
#endif
#ifdef snap_density
	vdensid = ncvid (ncid, 'dens', rcode)
#endif
#ifdef snap_vorticity
	vvortid = ncvid (ncid, 'pvort', rcode)
#endif
#ifdef snap_tau
	vtauxid = ncvid (ncid, 'tx' ,rcode)
	vtauyid = ncvid (ncid, 'ty' ,rcode)
#endif
#ifdef atmosphere
	vatmos_h_id = ncvid (ncid, 'h_atmos',rcode)
	vatmos_u_id = ncvid (ncid, 'u_atmos',rcode)
	vatmos_v_id = ncvid (ncid, 'v_atmos',rcode)
#endif
#ifdef snap_particle
	vpxid = ncvid (ncid, 'xp', rcode)
	vpyid = ncvid (ncid, 'yp', rcode)
#endif
#ifdef snap_inverse_model
	vthetaid = ncvid (ncid, 'Th',rcode)
#endif
	go to 101
	endif

c define dimension
	lonid = ncddef (ncid,'longitude',imt,rcode)
	latid = ncddef (ncid,'latitude',jmt,rcode)
	idepthid = ncddef (ncid,'depth',km,rcode)
	itimeid = ncddef (ncid,'time',ncunlim,rcode)
	pid = ncddef (ncid,'array',NPP,rcode)
	invid = ncddef (ncid,'array1',nn,rcode)
	lonid1 = ncddef (ncid1,'longitude',imt,rcode)
	latid1 = ncddef (ncid1,'latitude',jmt,rcode)
        idepthidmld = ncddef (ncid1,'depth',51,rcode)
	itimeid1 = ncddef (ncid1,'time',ncunlim,rcode)
#ifdef snap_mld_extended
	lonid2 = ncddef (ncid2,'longitude',imt,rcode)
	latid2 = ncddef (ncid2,'latitude',jmt,rcode)
        idepthidmld2 = ncddef (ncid2,'depth',51,rcode)
	itimeid2 = ncddef (ncid2,'time',ncunlim,rcode)
#endif
c store dimension id into arrays
	vardim(1) = lonid
	vardim(2) = idepthid
	vardim(3) = latid
	vardim(4) = itimeid
	
	varpdim (1) = pid
	varpdim (2) = itimeid

	varinvdim (1) = invid
	varinvdim (2) = itimeid

	vardim3(1) = lonid
	vardim3(2) = latid
	vardim3(3) = itimeid

	vardim2(1) = lonid
	vardim2(2) = latid

	vardim33(1) = lonid
	vardim33(2) = idepthid
	vardim33(3) = latid

	vardimmld(1) = lonid1
	vardimmld(2) = idepthidmld
	vardimmld(3) = latid1
	vardimmld(4) = itimeid1
	vardimmld3(1) = lonid1
	vardimmld3(2) = latid1
	vardimmld3(3) = itimeid1
#ifdef  snap_mld_extended
        vardimmldext(1) = lonid2
        vardimmldext(2) = idepthidmld2
        vardimmldext(3) = latid2
        vardimmldext(4) = itimeid2
#endif
c define variable id
	vlonid = ncvdef (ncid,'longitude',ncfloat,1,lonid,rcode)
	vlatid = ncvdef (ncid,'latitude',ncfloat,1,latid,rcode)
	vdepthid = ncvdef (ncid,'depth',ncfloat,1,idepthid,rcode)
	vtimeid = ncvdef (ncid,'time',ncfloat,1,itimeid,rcode)
	vlonid1 = ncvdef (ncid1,'longitude',ncfloat,1,lonid1,rcode)
	vlatid1 = ncvdef (ncid1,'latitude',ncfloat,1,latid1,rcode)
	vdepthidmld = ncvdef (ncid1,'depth',ncfloat,1,idepthidmld,rcode)
	vtimeid1 = ncvdef (ncid1,'time',ncfloat,1,itimeid1,rcode)
#ifdef snap_mld_extended
	vlonid2 = ncvdef (ncid2,'longitude',ncfloat,1,lonid2,rcode)
	vlatid2 = ncvdef (ncid2,'latitude',ncfloat,1,latid2,rcode)
	vdepthidmld2 = ncvdef (ncid2,'depth',ncfloat,1,idepthidmld2,rcode)
	vtimeid2 = ncvdef (ncid2,'time',ncfloat,1,itimeid2,rcode)
#endif
c define attribute
        call ncaptc(ncid,vlatid,'units',ncchar,13,
     $    'degrees_north',rcode)
	call ncaptc(ncid,vlonid,'units',ncchar,12,
     $    'degrees_east',rcode)
	 call ncaptc(ncid,vdepthid,'units',ncchar,6,
     $    'meters',rcode)
        call ncaptc(ncid1,vlatid1,'units',ncchar,13,
     $    'degrees_north',rcode)
	call ncaptc(ncid1,vlonid1,'units',ncchar,12,
     $    'degrees_east',rcode)
	 call ncaptc(ncid1,vdepthidmld,'units',ncchar,6,
     $    'meters',rcode)
#ifdef snap_mld_extended
        call ncaptc(ncid2,vlatid2,'units',ncchar,13,
     $    'degrees_north',rcode)
	call ncaptc(ncid2,vlonid2,'units',ncchar,12,
     $    'degrees_east',rcode)
	 call ncaptc(ncid2,vdepthidmld2,'units',ncchar,6,
     $    'meters',rcode)
#endif
#ifdef snap_ssh
	vtid = ncvdef (ncid,'h',ncfloat,4,vardim,rcode)
	vetaid = ncvdef (ncid,'eta',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vtid,'units',ncchar,1,'M',rcode)
	call ncaptc(ncid,vetaid,'units',ncchar,1,'M',rcode)
	call ncapt(ncid,vtid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vetaid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_velocity
	vuid = ncvdef (ncid,'u',ncfloat,4,vardim,rcode)
	vvid = ncvdef (ncid,'v',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vuid,'units',ncchar,3,'M/S',rcode)
	call ncaptc(ncid,vvid,'units',ncchar,3,'M/S',rcode)
	call ncapt(ncid,vuid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vvid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_particle
	vpxid = ncvdef (ncid,'xp',ncfloat,2,varpdim,rcode)
	vpyid = ncvdef (ncid,'yp',ncfloat,2,varpdim,rcode)
	call ncaptc (ncid,vpxid,'units',ncchar,6,'number',rcode)
	call ncaptc (ncid,vpyid,'units',ncchar,6,'number',rcode)
#endif
#ifdef snap_inverse_model
	vthetaid = ncvdef (ncid,'Th',ncdouble,2,varinvdim,rcode)
	call ncaptc (ncid,vthetaid,'units',ncchar,4,'flux',rcode)
#endif
#ifdef snap_tau
	vtauxid = ncvdef (ncid,'tx',ncfloat,3,vardim3,rcode)
	vtauyid = ncvdef (ncid,'ty',ncfloat,3,vardim3,rcode)
	call ncaptc(ncid,vtauxid,'units',ncchar,9,'Newt/M2',rcode)
	call ncaptc(ncid,vtauyid,'units',ncchar,9,'Newt/M2',rcode)
	call ncapt(ncid,vtauxid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vtauyid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef atmosphere
	vatmos_h_id = ncvdef (ncid,'h_atmos',ncfloat,3,vardim3,rcode)
	vatmos_u_id = ncvdef (ncid,'u_atmos',ncfloat,3,vardim3,rcode)
	vatmos_v_id = ncvdef (ncid,'v_atmos',ncfloat,3,vardim3,rcode)
	call ncaptc(ncid,vatmos_h_id,'units',ncchar,5,'meter',rcode)
	call ncaptc(ncid,vatmos_u_id,'units',ncchar,3,'m/s',rcode)
	call ncaptc(ncid,vatmos_v_id,'units',ncchar,3,'m/s',rcode)
	call ncapt(ncid,vatmos_h_id,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vatmos_u_id,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid,vatmos_v_id,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_tracer
	do nt=1,nn
	vtempid(nt) = ncvdef (ncid,Tr(nt),ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vtempid(nt),'units',ncchar,5,'color',
     $  rcode)
	call ncapt(ncid,vtempid(nt),'missing_value',ncfloat
     $    ,1,undef,rcode)
	enddo
#endif
#ifdef snap_mld
	vmldtempid = ncvdef (ncid1,'temp',ncfloat,4,vardimmld,rcode)
	vmldsaltid = ncvdef (ncid1,'salt',ncfloat,4,vardimmld,rcode)
	vmlduvelid = ncvdef (ncid1,'uvel',ncfloat,4,vardimmld,rcode)
	vmldvvelid = ncvdef (ncid1,'vvel',ncfloat,4,vardimmld,rcode)
	vmldmldid = ncvdef (ncid1,'mld',ncfloat,3,vardimmld3,rcode)     
        vmldshid = ncvdef (ncid1,'SH', ncdouble,4,vardimmld,rcode)
        vmldsmid = ncvdef (ncid1,'SM', ncdouble,4,vardimmld,rcode)
        vmldqnetid = ncvdef (ncid1,'diag', ncfloat,4,vardimmld,rcode)
#ifdef snap_mld_extended
        write(*,*) ncid2, ncid1, ncid
        vmld4did1 = ncvdef (ncid2, 'tke', ncfloat,4,vardimmldext, rcode)
        vmld4did2 = ncvdef (ncid2, 'Rif', ncfloat,4,vardimmldext, rcode)
        vmld4did3 = ncvdef (ncid2, 'mlen', ncfloat,4,vardimmldext, rcode)
        vmld4did4 = ncvdef (ncid2, 'St_h', ncfloat,4,vardimmldext, rcode)
        vmld4did5 = ncvdef (ncid2, 'St_m', ncfloat,4,vardimmldext, rcode)
        vmld4did6 = ncvdef (ncid2, 'diag_ext', ncfloat,4,vardimmldext, rcode)
#endif
	call ncaptc(ncid1,vmldtempid,'units',ncchar,5,'deg-C',
     $  rcode)
	call ncaptc(ncid1,vmldsaltid,'units',ncchar,3,'psu',
     $  rcode)
	call ncaptc(ncid1,vmlduvelid,'units',ncchar,3,'m/s',
     $  rcode)
	call ncaptc(ncid1,vmldvvelid,'units',ncchar,3,'m/s',
     $  rcode)
	call ncaptc(ncid1,vmldmldid,'units',ncchar,1,'m',
     $  rcode)
	call ncaptc(ncid1,vmldshid,'units',ncchar,4,'m2/s',
     $  rcode)
	call ncaptc(ncid1,vmldsmid,'units',ncchar,4,'m2/s',
     $  rcode)
	call ncaptc(ncid1,vmldqnetid,'units',ncchar,4,'misc',
     $  rcode)
#ifdef snap_mld_extended
        call ncaptc(ncid2,vmld4did1,'units', ncchar,4,'tke ',rcode)
        call ncaptc(ncid2,vmld4did2,'units', ncchar,4,'Rif ',rcode)
        call ncaptc(ncid2,vmld4did3,'units', ncchar,4,'mlen',rcode)
        call ncaptc(ncid2,vmld4did4,'units', ncchar,4,'St_h',rcode)
        call ncaptc(ncid2,vmld4did5,'units', ncchar,4,'St_m',rcode)
        call ncaptc(ncid2,vmld4did6,'units', ncchar,4,'misc',rcode)
#endif
	call ncapt(ncid1,vmldtempid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmldsaltid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmlduvelid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmldvvelid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmldmldid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmldshid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmldsmid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	call ncapt(ncid1,vmldqnetid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#ifdef snap_mld_extended
        call ncapt(ncid2,vmld4did1, 'missing_value', ncfloat
     $    ,1,undef,rcode)
        call ncapt(ncid2,vmld4did2, 'missing_value', ncfloat
     $    ,1,undef,rcode)
        call ncapt(ncid2,vmld4did3, 'missing_value', ncfloat
     $    ,1,undef,rcode)
        call ncapt(ncid2,vmld4did4, 'missing_value', ncfloat
     $    ,1,undef,rcode)
        call ncapt(ncid2,vmld4did5, 'missing_value', ncfloat
     $    ,1,undef,rcode)
        call ncapt(ncid2,vmld4did6, 'missing_value', ncfloat
     $    ,1,undef,rcode)
#endif
#endif
#ifdef snap_we
	vsaltid = ncvdef (ncid,'we',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vsaltid,'units',ncchar,3,'m/s',rcode)
	call ncapt(ncid,vsaltid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_density
	vdensid = ncvdef (ncid,'dens',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vdensid,'units',ncchar,3,'psu',rcode)
	call ncapt(ncid,vdensid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
#ifdef snap_vorticity
	vvortid = ncvdef (ncid,'pvort',ncfloat,4,vardim,rcode)
	call ncaptc(ncid,vvortid,'units',ncchar,4,'/sec',rcode)
	call ncapt(ncid,vvortid,'missing_value',ncfloat
     $    ,1,undef,rcode)
#endif
	vmaskid = ncvdef (ncid,'mask',ncfloat,2,vardim2,rcode)
	call ncaptc (ncid,vmaskid,'units',ncchar,4,'mask',rcode)
	call ncapt(ncid,vmaskid,'missing_value',ncfloat
     $    ,1,undef,rcode)

	vdxuid = ncvdef (ncid,'dxu',ncfloat,2,vardim2,rcode)
	call ncaptc (ncid,vdxuid,'units',ncchar,4,'dxu',rcode)
	call ncapt(ncid,vdxuid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	
	vdyvid = ncvdef (ncid,'dyv',ncfloat,2,vardim2,rcode)
	call ncaptc (ncid,vdyvid,'units',ncchar,4,'dyv',rcode)
	call ncapt(ncid,vdyvid,'missing_value',ncfloat
     $    ,1,undef,rcode)
	
	
#ifdef monthly_climatology
        call ncaptc(ncid,vtimeid,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
        call ncaptc(ncid1,vtimeid1,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid1,vtimeid1,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
#ifdef snap_mld_extended
        call ncaptc(ncid2,vtimeid2,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid2,vtimeid2,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
#endif
#endif
#ifdef monthly_wind
        call ncaptc(ncid,vtimeid,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
        call ncaptc(ncid1,vtimeid1,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid1,vtimeid1,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
#ifdef snap_mld_extended
        call ncaptc(ncid2,vtimeid2,'units',ncchar,21,
     $     stamp_day(month_start_snap),rcode)
        call ncaptc(ncid2,vtimeid2,'time_origine',ncchar,20,
     $     stamp_time(month_start_snap),rcode)
#endif
#endif
#ifdef daily_wind
        call ncaptc(ncid,vtimeid,'units',ncchar,21,
     $     stamp_day_daily(iday_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'time_origine',ncchar,20,
     $     stamp_time_daily(iday_start_snap),rcode)
        call ncaptc(ncid,vtimeid,'calendar',ncchar,6,
     $     'noleap',rcode)
        call ncaptc(ncid1,vtimeid1,'units',ncchar,21,
     $     stamp_day_daily(iday_start_snap),rcode)
        call ncaptc(ncid1,vtimeid1,'time_origine',ncchar,20,
     $     stamp_time_daily(iday_start_snap),rcode)
        call ncaptc(ncid1,vtimeid1,'calendar',ncchar,6,
     $     'noleap',rcode)
#ifdef snap_mld_extended
        call ncaptc(ncid2,vtimeid2,'units',ncchar,21,
     $     stamp_day_daily(iday_start_snap),rcode)
        call ncaptc(ncid2,vtimeid2,'time_origine',ncchar,20,
     $     stamp_time_daily(iday_start_snap),rcode)
        call ncaptc(ncid2,vtimeid2,'calendar',ncchar,6,
     $     'noleap',rcode)
#endif
#endif
	
c
c Global attributes.
c
c        call ncaptc(ncid,nc_global,'title',ncchar,12,
c     $    'Model output',rcode)
c        call ncaptc(ncid,nc_global,'label',ncchar,33,
c     $    '1&1/2-layer_reduced_gravity model',rcode)
c        call ncaptc(ncid,nc_global,'history',ncchar,28,
c     $    'generated by module: snaps.F',rcode)
	
c leaving in define mode
	call ncendf (ncid, rcode)
	call ncendf (ncid1, rcode)
#ifdef snap_mld_extended
	call ncendf (ncid2, rcode)
#endif

c puting values
	call ncvpt (ncid,vlonid,1,imt,rlon,rcode)
	call ncvpt (ncid,vlatid,1,jmt,rlat,rcode)
	call ncvpt (ncid,vdepthid,1,km,rdepth,rcode)
	call ncvpt (ncid,vmaskid,start2,count2,mask,rcode)
	call ncvpt (ncid,vdxuid,start2,count2,dxu,rcode)
	call ncvpt (ncid,vdyvid,start2,count2,dyv,rcode)

	call ncvpt (ncid1,vlonid1,1,imt,rlon,rcode)
	call ncvpt (ncid1,vlatid1,1,jmt,rlat,rcode)
	call ncvpt (ncid1,vdepthidmld,1,51,rdepth_mld,rcode)
#ifdef snap_mld_extended
	call ncvpt (ncid2,vlonid2,1,imt,rlon,rcode)
	call ncvpt (ncid2,vlatid2,1,jmt,rlat,rcode)
	call ncvpt (ncid2,vdepthidmld2,1,51,rdepth_mld,rcode)
#endif
 101    continue
	 if ( loop_index .gt. 1) call ncendf (ncid, rcode)
	 if ( loop_index .gt. 1) call ncendf (ncid1, rcode)
#ifdef snap_mld_extended
	 if ( loop_index .gt. 1) call ncendf (ncid2, rcode)
#endif
#ifdef snap_tau
	call ncvpt (ncid,vtauxid,start3,count3,tauxx,rcode)
	call ncvpt (ncid,vtauyid,start3,count3,tauyy,rcode)
#endif
#ifdef atmosphere
	call ncvpt (ncid,vatmos_h_id,start3,count3,h_atmos_snap,rcode)
	call ncvpt (ncid,vatmos_u_id,start3,count3,u_atmos_snap,rcode)
	call ncvpt (ncid,vatmos_v_id,start3,count3,v_atmos_snap,rcode)
#endif
	call ncvpt (ncid,vtimeid,loop_index,1,rtime(1),rcode)
	call ncvpt (ncid1,vtimeid1,loop_index,1,rtime(1),rcode)
#ifdef snap_mld_extended
	call ncvpt (ncid2,vtimeid2,loop_index,1,rtime(1),rcode)
#endif
#ifdef snap_ssh
	call ncvpt (ncid,vtid,start4,count4,hh,rcode)
	call ncvpt (ncid,vetaid,start4,count4,reta,rcode)
#endif
#ifdef snap_velocity
	call ncvpt (ncid,vuid,start4,count4,uu,rcode)
	call ncvpt (ncid,vvid,start4,count4,vv,rcode)
#endif
#ifdef snap_tracer
        do nt=1,nn
        do i=1,imt
        do k=1,km
        do j=1,jmt
         temp_snap(i,k,j) = t(i,k,j,nt,taun)!*rrkmt(i,j)
        enddo
        enddo
        enddo
#ifdef output_average
        do i=1,imt
        do k=1,km
        do j=1,jmt
         if (nt .eq. 1) temp_snap(i,k,j) = tr01_output(i,k,j)!*rrkmt(i,j)
         if (nt .eq. 2) temp_snap(i,k,j) = tr02_output(i,k,j)!*rrkmt(i,j)
        enddo
        enddo
        enddo
#endif
        call ncvpt (ncid,vtempid(nt),start4,count4,temp_snap,rcode)
        enddo
#endif
#ifdef snap_mld
        do i=2,imt-1
        do j=2,jmt-1
        do k=1,51
         temp_mld(i,k,j) = temp(i,k,j,1) * rkmh(i,j)
         salt_mld(i,k,j) = salt(i,k,j,1) * rkmh(i,j)
         uvel_mld(i,k,j) = uvel(i,k,j,1) * rkmh(i,j)
         vvel_mld(i,k,j) = vvel(i,k,j,1) * rkmh(i,j)
         SMCoeff_mld(i,k,j) = SMCoeff(i,k,j) * rkmh(i,j)
         SHCoeff_mld(i,k,j) = SHCoeff(i,k,j) * rkmh(i,j)
#ifdef snap_mld_extended
         diag_ext1_mld(i,k,j) = diag_ext1(i,k,j) * rkmh(i,j)
         diag_ext2_mld(i,k,j) = diag_ext2(i,k,j) * rkmh(i,j)
         diag_ext3_mld(i,k,j) = diag_ext3(i,k,j) * rkmh(i,j)
         diag_ext4_mld(i,k,j) = diag_ext4(i,k,j) * rkmh(i,j)
         diag_ext5_mld(i,k,j) = diag_ext5(i,k,j) * rkmh(i,j)
         diag_ext6_mld(i,k,j) = diag_ext6(i,k,j) * rkmh(i,j)
#endif
#ifdef output_average
        temp_mld(i,k,j) = temp_output(i,k,j) * rkmh(i,j)
        salt_mld(i,k,j) = salt_output(i,k,j) * rkmh(i,j)
        uvel_mld(i,k,j) = uvel_output(i,k,j) * rkmh(i,j)
        vvel_mld(i,k,j) = vvel_output(i,k,j) * rkmh(i,j)
        SMCoeff_mld(i,k,j) = SMCoeff_output(i,k,j) * rkmh(i,j)
        SHCoeff_mld(i,k,j) = SHCoeff_output(i,k,j) * rkmh(i,j)
#ifdef snap_mld_extended
        diag_ext1_mld(i,k,j) = diag_ext1_output(i,k,j) * rkmh(i,j)
        diag_ext2_mld(i,k,j) = diag_ext2_output(i,k,j) * rkmh(i,j)
        diag_ext3_mld(i,k,j) = diag_ext3_output(i,k,j) * rkmh(i,j)
        diag_ext4_mld(i,k,j) = diag_ext4_output(i,k,j) * rkmh(i,j)
        diag_ext5_mld(i,k,j) = diag_ext5_output(i,k,j) * rkmh(i,j)
        diag_ext6_mld(i,k,j) = diag_ext6_output(i,k,j) * rkmh(i,j)
#endif
#endif
        enddo
         mld_mld(i,j) = mld_mld(i,j) * rkmh(i,j)
        enddo
        enddo
	call ncvpt (ncid1,vmldtempid,start4,count4mld,temp_mld,rcode)
	call ncvpt (ncid1,vmldsaltid,start4,count4mld,salt_mld,rcode)
	call ncvpt (ncid1,vmlduvelid,start4,count4mld,uvel_mld,rcode)
	call ncvpt (ncid1,vmldvvelid,start4,count4mld,vvel_mld,rcode)
	call ncvpt (ncid1,vmldmldid,start3,count3,mld_mld,rcode)
        call ncvpt (ncid1,vmldshid,start4,count4mld,SHCoeff_mld,rcode)
        call ncvpt (ncid1,vmldsmid,start4,count4mld,SMCoeff_mld,rcode)
#ifdef snap_mld_extended
        call ncvpt (ncid2,vmld4did1,start4,count4mld,diag_ext1_mld,rcode)
        call ncvpt (ncid2,vmld4did2,start4,count4mld,diag_ext2_mld,rcode)
        call ncvpt (ncid2,vmld4did3,start4,count4mld,diag_ext3_mld,rcode)
        call ncvpt (ncid2,vmld4did4,start4,count4mld,diag_ext4_mld,rcode)
        call ncvpt (ncid2,vmld4did5,start4,count4mld,diag_ext5_mld,rcode)
        call ncvpt (ncid2,vmld4did6,start4,count4mld,diag_ext6_mld,rcode)
#endif
#if defined output_average
        call ncvpt (ncid1,vmldqnetid,start4,count4mld,rmld_misc_output,rcode)
#else
        call ncvpt (ncid1,vmldqnetid,start4,count4mld,rmld_misc,rcode)
#endif
#endif
#ifdef snap_we
	call ncvpt (ncid,vsaltid,start4,count4,salt_snap,rcode)
#endif
#ifdef snap_density
#if defined output_average
	call ncvpt (ncid,vdensid,start4,count4,dens_output,rcode)
#else
	call ncvpt (ncid,vdensid,start4,count4,denss,rcode)
#endif
#endif
#ifdef snap_vorticity
	call ncvpt (ncid,vvortid,start4,count4,pvorticity,rcode)
#endif
#ifdef snap_particle
	call ncvpt (ncid,vpxid,startp,countp,xp,rcode)
	call ncvpt (ncid,vpyid,startp,countp,yp,rcode)
#endif
#ifdef snap_inverse_model
	call ncvpt (ncid,vthetaid,startinv,countinv,Thetainv,rcode)
#endif

c close netcdf file
	call ncclos (ncid, rcode)
	call ncclos (ncid1, rcode)
#ifdef snap_mld_extended
	call ncclos (ncid2, rcode)
#endif
	return
	end
