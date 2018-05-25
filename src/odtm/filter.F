	subroutine filter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	an implimentation of Asselin-Robert Filter
c	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod, only : eta, h, i, j, rkmh, rkmu, rkmv, t, u, v, k
        use size_mod, only : loop, taum, taun, taup, imt, jmt, km, nn
        use param_mod, only : alpha, dt
        use smooth_hanning_mod, only : smooth
        implicit none
	    real vel_lim(km), t_lim(k), t_llim(km)
        integer :: nt, idim5, idim4, idim3, idim2, idim1, ii, jj, kk, iii
        integer :: imte, imts, jjj, jmte, jmts, linear_switch
        real :: t_switch
	
        do k=1,km
	vel_lim(k) = 1.0
        enddo
        t_lim(1) = 27.0
        t_lim(2) = 20.0
        t_lim(3) = 14.0
        t_lim(4) = 14.0
        t_llim(1) = 07.0
        t_llim(2) = 10.0
        t_llim(3) = 05.0
        t_llim(4) = 01.0

	do i=1,imt
	do k=1,km-1
	do j=1,jmt

	 u(i,j,k,taun) = u(i,j,k,taun)+ alpha*rkmu(i,j)*(
     &    u(i,j,k,taup) -2*u(i,j,k,taun) + u(i,j,k,taum) )
	
         v(i,j,k,taun) = v(i,j,k,taun)+ alpha*rkmv(i,j)*(
     &    v(i,j,k,taup) -2*v(i,j,k,taun) + v(i,j,k,taum) )
	
	 do nt=1,nn
         t(i,j,k,nt,taun) = t(i,j,k,nt,taun)+ alpha*(
     &    t(i,j,k,nt,taup) -2*t(i,j,k,nt,taun) + t(i,j,k,nt,taum) )
	 enddo
	
	
        h(i,j,k,taun) = h(i,j,k,taun)+ alpha*rkmh(i,j)*(
     &    h(i,j,k,taup) -2*h(i,j,k,taun) + h(i,j,k,taum) )

#ifdef atmosphere
        u_atmos(i,j,taun) = u_atmos(i,j,taun)+ alpha*(
     &    u_atmos(i,j,taup) -2*u_atmos(i,j,taun) + u_atmos(i,j,taum) )
        v_atmos(i,j,taun) = v_atmos(i,j,taun)+ alpha*(
     &    v_atmos(i,j,taup) -2*v_atmos(i,j,taun) + v_atmos(i,j,taum) )
        h_atmos(i,j,taun) = h_atmos(i,j,taun)+ alpha*(
     &    h_atmos(i,j,taup) -2*h_atmos(i,j,taun) + h_atmos(i,j,taum) )
#endif
	
	u(i,j,k,taum) = u(i,j,k,taun)
	u(i,j,k,taun) = u(i,j,k,taup)
	u(i,j,k,taup) = 0.0

	
	v(i,j,k,taum) = v(i,j,k,taun)
	v(i,j,k,taun) = v(i,j,k,taup)
	v(i,j,k,taup) = 0.0

	h(i,j,k,taum) = h(i,j,k,taun)
	h(i,j,k,taun) = h(i,j,k,taup)
	h(i,j,k,taup) = 0.0

	eta(i,j,k,taum) = eta(i,j,k,taun)
	eta(i,j,k,taun) = eta(i,j,k,taup)
	eta(i,j,k,taup) = 0.0
	
	do nt=1,nn
	t(i,j,k,nt,taum) = t(i,j,k,nt,taun)
	t(i,j,k,nt,taun) = t(i,j,k,nt,taup)
	t(i,j,k,nt,taup) = 0.0

        if (k.gt.1) then !Prajeesh  
        if (t(i,j,k,1,taum).gt.t(i,j,k-1,1,taum))then
        t(i,j,k,1,taum) = (t(i,j,k-1,1,taum)+ t(i,j,k,1,taum))/2.0
        t(i,j,k-1,1,taum) = t(i,j,k,1,taum)
        endif
        endif

        if (k.gt.1) then !Prajeesh  
        if (t(i,j,k,1,taun).gt.t(i,j,k-1,1,taun))then
        t(i,j,k,1,taun) = (t(i,j,k-1,1,taun)+ t(i,j,k,1,taun))/2.0
        t(i,j,k,1,taun) = t(i,j,k-1,1,taun)
        endif
        endif
	enddo
	
#ifdef atmosphere
	u_atmos(i,j,taum) = u_atmos(i,j,taun)
        u_atmos(i,j,taun) = u_atmos(i,j,taup)
        v_atmos(i,j,taum) = v_atmos(i,j,taun)
        v_atmos(i,j,taun) = v_atmos(i,j,taup)
        h_atmos(i,j,taum) = h_atmos(i,j,taun)
        h_atmos(i,j,taun) = h_atmos(i,j,taup)
#endif
	
	enddo
	enddo
	enddo

#ifdef prescribe_S_boundary
        do i=1,imt
        do j=1,2
        do k=1,km
         u(i,j,k,taun) = 0.0 * rkmu(i,j) * rkmv(i,j)
        enddo
        enddo
        enddo
        do i=1,imt-1
        do j=1,2
        do k=1,km
         v(i,j,k,taun) = 0.0 * rkmv(i,j) * rkmu(i,j)
        enddo
        enddo
        enddo
        
#endif


#ifdef apply_spatial_filter
c
c 5-point smoother for u, v and h
c
	
	idim5 = 5  ! index for 5 point smoother, default =5
	idim4 = 0  ! index for t, if =4 it smooths
	idim3 = 3  ! index for h, if =3 it smooths
	idim2 = 2  ! index for v, if =2 it smooths
	idim1 = 1  ! index for u, if =1 it smooths
	
	do ii=1,imt
	do jj=1,jmt
	do kk=1,km
         u(ii,jj,kk,4) = u(ii,jj,kk,2)
         v(ii,jj,kk,4) = v(ii,jj,kk,2)
         h(ii,jj,kk,4) = h(ii,jj,kk,2)
         u(ii,jj,kk,3) = u(ii,jj,kk,1)
         v(ii,jj,kk,3) = v(ii,jj,kk,1)
         h(ii,jj,kk,3) = h(ii,jj,kk,1)
#ifdef trace
        do nt=1,nn
         t(ii,jj,kk,nt,4) = t(ii,jj,kk,nt,2)
         t(ii,jj,kk,nt,3) = t(ii,jj,kk,nt,1)
        enddo
#endif
        enddo
        enddo
        enddo

c        pvort(:,:,:) = 0.0
	do ii=2,imt-1
	do jj=2,jmt-1
	do kk=1,km-1
	
        t_switch = 1.0
        linear_switch = 1.0
        do iii=ii-1,ii+1
        do jjj=jj-1,jj+1
         if(rkmu(iii,jjj) .eq. 0.0) 
     & linear_switch = 1.0 *mod(loop,int(86400/dt)*1)
         if(rkmv(iii,jjj) .eq. 0.0) 
     & linear_switch = 1.0 *mod(loop,int(86400/dt)*1)
        enddo
        enddo
        t_switch = 1.0 *mod(loop,int(86400/dt)*10)

	
	 if (abs(u(ii,jj,kk,taun)) .ge. vel_lim(kk) .or. 
     &    abs(v(ii,jj,kk,taun))  .ge. vel_lim(kk) .or.
     &    linear_switch .eq. 0 .or. t_switch .eq. 0)  then


	  imts = ii-1
	  imte = ii+1
	  jmts = jj-1
	  jmte = jj+1
	call smooth (idim1, idim2, idim3, idim4, idim5,
     &  imts, jmts, kk, imte, jmte, kk)
	 endif

        enddo
        enddo
        enddo


#endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef apply_spatial_filter_mld
c
c 5-point smoother for u, v and h
c
	
	idim5 = 5  ! index for 5 point smoother, default =5
	idim4 = 4  ! index for s, if =4 it smooths
	idim3 = 3  ! index for t, if =3 it smooths
	idim2 = 2  ! index for v, if =2 it smooths
	idim1 = 1  ! index for u, if =1 it smooths
	
        kmax = 51

        if (mod (loop, 10000000) .eq. 0) then

	do ii=1,imt
	do jj=1,jmt
	do kk=1,kmax
         temp(ii,jj,kk,2) = 0.0
         salt(ii,jj,kk,2) = 0.0
         uvel(ii,jj,kk,2) = 0.0
         vvel(ii,jj,kk,2) = 0.0
         tmask(ii,jj) = 1.0
         smask(ii,jj) = 1.0
         umask(ii,jj) = 1.0
         vmask(ii,jj) = 1.0
        if(temp(ii,jj,kk,1).eq. 0.0)tmask(ii,jj) = 0.0
        if(salt(ii,jj,kk,1).eq. 0.0)smask(ii,jj) = 0.0
        if(uvel(ii,jj,kk,1).eq. 0.0)umask(ii,jj) = 0.0
        if(vvel(ii,jj,kk,1).eq. 0.0)vmask(ii,jj) = 0.0
        enddo
        enddo
        enddo

	do ii=2,imt-1
	do jj=2,jmt-1
	
	  imts = ii-1
	  imte = ii+1
	  jmts = jj-1
	  jmte = jj+1

	call smooth_mld(idim3, idim4, idim5,
     &  imts, jmts, 1, imte, jmte, kmax)

	enddo
	enddo

	do ii=2,imt-1
	do jj=2,jmt-1
	do kk=1,kmax
         if (temp(ii,jj,kk,2) .ne. 0.0 ) 
     & temp(ii,jj,kk,1) = temp(ii,jj,kk,2) 
         if (salt(ii,jj,kk,2) .ne. 0.0 ) 
     & salt(ii,jj,kk,1) = salt(ii,jj,kk,2)
        enddo
        enddo
        enddo

        endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif


	return
	
	end subroutine filter
