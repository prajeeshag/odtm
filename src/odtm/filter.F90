module filter_mod
    use size_mod, only : eta, h, i, j, rkmh, rkmu, rkmv, t, u, v, k
    use size_mod, only : loop, taum, taun, taup, km, nn
    use size_mod, only : isc, iec, jsc, jec, dau, dav, dah
    use size_mod, only : isd, ied, jsd, jed, temp, salt, kmaxMYM
    use param_mod, only : alpha, dt, day2sec
    use mpp_domains_mod, only : domain2d, mpp_update_domains
    use mpp_mod, only : mpp_error, NOTE, WARNING, FATAL, mpp_sum

    implicit none
    
    contains
            
    subroutine filter(domain)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c
    !c	an implimentation of Asselin-Robert Filter
    !c	
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        type(domain2d) :: domain
        integer :: nt, idim5, idim4, idim3, idim2, idim1, ii, jj, kk, iii
        integer :: imte, imts, jjj, jmte, jmts, ncount
        real :: t_switch, vel_lim(km)
        logical :: lmask(isd:ied,jsd:jed,km-1), linear_switch

        vel_lim(1:km) = 1.   
 
#ifdef apply_spatial_filter
        
        linear_switch = mod(loop,int(day2sec/dt)*10)==0

        lmask = .true.
        if (.not.linear_switch) then
            lmask = .false.
            do kk = 1, km-1 
                lmask(isc:iec,jsc:jec,kk) = abs(u(isc:iec,jsc:jec,kk,taup)) > vel_lim(kk) &
                                        .or.abs(v(isc:iec,jsc:jec,kk,taup)) > vel_lim(kk)
            enddo
            ncount = count(lmask)
            call mpp_sum(ncount)
        endif

        if (linear_switch .or. ncount > 0) then
            call mpp_update_domains(u(:,:,:,taup),domain)
            call mpp_update_domains(v(:,:,:,taup),domain)
            call smooth_hanning(u(:,:,:,taup),dmask=rkmu,mask=lmask,area=dau)
            call smooth_hanning(v(:,:,:,taup),dmask=rkmv,mask=lmask,area=dav)
        endif
        ! Added by Vinu 28-05-2018
        if (linear_switch .or. ncount > 0) then
            call mpp_update_domains(u(:,:,:,taun),domain)
            call mpp_update_domains(v(:,:,:,taun),domain)
            call smooth_hanning(u(:,:,:,taun),dmask=rkmu,mask=lmask,area=dau)
            call smooth_hanning(v(:,:,:,taun),dmask=rkmv,mask=lmask,area=dav)
        endif
        ! Added by Vinu 28-05-2018
        ! Added by Vinu 29-05-2018
        if (linear_switch ) then
            call mpp_update_domains(temp(:,:,:,1),domain)
            call mpp_update_domains(salt(:,:,:,1),domain)
            call smooth_hanning(temp(:,:,:,1),dmask=rkmh,area=dah)
            call smooth_hanning(salt(:,:,:,1),dmask=rkmh,area=dah)
        endif
        ! Added by Vinu 29-05-2018
#endif
    
        do i=isc, iec
            do j=jsc, jec
                do k=1,km-1
                    u(i,j,k,taun) = u(i,j,k,taun)+ alpha*rkmu(i,j) &
                     *(u(i,j,k,taup) -2*u(i,j,k,taun) + u(i,j,k,taum) )
             
                    v(i,j,k,taun) = v(i,j,k,taun)+ alpha*rkmv(i,j) & 
                    *(v(i,j,k,taup) -2*v(i,j,k,taun) + v(i,j,k,taum) )
    
                    do nt=1,nn
                        t(i,j,k,nt,taun) = t(i,j,k,nt,taun)+ alpha &
                        *(t(i,j,k,nt,taup) -2*t(i,j,k,nt,taun) + t(i,j,k,nt,taum) ) 
                    enddo 
    
                    h(i,j,k,taun) = h(i,j,k,taun)+ alpha*rkmh(i,j) &
                    *(h(i,j,k,taup) -2*h(i,j,k,taun) + h(i,j,k,taum) )
    
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
                                t(i,j,k,1,taum) = (t(i,j,k-1,1,taum) &
                                                + t(i,j,k,1,taum)) * 0.5
                                t(i,j,k-1,1,taum) = t(i,j,k,1,taum)
                            endif
                        endif
    
                        if (k.gt.1) then !Prajeesh  
                            if (t(i,j,k,1,taun).gt.t(i,j,k-1,1,taun))then
                                t(i,j,k,1,taun) = (t(i,j,k-1,1,taun) &
                                                + t(i,j,k,1,taun)) * 0.5
                                t(i,j,k,1,taun) = t(i,j,k-1,1,taun)
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo
    
#ifdef prescribe_S_boundary
        do i=isc, iec
            do j=1,2
                do k=1,km
                    u(i,j,k,taun) = 0.0 * rkmu(i,j) * rkmv(i,j)
                enddo
            enddo
        enddo
        do i=isc, iec
            do j=1,2
                do k=1,km
                    v(i,j,k,taun) = 0.0 * rkmv(i,j) * rkmu(i,j)
                enddo
            enddo
        enddo
#endif
    
        return
    end subroutine filter

    subroutine smooth_hanning(fld, dmask, mask, area)

        real, intent(inout) :: fld(:,:,:)
        real, intent(in) :: dmask(:,:)
        logical, intent(in), optional :: mask(:,:,:)
        real, intent(in), optional :: area(:,:)
    
        real :: fld1(size(fld,1),size(fld,2),size(fld,3)), rdiv1, rdiv2
        logical :: lmask(size(fld,1),size(fld,2),size(fld,3))
        integer :: is, ie, js, je, ks, ke
        integer :: i, j, k, im, jm, ip, jp

       
        is = 2; ie = size(fld,1) - 1
        js = 2; je = size(fld,2) - 1
        ks = 1; ke = size(fld,3)

        lmask = .true.
        if (present(mask)) lmask = mask

        if (all(.not.lmask)) return
        
        do k = ks, ke
            fld1(:,:,k) = fld(:,:,k) * dmask(:,:)
        enddo

        if (present(area)) then
            do k = ks, ke
                fld1(:,:,k) = fld1(:,:,k) * area(:,:)
            enddo
        endif
 
        do i = is, ie
            im = i - 1; ip = i + 1

            do j = js, je
                jm = j - 1; jp = j + 1

                rdiv1 = dmask(im,j)+dmask(ip,j)+dmask(i,jm)+dmask(i,jp)
                rdiv2 = dmask(im,jm)+dmask(ip,jm)+dmask(im,jp)+dmask(ip,jp)

                do k = ks, ke   
                    if (.not.lmask(i,j,k)) cycle
                    if (dmask(i,j)==0.) cycle
                    fld(i,j,k) = 0.25 * fld1(i,j,k) &
                                 + 0.5 * ( fld1(im,j,k) + fld1(ip,j,k) &
                                 + fld1(i,jm,k) + fld1(i,jp,k))/max(1.0,rdiv1) &
                                 + 0.25 * (fld1(im,jm,k) + fld1(ip,jm,k) &
                                 + fld1(im,jp,k) + fld1(ip,jp,k))/max(1.0,rdiv2)

                    if (present(area)) then
                        fld(i,j,k) = fld(i,j,k)/area(i,j)
                    endif

                enddo
            enddo
        enddo

    end subroutine smooth_hanning

end module filter_mod 
