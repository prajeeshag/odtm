module filter_mod
    use size_mod, only : eta, h, i, j, rkmh, rkmu, rkmv, t, u, v, k
    use size_mod, only : loop, taum, taun, taup, km, nn
    use size_mod, only : isc, iec, jsc, jec
    use param_mod, only : alpha, dt
    
    implicit none
    
    contains
            
    subroutine filter
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !c
    !c
    !c	an implimentation of Asselin-Robert Filter
    !c	
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit none
        integer :: nt, idim5, idim4, idim3, idim2, idim1, ii, jj, kk, iii
        integer :: imte, imts, jjj, jmte, jmts, linear_switch
        real :: t_switch
    
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
    
    
#ifdef apply_spatial_filter
    
       ! call smooth_hanning(u(:,:,1:km-1,taum),rkmu)
       ! call smooth_hanning(u(:,:,1:km-1,taun),rkmu)
       ! call smooth_hanning(v(:,:,1:km-1,taum),rkmv)
       ! call smooth_hanning(v(:,:,1:km-1,taun),rkmv)
       ! call smooth_hanning(h(:,:,1:km-1,taum),rkmh)
       ! call smooth_hanning(h(:,:,1:km-1,taun),rkmh)
    
#endif
    
        return
    end subroutine filter

    subroutine smooth_hanning(fld,mask,area)

        real, intent(inout) :: fld(:,:,:)
        real, intent(in) :: mask(:,:)
        real, intent(in), optional :: area(:,:)
    
        real :: fld1(size(fld,1),size(fld,2),size(fld,3)), rdiv1, rdiv2
        integer :: is, ie, js, je, ks, ke
        integer :: i, j, k, im, jm, ip, jp

        is = 2; ie = size(fld,1) - 1
        js = 2; je = size(fld,2) - 1
        ks = 1; ke = size(fld,3)

        do k = ks, ke
            fld1(:,:,k) = fld(:,:,k) * mask(:,:)
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
                rdiv1 = mask(im,j)+mask(ip,j)+mask(i,jm)+mask(i,jp)
                rdiv2 = mask(im,jm)+mask(ip,jm)+mask(im,jp)+mask(ip,jp)
                do k = ks, ke   
                   fld(i,j,k) = 0.25 * fld1(i,j,k) &
                                 + 0.5 * ( fld1(im,j,k) + fld1(ip,j,k) &
                                 + fld1(i,jm,k) + fld1(i,jp,k))/max(1.0,rdiv1) &
                                 + 0.25 * (fld1(im,jm,k) + fld1(ip,jm,k) &
                                 + fld1(im,jp,k) + fld1(ip,jp,k))/max(1.0,rdiv2)
                enddo
            enddo
        enddo

        if (present(area)) then
            do k = ks, ke
                fld(is:ie,js:je,k) = fld(is:ie,js:je,k)/area(is:ie,js:je)
            end do
        endif

    end subroutine smooth_hanning

end module filter_mod 
