
subroutine clinic(domain)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c   subroutine to solve baroclinic pressure gradient
!c
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use size_mod, only : dz, dz_max, gdx, gdy, h, i, j, loop, dz_min
    use size_mod, only : rdx, rkmh, rdy, taum1, taun1, taup1, u, v, wd
    use size_mod, only : we, we_mld, k, taum, taun, taup, theta, rdxh, rdyh
    use size_mod, only : isc, iec, jsc, jec, km
    use size_mod, only : isd, ied, jsd, jed
    use param_mod, only : day2sec, diffuse_h, diffuse_th, dphi, dt, dthe
    use param_mod, only : dtts, pi, re, tropdubdx, tropdvbdy

    use advec_mod, only : sum_2pt
    use mpp_domains_mod, only : mpp_update_domains, domain2d

    implicit none
    type(domain2d), intent(inout) :: domain
    real tx,ty, Lv
    integer :: ip, im, jp, jm,ii, jj
    real :: rtemp1, rtemp2, rtheta_vu, rtheta_vd, rtemp3, rtemp4, rmask_right
    real :: rmask_left, rmask_top,rmask_bot, rthetapone, rthetamone, rtemp5
    real :: rtemp6, x0, y0, A0, rexp, rsig, rsig1, dtts_back 
    logical :: NaN_check

    taum1 = 1
    taun1 = 2
    taup1 = 3
    
    do i=isc,iec
        do j=jsc,jec
            if (rkmh(i,j) == 0.) cycle
            do k=1,km-1
                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1

                tropdubdx = 0.0
                tropdvbdy = 0.0
    
                rtemp1 = u(ip,j,k,taun)* sum_2pt (3,ip,j,i,j) / sum_2pt (31,ip,j,i,j)
                rtemp2 = u(i,j,k,taun)*  sum_2pt (3,im,j,i,j) / sum_2pt (31,im,j,i,j)
                dphi = (rdx(ip)+rdx(i))/2.0
                tropdubdx = (rtemp1 - rtemp2)*rdxh(i,j)

                rtheta_vu = cos((theta(j)+theta(j+1))*0.5)
                rtheta_vd = cos((theta(j)+theta(j-1))*0.5)
                dthe = (rdy(j) + rdy(j+1))/2.0

                rtemp3 = v(i,jp,k,taun) * sum_2pt (3,i,j,i,jp) *rtheta_vu / sum_2pt (31,i,j,i,jp)
                rtemp4 = v(i,j,k,taun) * sum_2pt (3,i,j,i,jm)  *rtheta_vd/ sum_2pt (31,i,j,i,jm)
                tropdvbdy = ( rtemp3 - rtemp4 ) * rdyh(i,j) / cos(theta(j))

                if (diffuse_th/=0) then
                    rmask_right = 1.0
                    rmask_left = 1.0
                    rmask_top = 1.0
                    rmask_bot = 1.0
                    
                    rthetapone = (theta(j+1)+theta(j) )*0.5
                    rthetamone = (theta(j-1)+theta(j) )*0.5
                    rtemp1 = Re*Re*cos(theta(j))*cos(theta(j)) 
                    rtemp1 = 1.0/rtemp1
                     if (rkmh(ip,j) .eq. 0.0) rmask_right = 0.0
                     if (rkmh(im,j) .eq. 0.0) rmask_left = 0.0 !Prajeesh   
                     if (rkmh(i,j+1) .eq. 0.0) rmask_top = 0.0
                     if (rkmh(i,jm) .eq. 0.0) rmask_bot = 0.0 !Prajeesh 
 
                    dphi = (rdx(ip)+rdx(i))/2.0
                    rtemp2 = rmask_right*(h(ip,j,k,taum) - h(i,j,k,taum) )/dphi
                    dphi = (rdx(i)+rdx(i-1))/2.0
                    rtemp3 = rmask_left*(h(i,j,k,taum) - h(im,j,k,taum) )/dphi
                    dphi = rdx(i)
                    rtemp4 = (rtemp2 - rtemp3)/dphi
                    rtemp5 = rtemp1*rtemp4  ! reserve
                    rtemp1 = Re*Re*cos(theta(j))
                    rtemp1 = 1.0/rtemp1
                    dthe = (rdy(j+1) + rdy(j))/2.0
                    rtemp2 = rmask_top*(h(i,j+1,k,taun) - h(i,j,k,taun))/dthe
                    dthe = (rdy(j-1) + rdy(j))/2.0
                    rtemp3 = rmask_bot*(h(i,j,k,taun) - h(i,jm,k,taun))/dthe
                    dthe = rdy(j)
                    rtemp4 = (cos(rthetapone)*rtemp2 - cos(rthetamone)*rtemp3)/dthe
                    rtemp6 = rtemp1*rtemp4  ! reserve
                    
                    diffuse_h = rtemp5 + rtemp6
                endif

                dtts_back = dtts
        
                we_mld(i,j,k) = (tropdubdx + tropdvbdy) + we(i,j,k) 

                h(i,j,k,taup) = h(i,j,k,taum)*rkmh(i,j) + ( &  
                              - (tropdubdx + tropdvbdy) &
                              +  we(i,j,k)*rkmh(i,j)*dtts_back/dtts &
                              +  wd(i,j,k)*rkmh(i,j)*dtts_back/dtts &
                              +  diffuse_th * diffuse_h &
                              ) *dtts*rkmh(i,j)

            enddo
        enddo
    enddo

    call mpp_update_domains(h(:,:,:,taup),domain)

    do i = isd, ied
        do j = jsd, jed
            if (rkmh(i,j)==0.) cycle
            do k = 1, km-1
                im = i - 1
                if (im<isc) im=isc
                jm = j - 1
                if (jm<isc) jm=jsc
                ip = i + 1
                if (ip>iec) ip=iec
                jp = j + 1
                if (jp>jec) jp=jec
                
                if ( h(i,j,k,taup) .gt. dz_max(k)) then
                    do ii=im,ip
                        do jj=jm,jp
                            h(ii,jj,k,taup) =  h(ii,jj,k,taup) &
                            - (10**(-1.5))*  h(ii,jj,k,taun)*rkmh(ii,jj)
                        enddo
                    enddo
                endif

                if ( h(i,j,k,taup) .lt. dz_min(k)) then
                    do ii=im,ip
                        do jj=jm,jp
                            h(ii,jj,k,taup) =  h(ii,jj,k,taup) &
                            + (10**(-1.5))* h(ii,jj,k,taun)*rkmh(ii,jj)
                        enddo
                    enddo
                endif
            enddo
        enddo
    enddo

    call mpp_update_domains(h(:,:,:,taup),domain)

    return
end subroutine clinic
