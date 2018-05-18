        subroutine output_manager (no_of_snap)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Average the output over the sampling period and write output
c
c
c
c
c       Generated on : 04-June-2017     author: valsala@tropmet.res.in
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef output_average
            use size_mod
            use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

        do i=1,imt
        do j=1,jmt
        do k=1,km
        h_output(i,k,j) = h_output(i,k,j) + h(i,k,j,taun)*rkmh(i,j)
        u_output(i,k,j) = u_output(i,k,j) + u(i,k,j,taun)*rkmu(i,j)
        v_output(i,k,j) = v_output(i,k,j) + v(i,k,j,taun)*rkmv(i,j)
        eta_output(i,k,j) = eta_output(i,k,j) + eta(i,k,j,taun)
        tr01_output(i,k,j) = tr01_output(i,k,j) + t(i,k,j,1,taun)
        tr02_output(i,k,j) = tr02_output(i,k,j) + t(i,k,j,2,taun)
        we_output(i,k,j) = we_output(i,k,j) + we(i,k,j)
        pvort_output(i,k,j) = pvort_output(i,k,j) + pvort(i,k,j)
        dens_output(i,k,j) = dens_output(i,k,j) + denss(i,k,j)
        enddo
        do k=1,51
        temp_output(i,k,j) = temp_output(i,k,j) + temp(i,k,j,1)
        salt_output(i,k,j) = salt_output(i,k,j) + salt(i,k,j,1)
        uvel_output(i,k,j) = uvel_output(i,k,j) + uvel(i,k,j,1)
        vvel_output(i,k,j) = vvel_output(i,k,j) + vvel(i,k,j,1)
        SMCoeff_output(i,k,j) = SMCoeff_output(i,k,j) + SMCoeff(i,k,j)
        SHCoeff_output(i,k,j) = SHCoeff_output(i,k,j) + SHCoeff(i,k,j)
        rmld_misc_output(i,k,j) = rmld_misc_output(i,k,j) + rmld_misc(i,k,j)
#ifdef snap_mld_extended
        diag_ext1_output(i,k,j)=diag_ext1_output(i,k,j)+diag_ext1(i,k,j)
        diag_ext2_output(i,k,j)=diag_ext2_output(i,k,j)+diag_ext2(i,k,j)
        diag_ext3_output(i,k,j)=diag_ext3_output(i,k,j)+diag_ext3(i,k,j)
        diag_ext4_output(i,k,j)=diag_ext4_output(i,k,j)+diag_ext4(i,k,j)
        diag_ext5_output(i,k,j)=diag_ext5_output(i,k,j)+diag_ext5(i,k,j)
        diag_ext6_output(i,k,j)=diag_ext6_output(i,k,j)+diag_ext6(i,k,j)
#endif
        enddo
        tx_output(i,j) = tx_output(i,j) + taux_force(i,j)
        ty_output(i,j) = ty_output(i,j) + tauy_force(i,j)
        enddo
        enddo
        
        if ( mod (loop,loop_total/no_of_snap) .eq. 0) then

        h_output(:,:,:) = h_output(:,:,:) / (loop_total/no_of_snap)
        u_output(:,:,:) = u_output(:,:,:) / (loop_total/no_of_snap)
        v_output(:,:,:) = v_output(:,:,:) / (loop_total/no_of_snap)
        eta_output(:,:,:) = eta_output(:,:,:) / (loop_total/no_of_snap)
        tr01_output(:,:,:) = tr01_output(:,:,:) /(loop_total/no_of_snap)
        tr02_output(:,:,:) = tr02_output(:,:,:) /(loop_total/no_of_snap)
        we_output(:,:,:) = we_output(:,:,:) / (loop_total/no_of_snap)
        pvort_output(:,:,:) = pvort_output(:,:,:) /(loop_total/no_of_snap)
        dens_output(:,:,:) = dens_output(:,:,:)/(loop_total/no_of_snap)
        tx_output(:,:) = tx_output(:,:) / (loop_total/no_of_snap)
        ty_output(:,:) = ty_output(:,:) / (loop_total/no_of_snap)

        temp_output(:,:,:) = temp_output(:,:,:)/(loop_total/no_of_snap)
        salt_output(:,:,:) = salt_output(:,:,:)/(loop_total/no_of_snap)
        uvel_output(:,:,:) = uvel_output(:,:,:)/(loop_total/no_of_snap)
        vvel_output(:,:,:) = vvel_output(:,:,:)/(loop_total/no_of_snap)
        SMCoeff_output(:,:,:)=SMCoeff_output(:,:,:)/(loop_total/no_of_snap)
        SHCoeff_output(:,:,:)=SHCoeff_output(:,:,:)/(loop_total/no_of_snap)
        rmld_misc_output(:,:,:)=
     &             rmld_misc_output(:,:,:)/(loop_total/no_of_snap)
#ifdef snap_mld_extended
        diag_ext1_output(:,:,:) =
     &   diag_ext1_output(:,:,:)/(loop_total/no_of_snap)
        diag_ext2_output(:,:,:) =
     &   diag_ext2_output(:,:,:)/(loop_total/no_of_snap)
        diag_ext3_output(:,:,:) =
     &   diag_ext3_output(:,:,:)/(loop_total/no_of_snap)
        diag_ext4_output(:,:,:) =
     &   diag_ext4_output(:,:,:)/(loop_total/no_of_snap)
        diag_ext5_output(:,:,:) =
     &   diag_ext5_output(:,:,:)/(loop_total/no_of_snap)
        diag_ext6_output(:,:,:) =
     &   diag_ext6_output(:,:,:)/(loop_total/no_of_snap)
#endif
        endif
                             
#endif

        return
        end

        subroutine set_averages_to_zero
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod
            use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

#ifdef output_average
        u_output(:,:,:) = 0.0
        v_output(:,:,:) = 0.0
        h_output(:,:,:) = 0.0
        eta_output(:,:,:) = 0.0
        tx_output(:,:) = 0.0
        ty_output(:,:) = 0.0
        tr01_output(:,:,:) = 0.0
        tr02_output(:,:,:) = 0.0
        dens_output(:,:,:) = 0.0
        we_output(:,:,:) = 0.0
        pvort_output(:,:,:) = 0.0
        temp_output(:,:,:) = 0.0
        salt_output(:,:,:) = 0.0
        uvel_output(:,:,:) = 0.0
        vvel_output(:,:,:) = 0.0
        SHCoeff_output(:,:,:) = 0.0
        SMCoeff_output(:,:,:) = 0.0
        rmld_misc_output(:,:,:) = 0.0
        diag_ext1_output(:,:,:) = 0.0
        diag_ext2_output(:,:,:) = 0.0
        diag_ext3_output(:,:,:) = 0.0
        diag_ext4_output(:,:,:) = 0.0
        diag_ext5_output(:,:,:) = 0.0
        diag_ext6_output(:,:,:) = 0.0

#endif

        return
        end
