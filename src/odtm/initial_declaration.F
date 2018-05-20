	subroutine initial_declaration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	Declare all initial values related to arrays.
c	All maximum and minimum limits are declared here.
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : dz, k, rho, rho_max, rho_min, zdz
            use size_mod, only : tr01_max, tr01_min, tr02_max, tr02_min
            use size_mod, only : dz_max, hd, he, dz_min, km
        implicit none
        real :: rsum
!c      depth of the active layer
        
        dz(1) = 50.0     ! layer-1
        dz(2) = 25.0     ! layer-2
        dz(3) = 25.0     ! layer-3
        dz(4) = 25.0
        dz(5) = 25.0
        dz(6) = 25.0
        dz(7) = 25.0
        dz(8) = 25.0
        dz(9) = 25.0
        dz(10) = 50.0
        dz(11) = 100.0
        dz(12) = 2400.0

        zdz(1) = dz(1) - dz(1)/2.0
        rsum = dz(1)
        do k=2,km
        rsum = rsum + dz(k)
        zdz(k) = rsum - dz(k)/2.0
        enddo


        rho(1) = 1022.790082
        rho(2) = 1023.474698
        rho(3) = 1024.153075
        rho(4) = 1024.775912
        rho(5) = 1025.279050
        rho(6) = 1025.647721
        rho(7) = 1025.922296
        rho(8) = 1026.129476
        rho(9) = 1026.287885
        rho(10) = 1026.461588
        rho(11) = 1026.679178
        rho(12) = 1027.287534



        rho_max(1) = 1028.354643
        rho_max(2) = 1028.272647
        rho_max(3) = 1028.218280
        rho_max(4) = 1028.346995
        rho_max(5) = 1028.388880
        rho_max(6) = 1028.457778
        rho_max(7) = 1028.507549
        rho_max(8) = 1028.545115
        rho_max(9) = 1028.555075
        rho_max(10) = 1028.560275
        rho_max(11) = 1028.583104
        rho_max(12) = 1028.607435

	rho_min(1) = 1017.487049
        rho_min(2) = 1020.429045
        rho_min(3) = 1021.431419
        rho_min(4) = 1022.281136
        rho_min(5) = 1023.290114
        rho_min(6) = 1023.582514
        rho_min(7) = 1024.102586
        rho_min(8) = 1024.76888
        rho_min(9) = 1025.21475
        rho_min(10) = 1025.531348
        rho_min(11) = 1025.985107
        rho_min(12) = 1026.234502

	

	tr01_max(1) = 34.93
        tr01_max(2) = 30.30
        tr01_max(3) = 29.40
        tr01_max(4) = 27.90
        tr01_max(5) = 25.53
        tr01_max(6) = 24.56
        tr01_max(7) = 22.90
        tr01_max(8) = 22.17
        tr01_max(9) = 22.06
        tr01_max(10) = 21.93
        tr01_max(11) = 21.79
        tr01_max(12) = 21.77



        tr01_min(1) = 16.78
        tr01_min(2) = 15.79
        tr01_min(3) = 14.62
        tr01_min(4) = 13.39
        tr01_min(5) = 12.68
        tr01_min(6) = 12.06
        tr01_min(7) = 11.24
        tr01_min(8) = 10.94
        tr01_min(9) = 10.53
        tr01_min(10) = 09.79
        tr01_min(11) = 08.21
        tr01_min(12) = 02.71

	tr02_max(1) = 40.32
        tr02_max(2) = 40.27
        tr02_max(3) = 40.38
        tr02_max(4) = 40.41
        tr02_max(5) = 40.49
        tr02_max(6) = 40.53
        tr02_max(7) = 40.56
        tr02_max(8) = 40.58
        tr02_max(9) = 40.58
        tr02_max(10) = 40.58
        tr02_max(11) = 40.59
        tr02_max(12) = 40.61




        tr02_min(1) = 28.83
        tr02_min(2) = 32.52
        tr02_min(3) = 33.32
        tr02_min(4) = 33.95
        tr02_min(5) = 34.26
        tr02_min(6) = 34.39
        tr02_min(7) = 34.41
        tr02_min(8) = 34.42
        tr02_min(9) = 34.38
        tr02_min(10) = 34.38
        tr02_min(11) = 34.37
        tr02_min(12) = 34.36


! he is the threshold height of layer at and below which entrainment starts

        he(1) =  35.0
        he(2) =  15.0
        he(3) =  15.0
        he(4) =  15.0
        he(5) =  15.0
        he(6) =  15.0
        he(7) =  15.0
        he(8) =  15.0
        he(9) =  15.0
        he(10) =  35.0
        he(11) =  50.0
        he(12) = 800.0

! hd is the ht of layer, h>=hd- detrainment present
        hd(1) = 100.0! what is hd?
        hd(2) = 50.0
        hd(3) = 50.0
        hd(4) = 50.0
        hd(5) = 50.0
        hd(6) = 50.0
        hd(7) = 50.0
        hd(8) = 50.0
        hd(9) = 50.0
        hd(10) = 100.0
        hd(11) = 200.0
        hd(12) = 1600.0

        dz_max(1) = 3*dz(1) ! what is dz_max?
        dz_max(2) = 3*dz(2)
        dz_max(3) = 3*dz(3)
        dz_max(4) = 3*dz(4)
        dz_max(5) = 3*dz(5)
        dz_max(6) = 3*dz(6)
        dz_max(7) = 3*dz(7)
        dz_max(8) = 3*dz(8)
        dz_max(9) = 3*dz(9)
        dz_max(10) = 3*dz(10)
        dz_max(11) = 3*dz(11)
        dz_max(12) = 3*dz(12)

        dz_min(1) = 20 
        dz_min(2) = 10 
        dz_min(3) = 10
        dz_min(4) = 10
        dz_min(5) = 10
        dz_min(6) = 10
        dz_min(7) = 10
        dz_min(8) = 10
        dz_min(9) = 10
        dz_min(10) = 20
        dz_min(11) = 20
        dz_min(12) = 500



        return
        end
