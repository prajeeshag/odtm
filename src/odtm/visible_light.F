        function rvisible_extinct (rchla, kk)
        
        implicit none
c
c table of coefficients from Morel and Amtoine (1994)
c J. Physical Oceanography, Vol. 24, 1652-1664
c to find the attenuation of PAR from chlorophyll-a
c concentrations.
c
        real :: rvisible_extinct

        real chlz(51)

       real  rVA (2,6), rZA(2,6) 

        real :: chl_X, chl_X2, chl_X3, chl_X4, chl_X5, rchla
        real :: rtempV1, rtempV2, rtempZ1, rtempZ2, rtemp1, rtemp2
        integer :: kk


        chlz(1) =  1.08958944
        chlz(2) =  1.29577000
        chlz(3) =  1.52888000
        chlz(4) =  1.93672000
        chlz(5) =  2.54879000
        chlz(6) =  3.46583000
        chlz(7) =  4.71489000
        chlz(8) =  6.49442000
        chlz(9) =  8.71933000
        chlz(10) = 10.82823000
        chlz(11) = 12.50114000
        chlz(12) = 13.54078000
        chlz(13) = 13.89062000
        chlz(14) = 13.60208000
        chlz(15) = 12.42530000
        chlz(16) = 10.99515000
        chlz(17) =  9.02130000
        chlz(18) =  6.94887000
        chlz(19) =  4.90519000
        chlz(20) =  3.64370000
        chlz(21) =  2.00979000
        chlz(22) =  1.16200000
        chlz(23) =  0.54000000
        chlz(24) =  0.54000000
        chlz(25) =  0.38880000
        chlz(26) =  0.14000000
        chlz(27) =  0.14000000
        chlz(28) =  0.14000000
        chlz(29) =  0.10220000
        chlz(30) =  0.04000000
        chlz(31) =  0.04000000
        chlz(32) =  0.04000000
        chlz(33) =  0.04000000
        chlz(34) =  0.04000000
        chlz(35) =  0.04000000
        chlz(36) =  0.03622000
        chlz(37) =  0.03000000
        chlz(38) =  0.03000000
        chlz(39) =  0.03000000
        chlz(40) =  0.03000000
        chlz(41) =  0.03000000
        chlz(42) =  0.03000000
        chlz(43) =  0.03000000
        chlz(44) =  0.03000000
        chlz(45) =  0.03000000
        chlz(46) =  0.02122000
        chlz(47) =  0.01000000
        chlz(48) =  0.01000000
        chlz(49) =  0.01000000
        chlz(50) =  0.01000000
        chlz(51) =  0.01000000


        rVA(1,1) =  0.353
        rVA(1,2) = -0.047
        rVA(1,3) =  0.083
        rVA(1,4) =  0.047
        rVA(1,5) = -0.011
        rVA(1,6) = -0.009

        rVA(2,1) =  0.647
        rVA(2,2) =  0.047
        rVA(2,3) = -0.083
        rVA(2,4) = -0.047
        rVA(2,5) =  0.011
        rVA(2,6) =  0.009

        rZA(1,1) =  1.662
        rZA(1,2) = -0.605
        rZA(1,3) =  0.128
        rZA(1,4) = -0.033
        rZA(1,5) = -0.051
        rZA(1,6) = -0.004

        rZA(2,1) =  8.541
        rZA(2,2) = -8.924
        rZA(2,3) =  4.020
        rZA(2,4) = -0.077
        rZA(2,5) = -0.536
        rZA(2,6) =  0.055
 
c
c special treatment for solar radiation damping with depth
c
        ! set a minimum chlorophyl value

        !Cution, above unit for chlorophyll_a is mg/m3

         chl_X = log10(rchla)
         chl_X2= chl_X**2
         chl_X3= chl_X**3
         chl_X4= chl_X**4
         chl_X5= chl_X**5

         rtempV1 = rVA(1,1) + rVA(1,2)*chl_X + rVA(1,3)*chl_X2
     &                      + rVA(1,4)*chl_X3+ rVA(1,5)*chl_X4
     &                      + rVA(1,6)*chl_X5

         rtempV2 = rVA(2,1) + rVA(2,2)*chl_X + rVA(2,3)*chl_X2
     &                      + rVA(2,4)*chl_X3+ rVA(2,5)*chl_X4
     &                      + rVA(2,6)*chl_X5

         rtempZ1 = rZA(1,1) + rZA(1,2)*chl_X + rZA(1,3)*chl_X2
     &                      + rZA(1,4)*chl_X3+ rZA(1,5)*chl_X4
     &                      + rZA(1,6)*chl_X5

         rtempZ2 = rZA(2,1) + rZA(2,2)*chl_X + rZA(2,3)*chl_X2
     &                      + rZA(2,4)*chl_X3+ rZA(2,5)*chl_X4
     &                      + rZA(2,6)*chl_X5

         rtemp1 = rtempV1*exp(-1.0d0*kk*5.0/rtempZ1)
         rtemp2 = rtempV2*exp(-1.0d0*kk*5.0/rtempZ2)
         rvisible_extinct = rtemp1 + rtemp2

        end
