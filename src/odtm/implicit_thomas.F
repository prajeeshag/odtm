        module implicit_thomas_mod

        contains 

        subroutine implicit_vmix (KH, dz_in, dt_in, km, t1, t2)
        
        implicit none
    
        integer :: i, k, km

        real t1(km), t2(km), KH(km)
	real t(0:km+1), dz(0:km+1)
	real dtts, tmask(km)
	real rD1(km),  rD2(km)
	real rDDu(km), rDDl(km)
        real dzfact(km,km), dzfact_rev(km,km)
        real rheat1, rheat2
        real a(km), b(km), c(km), d(km)
        real bt_temp
        real diffuse(0:km+1)
        real dz_in, dt_in
        real dz1, dt
        
	
c
c initialize data
c
        do k=1,km
        t2(k) = t1(k)
        enddo

         dz1 = dz_in
         dt  = dt_in

         dtts = dt
        
         do i=1,km      
         dz(i) = dz1
         t(i) = t1(i)
         diffuse(i) = KH(i) 
         enddo
        
         diffuse(0) = 0
         diffuse(km+1) = 0
         dz(0) = dz(1)
         dz(km+1) = dz(km)

	 do i=1,km
	 rDDu(i) =  0.5*dtts/((dz(i-0)+dz(i)))
	 rDDl(i) =  0.5*dtts/((dz(i)+dz(i+0)))
         if (i-1 .eq. 0) then
         rDDu(i) = 0.5*dtts /(dz(i-0)+dz(i))
         rDDl(i) = 0.5*dtts /(dz(i)+dz(i+0))
         endif
         rD1(i) = diffuse(i-1)/((dz(i-1)+dz(i))*0.5)
         rD2(i) = diffuse(i)/((dz(i)+dz(i+1))*0.5)
	 enddo
	
         rD1(1) = 0.0
         rD2(km) = 0.0
	
c
c form the matrix
c	        
         t(0) = 0
         t(km+1) = 0
	 do k=1,km
         a(k) = -rDDu(k)*rD1(k)
         b(k) = (1.0 + (rDDu(k)*rD1(k)+rDDl(k)*rD2(k)) )
         bt_temp = (1.0 - (rDDu(k)*rD1(k)+rDDl(k)*rD2(k)) )
         c(k) = -rDDl(k)*rD2(k)
         d(k) = -c(k)*t(k+1) + bt_temp*t(k) - a(k)*t(k-1)

         if (abs(b(k)) .lt. abs(a(k)) + abs(c(k)) ) goto 101

        enddo


         call thomas(km, a, b, c, d)
        
        do k=1,km
        t2(k) = d(k)
        enddo

 101    continue

        return
	end subroutine 
	
        subroutine thomas (km, a,b,c,d)

        implicit none
        integer km, j

        real, intent(inout) :: a(km),b(km)
        real, intent(inout) :: c(km),d(km)

!combined decomposition and forward substitution
        do j = 2,km 
            a(j) = a(j)/b(j-1)
            b(j) = b(j)-a(j)*c(j-1)
            d(j) = d(j)-a(j)*d(j-1)
          end do

!back substitution
          d(km) = d(km)/b(km)
          do j = km-1,1,-1
            d(j) = (d(j)-c(j)*d(j+1))/b(j)
          end do
        
        return
        end subroutine

        end module implicit_thomas_mod 
