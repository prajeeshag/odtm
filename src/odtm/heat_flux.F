       
            module heat_flux_mod

            contains
         
         subroutine heat_flux( hm, airt_in, sst, uwndin, vwndin, 
     &   cldin,qnet_out, qs, ql, qb, tx, ty, Lv)
! hm - humidity, airt- air temperature, c= fractional cloudiness
! qnet - to return net flux
        implicit none

        real  hm, airt_in, sst, uwndin, vwndin, c
        real  rho_a, cd, cp, Lv, c1
        real  emmis, tx,ty
! empirical coefficent
        real  cs, cl
        real  qs, ql, qm, qb, qnet_out
        real TT, TTt, ew,Pa
        real sigma
        real :: cldin, coe, del_t, dq, sb, tau
               
        rho_a= 1.175
        cd=2.8*10**(-3)
        cp=1004.
        Lv=2.5008991412E+06 
     &    -2.3691806479E+03 * sst
     &    +2.6776439436E-01 * sst * sst
     &    -8.1027544602E-03 * sst * sst * sst
     &    -2.0799346624E-05 * sst * sst * sst * sst

        c1= 0.55
        emmis=0.97      
        sigma= 5.6e-8
        cs= 1.0e-3
        cl=1.5e-3
        Pa =1e3 ! presure in mb
        
        tau= sqrt(uwndin**2+vwndin**2)
        del_t= airt_in - sst ! +ve- air is warmer-> heat in ocean

! DRAG VALUES FROM SMITH (1988), JGR, TABLE-2. FOR T > 20 DEGREES
! AND T < 20
        if (sst .ge. 20) then
        if (tau .ge.  0.0 .and. tau .lt.  2.0) cs = 2.07e-3
        if (tau .ge.  2.0 .and. tau .lt.  5.0) cs = 1.48e-3
        if (tau .ge.  5.0 .and. tau .lt. 10.0) cs = 1.20e-3
        if (tau .ge. 10.0 .and. tau .lt. 15.0) cs = 1.10e-3
        if (tau .ge. 15.0 .and. tau .lt. 20.0) cs = 1.06e-3
        if (tau .ge. 20.0 .and. tau .lt. 25.0) cs = 1.04e-3
        cl = cs * 1.20
        endif
        if (sst .lt. 20) then
        if (tau .ge.  0.0 .and. tau .lt.  2.0) cs = 1.96e-3
        if (tau .ge.  2.0 .and. tau .lt.  5.0) cs = 1.42e-3
        if (tau .ge.  5.0 .and. tau .lt. 10.0) cs = 1.17e-3
        if (tau .ge. 10.0 .and. tau .lt. 15.0) cs = 1.08e-3
        if (tau .ge. 15.0 .and. tau .lt. 20.0) cs = 1.05e-3
        if (tau .ge. 20.0 .and. tau .lt. 25.0) cs = 1.03e-3
        cl = cs * 1.20
        endif
        if (sst .ge. 20) then
        if (tau .ge.  0.0 .and. tau .lt.  2.0) cd = 1.68e-3
        if (tau .ge.  2.0 .and. tau .lt.  5.0) cd = 1.45e-3
        if (tau .ge.  5.0 .and. tau .lt. 10.0) cd = 1.57e-3
        if (tau .ge. 10.0 .and. tau .lt. 15.0) cd = 1.74e-3
        if (tau .ge. 15.0 .and. tau .lt. 20.0) cd = 1.93e-3
        if (tau .ge. 20.0 .and. tau .lt. 25.0) cd = 2.14e-3
        endif
        if (sst .lt. 20) then
        if (tau .ge.  0.0 .and. tau .lt.  2.0) cd = 1.61e-3
        if (tau .ge.  2.0 .and. tau .lt.  5.0) cd = 1.40e-3
        if (tau .ge.  5.0 .and. tau .lt. 10.0) cd = 1.52e-3
        if (tau .ge. 10.0 .and. tau .lt. 15.0) cd = 1.71e-3
        if (tau .ge. 15.0 .and. tau .lt. 20.0) cd = 1.91e-3
        if (tau .ge. 20.0 .and. tau .lt. 25.0) cd = 2.12e-3
        endif

! TX WIND STRESS
        tx= rho_a*cd*uwndin*abs(uwndin)
        ty= rho_a*cd*vwndin*abs(vwndin)

! SENSIBLE HEAT FLUX +ve-> heat in ocean
        qs= rho_a*cp*cs*tau*del_t         


! LATENT HEAT FLUX, minus of what is in paper
! +ve -> heat is in ocean
!        dq=hm*qsat(airt_in, Pa) - 0.98*qsat(sst, Pa)
        dq=hm - 0.98*qsat(sst, Pa)
        ql= rho_a*cl*Lv*tau*dq
! LONGWAVE RADIATION from the ocean surface, always to atmosphere hence
! -ve 
        sb= (sst+273.15)**3 ! Stephan boltzmann's law requires T in
!kelvin
        coe=(0.154*(1.0-c1*cldin)*(sst+273.15))-(4*del_t)
        qb=-1.0*emmis*sigma*sb*coe

        qb=-1.0*emmis*sigma*(sst+273.15)**4 +
     &      0.8*sigma*((airt_in+273.15)**4)

! KANTHA AND CLAYSON, Small scale processes, Page-425
        qb=-1.0*emmis*sigma*(sst+273.15)**4 + (1 - 0.045)*
     &      sigma*((airt_in+273.15)**4)* (1.0-0.261*exp((-7.77e-4)*
     &      (273.15 - (airt_in+273.15))**2))
! net outward flux = sensible heat + latent heat - longwave radiation
        qnet_out =  qs+ql+qb
 
        return
        end subroutine heat_flux

 
         function qsat(Ta,Pa)
        
            implicit none

! QSAT: computes specific humidity at saturation.
! q=QSAT(Ta) computes the specific humidity (kg/kg) at satuation at
! air temperature Ta (deg C). Dependence on air pressure, Pa, is
! small
! but is included as an optional input.
!
!    INPUT:   Ta - air temperature  [C]
!             Pa - (optional) pressure [mb]
!
!    OUTPUT:  q  - saturation specific humidity  [kg/kg]
     
! Version 1.0 used Tetens' formula for saturation vapor pressure
!! from Buck (1981), J. App. Meteor., 1527-1532.  This version
! follows the saturation specific humidity computation in the COARE
! Fortran code v2.5b.  This results in an increase of ~5% in
! latent heat flux compared to the calculation with version 1.0.
!     
!    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       3/8/97: version 1.0
!       4/7/99: version 1.2 (revised as above by AA)
!       8/5/99: version 2.0
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        real :: ew, pa, qsat, ta     
     
        ew = 6.1121*(1.0007+3.46e-6*Pa)*exp((17.502*Ta)/(240.97+Ta))
        qsat = 0.62197*(ew/(Pa-0.378*ew))  !% mb -> kg/kg        
   

        end function 
            end module heat_flux_mod 
