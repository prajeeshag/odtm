#ifdef smagorinsky_biharmonic
	subroutine biharmon
        use size_mod
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	an attempt to impliment bi-harmonic friciton from the DIFF_flux
c	of smagaronsky diffusion scheme.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!#include "size.h"
!#include "polar.h"
#include "param.h"
	
	rdel = min(dx,dy)
c =========================== B.d4u.dx4 ==================================
	rtemp1 = sqrt(Asmag(i,k,j))*(DIFFuu(i+1,k,j)-DIFFuu(i,k,j))/dx
	rtemp2 = sqrt(Asmag(i-1,k,j))*(DIFFuu(i,k,j)-DIFFuu(i-1,k,j))/dx
	rtemp3 = 1/dx
c ============================ B.d4u/dy4 =================================
	rtemp4=0.25*(sqrt(Asmag(i,k,j))+sqrt(Asmag(i,k,j+1))+
     &   sqrt(Asmag(i-1,k,j+1))
     &   +sqrt(Asmag(i-1,k,j)))/dy
	
	rtemp5 = DIFFvu(i,k,j+1)-DIFFvu(i,k,j)
	rtemp6=0.25*(sqrt(Asmag(i,k,j-1))+
     &     sqrt(Asmag(i,k,j))+sqrt(Asmag(i-1,k,j))
     &     +sqrt(Asmag(i-1,k,j-1)))/dy
	rtemp7 = DIFFvu(i,k,j)-DIFFvu(i,k,j-1)
	rtemp8 = 1/dy
	
	biharmon_u = rtemp3*(rtemp1-rtemp2)*rdel/sqrt(8.0)
     &      + (rtemp4*rtemp5 - rtemp6*rtemp7)*rtemp8*
     &     rdel/sqrt(8.0)

c =========================== A.d2v/dx2 =================================
	rtemp1 = 0.25*(sqrt(Asmag(i,k,j))+
     &    sqrt(Asmag(i+1,k,j))+sqrt(Asmag(i+1,k,j-1))
     &    +sqrt(Asmag(i,k,j-1)))/dx
	rtemp2 = DIFFuv(i+1,k,j)-DIFFuv(i,k,j)
	rtemp3 = 0.25*(sqrt(Asmag(i-1,k,j))+
     &    sqrt(Asmag(i,k,j))+sqrt(Asmag(i,k,j-1))
     &    +sqrt(Asmag(i-1,k,j-1)))/dx
	rtemp4 = DIFFuv(i,k,j)-DIFFuv(i-1,k,j)
	rtemp5 = 1/dx

	rtemp6 =sqrt(Asmag(i,k,j))*(DIFFvv(i,k,j+1)-DIFFvv(i,k,j))/dy
	rtemp7 =sqrt(Asmag(i,k,j-1))*(DIFFvv(i,k,j)-DIFFvv(i,k,j-1))/dy
	rtemp8 = 1/dy
	
	biharmon_v = ((rtemp1*rtemp2 - rtemp3*rtemp4)*rtemp5*
     &   rdel/sqrt(8.0))
     &   + (rtemp6 - rtemp7)*rtemp8*rdel/sqrt(8.0)
	

	return
	stop
	end

#endif
