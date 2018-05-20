#ifdef cyclic_boundary
	subroutine cyclic
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c	attempt to impliment cyclic boundary condition
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include "size.h"
!#include "polar.h"
#include "param.h"

	do j=2,9
        do k=1,km
	
	 u(1,k,j,taus) = u(imt-1,k,j,taup)*rkmu(1,j)
	 v(1,k,j,taus) = v(imt-1,k,j,taup)*rkmv(1,j)
	 h(1,k,j,taus) = h(imt-1,k,j,taup)*rkmt(1,j)
	 eta(1,k,j,taus) = eta(imt-1,k,j,taup)*rkmt(1,j)
	 
	 
	 u(imt,k,j,taus) = u(2,k,j,taup)*rkmu(imt,j)
	 v(imt,k,j,taus) = v(2,k,j,taup)*rkmv(imt,j)
	 h(imt,k,j,taus) = h(2,k,j,taup)*rkmt(imt,j)
	 eta(imt,k,j,taus) = eta(2,k,j,taup)*rkmt(imt,j)
	 enddo
	 enddo

	return
	stop
	end
#endif
