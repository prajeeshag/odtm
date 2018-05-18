	subroutine restart
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	a subroutine to write restart file for the consicutive runs.
c	format is binary.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod
            use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

	open (11, file='restart.swm',form='unformatted',
     &   access='sequential')

c	write (11) imt, jmt, dxd, dyd, dt, diffuse, diffuse_th,
c     &  diffuse_ah, tracer_switch
c	write (11) loop
c	write (11) month
c	write (11) lpm
c	write (11) lpd
c	write (11) month_wind
c	write (11) iday_wind
c	write (11) month_start_snap
	
	do i=1,imt
	write (11) ((u(i,k,j,taun),j=1,jmt),k=1,km) !ok
	write (11) ((v(i,k,j,taun),j=1,jmt),k=1,km) !ok
	write (11) ((h(i,k,j,taun),j=1,jmt),k=1,km) !ok
#ifdef trace
	do nt=1,nn
	write (11) ((t(i,k,j,nt,taun),j=1,jmt),k=1,km) !ok
	enddo
#endif
	write (11) ((u(i,k,j,taum),j=1,jmt),k=1,km)
	write (11) ((v(i,k,j,taum),j=1,jmt),k=1,km)
	write (11) ((h(i,k,j,taum),j=1,jmt),k=1,km)
#ifdef trace
	do nt=1,nn
	write (11) ((t(i,k,j,nt,taum),j=1,jmt),k=1,km)
	enddo
#endif
#ifdef trace
	write (11) ((temp(i,k,j,taun),j=1,jmt),k=1,kmaxMYM)
	write (11) ((salt(i,k,j,taun),j=1,jmt),k=1,kmaxMYM)
	write (11) ((uvel(i,k,j,taun),j=1,jmt),k=1,kmaxMYM)
	write (11) ((vvel(i,k,j,taun),j=1,jmt),k=1,kmaxMYM)
	write (11) ((temp(i,k,j,taum),j=1,jmt),k=1,kmaxMYM)
	write (11) ((salt(i,k,j,taum),j=1,jmt),k=1,kmaxMYM)
	write (11) ((uvel(i,k,j,taum),j=1,jmt),k=1,kmaxMYM)
	write (11) ((vvel(i,k,j,taum),j=1,jmt),k=1,kmaxMYM)
#endif
	enddo

	close (11)

	return
	stop
	end
