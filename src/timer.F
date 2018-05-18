	subroutine timer
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	an attempt to display time ellapsed in the model run
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

	call system_clock(itimer3,itimerrate,itimermax)
	rtimer = real(itimer3-itimer2)/real(itimerrate)
	total_time = rtimer*number_of_snap

	write(*,10) rtimer/60.0
	if ( timer_switch .eq. 0) rtimer1 = rtimer

 10	format (2x,f8.3,$)

c	call system_clock(itimer2,itimerrate,itimermax)
	time_switch = 0


	return
	stop
	end
