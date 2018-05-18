	subroutine stamp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	This subroutine is inherited to snaps.F. This is to make daily
c	time-stamp for writing netcdf output.
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod 
            use param_mod
!#include "size.h"
	character*3 month_name
	dimension month_name(12)
!#include "param.h"
	data month_name(1) /'JAN'/
	data month_name(2) /'FEB'/
	data month_name(3) /'MAR'/
	data month_name(4) /'APR'/
	data month_name(5) /'MAY'/
	data month_name(6) /'JUN'/
	data month_name(7) /'JUL'/
	data month_name(8) /'AUG'/
	data month_name(9) /'SEP'/
	data month_name(10) /'OCT'/
	data month_name(11) /'NOV'/
	data month_name(12) /'DEC'/

	iy = 0
	do mm = 1,12
	im = 0
        itemp = dpm(mm)
	do nnn = 1,itemp
	 iy = iy + 1
	 im = im + 1
	 write(stamp_day_daily(iy),'(a,i2.2,a,i2.2)')
     &     'days since 1992-',mm,'-',im
	 write(stamp_time_daily(iy),'(i2.2,a,a,a)')
     &     im,'-',month_name(mm),'-0000 00:00:00'
	 if (iy .eq. 365) then
	 write(stamp_day_daily(iy),'(a,i2.2,a,i2.2)')
     &     'days since 1992-',mm,'-',im
	 write(stamp_time_daily(iy),'(i2.2,a,a,a)')
     &     im,'-',month_name(mm),'-0000 00:00:00'
	 endif
c	 write(*,*) stamp_day_daily(iy),'   ', stamp_time_daily(iy)
	enddo
	enddo
	
	
c
c for tracer stamps
c
	do nt=1,nn
	 write(Tr(nt),'(a,i2.2)') 'Tr',nt
	enddo
	return
	stop
	end
