	subroutine stamp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	This subroutine is inherited to snaps.F. This is to make daily
c	time-stamp for writing netcdf output.
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : nn
            use param_mod, only : dpm, stamp_day_daily, stamp_time_daily, tr
        implicit none
   
	character(len=3) :: month_name(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', 'JUL','AUG','SEP','OCT','NOV','DEC'/)

        integer :: iy, mm, im, itemp, nnn, nt
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
	end subroutine stamp 
