program datefunc
  use time_manager_mod
  use stringfunc, only: to_lower

  implicit none

  type incdate_type
    character (len=256) :: calendar
    character (len=14) :: date
    integer :: inc(6)
  endtype
  
  type (incdate_type) :: incdateinp

  type(time_type) :: Time, Time_end
  character(len=256) :: carg, exename
  integer :: date(6), i, m, iostat
  double precision :: rtime, half=0.5


  call getarg(0,exename)
  call getarg(1,carg)
  carg=to_lower(carg)
  
  if (trim(carg)=='incdate') then
    call getarg(2,carg)
    carg=to_lower(carg)
    read (carg,*,iostat=iostat) incdateinp
    if (iostat/=0) then
      write(0,*) " "
      write(0,*) 'Error while reading input arguments'
      write(0,*) " "
      write(0,*) 'Usage:'
      write(0,*) trim(exename)//' incdate "calendartype date as (yyyymmddhhmnsc) &
            0  & increment as (yyyy,mm,dd,hh,mns,secs)" '
      write(0,*) " "
      write(0,*) 'Example:'
      write(0,*) trim(exename)//' incdate "julian 18000101000000 0,0,1,0,0,0" ' 
      stop 1
    endif
    
    call set_caltype(incdateinp%calendar)

    read(incdateinp%date,'(I4.4,5(I2.2))',iostat=iostat) date 
    if (iostat/=0) then
      write(0,*) 'Error while reading input date.'
      stop 1
    endif

    Time = set_date (date(1), date(2), date(3),  &
       date(4), date(5), date(6))
    Time_end = Time

    do m=1,incdateinp%inc(1)
      Time_end = Time_end + set_time(0,days_in_year(Time_end))
    end do

    do m=1,incdateinp%inc(2)
      Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do

    Time_end   = Time_end + set_time(incdateinp%inc(4)*3600 + &
                                     incdateinp%inc(5)*60 + &
                                     incdateinp%inc(6), &
                                     incdateinp%inc(3))
 
    call get_date(Time_end,date(1),date(2),date(3),date(4),date(5),date(6))

    write(*,'(I4.4,5(I2.2))') date(:)
    stop
  endif

  if (trim(carg)=='incdatemid')  then 
    call getarg(2,carg)
    carg=to_lower(carg)
    read (carg,*,iostat=iostat) incdateinp
    if (iostat/=0) then
      write(0,*) " "
      write(0,*) 'Error while reading input arguments'
      write(0,*) " "
      write(0,*) 'Usage:'
      write(0,*) trim(exename)//' incdatemid "calendartype date as (yyyymmddhhmnssecs) &
            0  & increment as (yyyy,mm,dd,hh,mns,secs)" '
      write(0,*) " "
      write(0,*) 'Example:'
      write(0,*) trim(exename)//' incdatemid "julian 18000101000000 0,0,1,0,0,0" ' 
      stop 1
    endif
    
    call set_caltype(incdateinp%calendar)
  
    read(incdateinp%date,'(I4.4,5(I2.2))',iostat=iostat) date 
    if (iostat/=0) then
      write(0,*) 'Error while reading input date.'
      stop 1
    endif

    Time = set_date (date(1), date(2), date(3),  &
       date(4), date(5), date(6))
    Time_end = Time

    do m=1,incdateinp%inc(1)
      Time_end = Time_end + set_time(0,days_in_year(Time_end))
    end do

    do m=1,incdateinp%inc(2)
      Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do

    Time_end   = Time_end + set_time(incdateinp%inc(4)*3600 + &
                                     incdateinp%inc(5)*60 + &
                                     incdateinp%inc(6), &
                                     incdateinp%inc(3))
    
    Time_end = real_to_time_type(real(time_type_to_real(Time_end-Time)*0.5))

    Time_end = Time + Time_end
 
    call get_date(Time_end,date(1),date(2),date(3),date(4),date(5),date(6))

    write(*,'(I4.4,5(I2.2))') date(:)
    stop
  endif

  write(0,*) 'Wrong or no operators given'
  write(0,*) 'Valid operators are:'
  write(0,*) 'incdate incdatemid'  
  contains
  
  subroutine set_caltype(caltype)
    character (len=*), intent (in) :: caltype
    
    select case(trim(caltype))
      case('thirty_day_months')
        call set_calendar_type(thirty_day_months)
      case('julian')
        call set_calendar_type(julian)
      case('noleap')
        call set_calendar_type(noleap)
      case default
        write(0,*) "Error: unrecognised calendar type: "//trim(caltype)
        write(0,*) "Available calendar types are:"
        write(0,*) "thirty_day_months, julian, noleap"
        stop 1
    end select
  end subroutine set_caltype

end program datefunc
