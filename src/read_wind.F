	subroutine read_wind (iread_index)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c
!c
!c	attempt to read wind forcing arrays adequate to model time stepts
!c
!c
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod
            use param_mod
        implicit none
        integer :: isteps, iread_index
	open (11,file='tau.swm',form='unformatted',access='sequential')
!c	open (11,file='/gpfs3/home/valsala/FLUXES/TROPFLUX/tau.swm',
!c     &           form='unformatted',access='sequential')
	do isteps = 1, iread_index
!c skip array until or iread_index
	if ( isteps .lt. iread_index) then
	 read (11)
	 read (11)
	else
	read (11) ((taux(i,j,1),i=1,imt),j=1,jmt) ! SHIKHA what is the third index
	read (11) ((tauy(i,j,1),i=1,imt),j=1,jmt) ! is taux(i,j,1) =taux(i,j,2)
	read (11) ((taux(i,j,2),i=1,imt),j=1,jmt) ! or reading the next
!record,in one loop 2 time steps are ebing read and analysed 
	read (11) ((tauy(i,j,2),i=1,imt),j=1,jmt)
	endif
	enddo
	
	close (11)
	wind_switch = 1

	return
	stop
	end
