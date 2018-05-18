	subroutine read_forcing (iread_index)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c	attempt to read wind forcing arrays adequate to model time stepts
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        use size_mod
        use param_mod
!#include "size.h"
!#include "polar.h"
!#include "param.h"

	open (11,file='hum.swm',form='unformatted',access='sequential')
	open (12,file='airt.swm',form='unformatted',access='sequential')
	open (13,file='uwnd.swm',form='unformatted',access='sequential')
	open (14,file='vwnd.swm',form='unformatted',access='sequential')
	open (15,file='ssw.swm',form='unformatted',access='sequential')
	open (16,file='cld.swm',form='unformatted',access='sequential')
	open (18,file='pme.swm',form='unformatted',access='sequential')
	open (19,file='chl.swm',form='unformatted',access='sequential')
	open (20,file='river.swm',form='unformatted',access='sequential')
	do isteps = 1, iread_index
c skip array until or iread_index
	if ( isteps .lt. iread_index) then
	 read (11)
	 read (12)
	 read (13)
	 read (14)
	 read (15)
	 read (16)
	 read (18)
	 read (19)
	 read (20)
	else
	read (11) ((sphm_read(i,j,1),i=1,imt),j=1,jmt)
	read (11) ((sphm_read(i,j,2),i=1,imt),j=1,jmt)
	read (12) ((airt_read(i,j,1),i=1,imt),j=1,jmt)
	read (12) ((airt_read(i,j,2),i=1,imt),j=1,jmt)
	read (13) ((uwnd_read(i,j,1),i=1,imt),j=1,jmt)
	read (13) ((uwnd_read(i,j,2),i=1,imt),j=1,jmt)
	read (14) ((vwnd_read(i,j,1),i=1,imt),j=1,jmt)
	read (14) ((vwnd_read(i,j,2),i=1,imt),j=1,jmt)
	read (15) ((ssw_read(i,j,1),i=1,imt),j=1,jmt)
	read (15) ((ssw_read(i,j,2),i=1,imt),j=1,jmt)
	read (16) ((cld_read(i,j,1),i=1,imt),j=1,jmt)
	read (16) ((cld_read(i,j,2),i=1,imt),j=1,jmt)
	read (18) ((pme_read(i,j,1),i=1,imt),j=1,jmt)
	read (18) ((pme_read(i,j,2),i=1,imt),j=1,jmt)
	read (19) ((chl_read(i,j,1),i=1,imt),j=1,jmt)
	read (19) ((chl_read(i,j,2),i=1,imt),j=1,jmt)
	read (20) ((rvr_read(i,j,1),i=1,imt),j=1,jmt)
	read (20) ((rvr_read(i,j,2),i=1,imt),j=1,jmt)
	endif
	enddo

	close (11)
	close (12)
	close (13)
	close (14)
	close (15)
	close (16)
	close (18)
	close (19)
	close (20)
	wind_switch = 1

	return
	stop
	end
