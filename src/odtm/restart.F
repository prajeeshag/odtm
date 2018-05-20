	subroutine restart
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c	a subroutine to write restart file for the consicutive runs.
c	format is binary.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            use size_mod, only : imt, jmt, km, i, j, k, u, v, h, t, temp, salt
            use size_mod, only : uvel, vvel, taun, taum, kmaxMYM, nn 
        implicit none
            integer :: nt

	open (11, file='restart.swm',form='unformatted',
     &   access='sequential')

	
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
	end subroutine 
