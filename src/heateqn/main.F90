
program main 
  use mpp_mod, only : mpp_init, mpp_npes, mpp_pe, mpp_exit 
  use fms_mod,  only : fms_init, fms_end, write_data
  use fms_io_mod, only : fms_io_init, fms_io_exit

  use mpp_domains_mod, only : domain2d, mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, CGRID_SW, mpp_update_domains

  implicit none
  integer :: nx = 1000, ny = 500, nt = 100
  integer :: i, j, t, domain_layout(2)
  real, allocatable :: u(:,:), u1(:,:)
  real :: dt=1., dx=1., dy=1., k = 0.5, cx, cy
  type(domain2d) :: domain
  integer:: isc, iec, jsc, jec
  integer:: isd, ied, jsd, jed
  integer:: js, je 
  call mpp_init()
  call fms_init()
  call fms_io_init()
 
  call mpp_define_layout((/1,nx,1,ny/),mpp_npes(),domain_layout)

  call mpp_define_domains((/1,nx,1,ny/), domain_layout, domain, xhalo=1, yhalo=1)
  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
  call mpp_get_data_domain(domain, isd, ied, jsd, jed)
   
  cx = dt*k/(dx*dx)
  cy = dt*k/(dy*dy)

  allocate(u(isd:ied,jsd:jed))

  allocate(u1(isc:iec,jsc:jec))

  u(:,:) = 0.

  if (isd==0 .and. jsd>=200 .and. jed<=300) then
    js = max(jsd,200); je = min(jed,300)
    u(0,js:je) = 10. !boundary condition
  endif


  do t = 1, nt

    do i = isc, iec
      do j = jsc, jec

        u1(i,j) = u(i,j) + cx * (u(i-1,j) - 2*u(i,j) + u(i+1,j)) &
                         + cy * (u(i,j-1) - 2*u(i,j) + u(i,j+1))

      end do
    end do

    u(isc:iec,jsc:jec) = u1(isc:iec,jsc:jec)

    call mpp_update_domains(u,domain)
  end do

  call write_data('heateqn_out', 'u', u(isc:iec,jsc:jec), domain=domain)

  call fms_io_exit()
  call fms_end()
  call mpp_exit()

end program main

