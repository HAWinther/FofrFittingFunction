program main
  use FittingFormula, only: pofk_enhancement
  implicit none
  integer, parameter :: dl = KIND(1.d0)
  real(dl) :: ratio, z, fofr0, k
  integer, parameter :: nk = 30
  real(dl), parameter :: kmin = 0.01, kmax = 10.0
  integer :: i, argc
  character(len=128) :: argv

  argc = iargc()
  if(argc < 2)then
    write(*,*) 'Run as ./test fofr0 redshift'
    stop
  endif
  
  call getarg(1,argv)
  read(argv,*) fofr0
  call getarg(2,argv)
  read(argv,*) z

  write(*,*) "For fR0 = ", fofr0, "and z = ", z, " we have the following (k, ratio(k)):"
  do i=1,nk
    k = exp(log(kmin) + log(kmax/kmin)*(i-1)/dble(nk-1))
    ratio = pofk_enhancement(z,fofr0,k)
    write(*,*) k, ratio
  enddo
end program main
