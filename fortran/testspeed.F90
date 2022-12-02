program testspeed
  use FreeFall, only : iprec, Re_iter, Re_bisection, Re_aersett
  
  implicit none

  real(kind=iprec) :: R, dumsum, dum, avgiter

  print*,'==============================================='
  print*,'SAMPLE VALUES'
  call testvalues(0., 4.)
  print*,'==============================================='
  print*,'RANGE 0. TO 1.'
  call calcspeed(0., 1., 10000000)
  print*,'==============================================='
  print*,'RANGE 1. TO 2.'
  call calcspeed(1., 2., 10000000)
  print*,'==============================================='
  print*,'RANGE 2. TO 3.'
  call calcspeed(2., 3., 10000000)
  print*,'==============================================='
  print*,'RANGE 3. TO 4.'
  call calcspeed(3., 4., 10000000)
  print*,'==============================================='

end program testspeed

subroutine testvalues(min_exp, max_exp)
  use FreeFall, only : iprec, Re_iter, Re_bisection, Re_aersett

  implicit none
  real(kind=iprec), intent(IN)  :: min_exp, max_exp
  real(kind=iprec) :: R
  integer          :: i
  do i=0,100
   R = 10** (min_exp+(i/100.)*(max_exp-min_exp))
   print*,i, R, Re_iter(R), Re_bisection(R), Re_aersett(R), &
        & (Re_aersett(R)-Re_bisection(R))/Re_bisection(R)*100.,(Re_iter(R)-Re_bisection(R))/Re_bisection(R)*100.
   enddo

end subroutine testvalues

subroutine calcspeed(min_exp, max_exp, testsize)
  use FreeFall, only : iprec, Re_iter, Re_bisection, Re_aersett

  implicit none
  
  integer, intent(in) :: testsize
  real(kind=iprec), intent(IN)  :: min_exp, max_exp
  real(kind=iprec), dimension(testsize) :: testvec
  real(kind=iprec)    :: dumsum, dum, avgiter
  integer      :: i, j, niter
  integer(kind=8) :: nitersum
  integer :: count_ini, count, count_max
  real(kind=8) :: count_rate

  ! Generate a large random vector of R between 1 and 1.e5

  call random_number(testvec)
  testvec=min_exp+(max_exp-min_exp)*testvec
  testvec=10**testvec
  
  dumsum=0.
  CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
  ! test speed of AerSett
  do i=1, testsize
     dumsum=dumsum+Re_aersett(testvec(i))
  enddo
  CALL SYSTEM_CLOCK(count, count_rate, count_max)
  print*, 'Aersett execution time (ns/call):', (count-count_ini)/count_rate*1.e+9/testsize
  print*,'dumsum',dumsum


  ! Count iterations of iter
  dumsum=0.
  nitersum=0
  do i=1, testsize
     dum=Re_iter(testvec(i), niter)
     nitersum=nitersum+niter
  enddo
  avgiter=(1.*nitersum)/(1.*testsize)

  ! test speed of iter
  dumsum=0.
  CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
  do i=1, testsize
     dumsum=dumsum+Re_iter(testvec(i))
  enddo
  CALL SYSTEM_CLOCK(count, count_rate, count_max)
  print*, 'Iter execution time (ns/call):', (count-count_ini)/count_rate*1.e+9/testsize
  print*,'dumsum',dumsum
  print*,'Number of iterations', avgiter

  ! Count iterations of bisection
  dumsum=0.
  nitersum=0
  do i=1, testsize
     dum=Re_bisection(testvec(i), niter)
     nitersum=nitersum+niter
  enddo
  avgiter=(1.*nitersum)/(1.*testsize)

  ! test speed of bisection
  dumsum=0.
  CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
  do i=1, testsize
     dumsum=dumsum+Re_bisection(testvec(i))
  enddo
  CALL SYSTEM_CLOCK(count, count_rate, count_max)
  print*, 'Bisection execution time (ns/call):', (count-count_ini)/count_rate*1.e+9/testsize
  print*,'dumsum',dumsum
  print*,'Number of iterations', avgiter

end subroutine calcspeed
