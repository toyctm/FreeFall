program testspeed
use FreeFall, only : iprec, Re_iter, Re_bisection, Re_aersett
  
implicit none

real(kind=iprec) :: R, dumsum
integer      :: i, j
integer, parameter :: testsize=10000000
real(kind=iprec), dimension(testsize) :: testvec
real(kind=iprec), parameter :: min_exp=1., max_exp=4.
integer :: count_ini, count, count_max
real(kind=8) :: count_rate
CALL SYSTEM_CLOCK(count, count_rate, count_max)
WRITE(*,*) count, count_rate, count_max

do i=0,100
   R = 10** (min_exp+(i/100.)*(max_exp-min_exp))
   print*,i, R, Re_iter(R), Re_bisection(R), Re_aersett(R), &
        & (Re_aersett(R)-Re_bisection(R))/Re_bisection(R)*100.,(Re_iter(R)-Re_bisection(R))/Re_bisection(R)*100.
enddo

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
print*, 'Aersett execution time (s):', (count-count_ini)/count_rate
print*,'dumsum',dumsum



! test speed of iter
dumsum=0.
CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
do i=1, testsize
   dumsum=dumsum+Re_iter(testvec(i))
enddo
CALL SYSTEM_CLOCK(count, count_rate, count_max)
print*, 'Iter execution time (s):', (count-count_ini)/count_rate
print*,'dumsum',dumsum

! test speed of bisection
dumsum=0.
CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
do i=1, testsize
   dumsum=dumsum+Re_bisection(testvec(i))
enddo
CALL SYSTEM_CLOCK(count, count_rate, count_max)
print*, 'Bisection execution time (s):', (count-count_ini)/count_rate
print*,'dumsum',dumsum

end program testspeed
