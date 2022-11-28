program testspeed
use FreeFall, only : iprec, delta_iter
  
implicit none

real(kind=iprec) :: R
integer      :: i

do i=1,1000
   R=0.1*i
   print*,i, delta_iter(R)
enddo

end program testspeed
