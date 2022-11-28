module FreeFall

  implicit none

  private
  
  public :: delta_iter, iprec
  
  integer, parameter :: iprec=4  ! 4 or 8


contains

  function f_cliftgauvin(Re)
    real(kind=iprec), intent(IN) :: Re
    real(kind=iprec) :: f_cliftgauvin
    
    f_cliftgauvin = 0.2415*Re**0.687+Re/12.*0.42/(1+19018/(Re**1.16))
  end function f_cliftgauvin
  
  function delta_iter(R)
    real(kind=iprec), intent(IN) :: R
    real(kind=iprec) :: delta_iter
    real(kind=iprec), parameter  :: eps=5.e-3
    real(kind=iprec)             :: d, dprev, error
    
  ! Evaluate the delta function as a function of the virtual Reynolds number Re
        d=0.
        error=1.
        do while(error>eps)
            dprev = d
            d = 1./(1+f_cliftgauvin(R*(1+d))) - 1.
            error = abs(d-dprev)
         enddo
        delta_iter=d
  end function delta_iter
    
  subroutine print_hello()

    print*,'Hello World'

  end subroutine print_hello
  
end module FreeFall
