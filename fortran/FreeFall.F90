module FreeFall

  implicit none

  private
  
  public :: Re_iter, Re_bisection, Re_aersett, iprec
  
  integer, parameter :: iprec=4  ! 4 or 8
  real(kind=iprec), parameter  :: eps=.5e-2

  interface Re_iter
     module procedure Re_iter_count, Re_iter_nocount
  end interface Re_iter

  interface Re_bisection
     module procedure Re_bisection_count, Re_bisection_nocount
  end interface Re_bisection

  
contains

  function f_cliftgauvin(Re)
    real(kind=iprec), intent(IN) :: Re
    real(kind=iprec) :: f_cliftgauvin
    
    f_cliftgauvin = 0.2415*Re**0.687+Re/12.*0.42/(1+19018/(Re**1.16))
  end function f_cliftgauvin

  function error(R, Re)
    real(kind=iprec) ::error
    real(kind=iprec), intent(in) :: r, Re

    error = abs(Re-R/(1+f_cliftgauvin(Re)))/Re
  end function error

  function Re_aersett(R)
    real(kind=iprec)  :: Re_aersett
    real(kind=iprec), intent(IN)  :: R
    
    Re_aersett = R * (1 - (1+(R/2.440)**(-0.4335))**(-1.905))
  
  end function Re_aersett

  function Re_bisection_nocount(R)
    real(kind=iprec)  :: Re_bisection_nocount
    real(kind=iprec), intent(IN)  :: R
    real(kind=iprec)  :: Re_lo, Re_hi, err
    
    Re_lo=0.
    Re_hi=R
    Re_bisection_nocount=R
    err=error(R,Re_bisection_nocount)
    do while(err>eps)
       Re_bisection_nocount=0.5*(Re_hi+Re_lo)
       if(R/(1+f_cliftgauvin(Re_bisection_nocount))<Re_bisection_nocount)then
          Re_hi=Re_bisection_nocount
       else
          Re_lo=Re_bisection_nocount
       endif
       err=(Re_hi-Re_lo)/Re_lo
    enddo
    
  end function Re_bisection_nocount
  
  function Re_bisection_count(R, niter)
    real(kind=iprec)  :: Re_bisection_count
    integer, intent(OUT)         :: niter
    real(kind=iprec), intent(IN)  :: R
    real(kind=iprec)  :: Re_lo, Re_hi, err
    
    Re_lo=0.
    Re_hi=R
    Re_bisection_count=R
    err=error(R,Re_bisection_count)
    niter=0
    do while(err>eps)
       niter=niter+1
       Re_bisection_count=0.5*(Re_hi+Re_lo)
       if(R/(1+f_cliftgauvin(Re_bisection_count))<Re_bisection_count)then
          Re_hi=Re_bisection_count
       else
          Re_lo=Re_bisection_count
       endif
       err=(Re_hi-Re_lo)/Re_lo
    enddo
    
  end function Re_bisection_count
  
  function Re_iter_nocount(R)
    real(kind=iprec), intent(IN) :: R
    real(kind=iprec) :: Re_iter_nocount
    real(kind=iprec)             :: d, dprev, err
    
  ! Evaluate the delta function as a function of the virtual Reynolds number Re
        d=0.
        err=error(R,R)
        do while(err>eps)
            dprev = d
            d = 1./(1+f_cliftgauvin(R*(1+d))) - 1.
            err = abs(d-dprev)/(1+d)
        enddo
        Re_iter_nocount=(1+d)*R
  end function Re_iter_nocount
      
  function Re_iter_count(R, niter)
    real(kind=iprec), intent(IN) :: R
    integer, intent(OUT)         :: niter
    real(kind=iprec) :: Re_iter_count
    real(kind=iprec)             :: d, dprev, err
    
  ! Evaluate the delta function as a function of the virtual Reynolds number Re
        d=0.
        err=error(R,R)
        niter=0
        do while(err>eps)
            niter=niter+1
            dprev = d
            d = 1./(1+f_cliftgauvin(R*(1+d))) - 1.
            err = abs(d-dprev)/(1+d)
        enddo
        Re_iter_count=(1+d)*R
  end function Re_iter_count
    
end module FreeFall
