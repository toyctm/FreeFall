module thermo
  use FreeFall, only : iprec

  implicit none
  
  private

  public :: rho, mu, mfp, R
  
  real(kind=iprec), parameter :: R=287.1

contains
  
  function rho(P, T)
    real(kind=iprec), intent(IN) :: P, T
    real(kind=iprec) :: rho

    rho = P / (R*T)
  end function rho

  function mu(P,T)
    real(kind=iprec), intent(IN) :: P, T
    real(kind=iprec) :: mu
    real(kind=iprec), parameter :: beta=1.458e-6
    real(kind=iprec), parameter :: S=110.4
    
    mu =  beta*T**1.5/(T+S)
 
  end function mu

  function mfp(P,T)
    real(kind=iprec), intent(IN) :: P, T
    real(kind=iprec) :: mfp
    real(kind=iprec), parameter :: pi = 3.14159265359
    mfp = sqrt(pi/8.)*mu(P,T)/0.4987445*1/sqrt(P*rho(P,T)) ! Jennings 1988
  end function MFP
 
end module thermo
