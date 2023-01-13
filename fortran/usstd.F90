module usstd
  use FreeFall, only : iprec
  use thermo, only : R
  implicit none

  private

  public :: usstd_profiles

  real(kind=iprec), parameter, dimension(8) :: zbot = [0.      , 11000. , 20000. , 32000. , 47000. , 51000.  , 71000. , 100000.]
  real(kind=iprec), parameter, dimension(8) :: lr = [-0.0065 , 0.0    , 0.001  , 0.0028 , 0.0    , -0.0028 , -0.002 , 0.]

  real(kind=iprec), parameter :: g=9.81
  
  
contains

  subroutine usstd_profiles(nz, zvals, pvals, tvals)
    integer, intent(IN) :: nz
    real(kind=iprec), intent(IN), dimension(nz) :: zvals
    real(kind=iprec), intent(OUT), dimension(nz) :: pvals, tvals

    real(kind=iprec) :: tprev, prev, zprev , lapser, z, dz, pcur, pprev, tcur
    integer :: islice, iz
    
    tprev= 288.15
    pprev= 101325.
    pvals(1)=pprev
    tvals(1) = tprev
    zprev = zvals(1)
    islice=0
    
    do iz=2,nz
       z=zvals(iz)
       if (z >= zbot(islice+1))then
          islice=islice+1
       endif
       dz=z-zprev
       tcur=tprev+lapser*dz
       pcur=pprev-g/(R*0.5*(tcur+tprev))*pprev*dz
       lapser=lr(islice)
       tvals(iz)=tcur
       pvals(iz)=pcur

       zprev=z
       tprev=tcur
       pprev=pcur
    enddo
  end subroutine usstd_profiles
  
end module usstd
