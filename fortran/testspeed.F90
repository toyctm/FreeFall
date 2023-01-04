program testspeed
  use FreeFall, only : iprec, Re_iter, Re_bisection, Re_aersett, vinfty_iter, vinfty_bisection, vinfty_aersett
  use usstd, only : usstd_profiles
  use thermo, only : rho, mu, mfp
  
  implicit none

  real(kind=iprec) :: R, dumsum, dum, avgiter, P, T
  integer, parameter :: nz=12001
  real(kind=iprec), dimension(nz) :: zvals, pvals, tvals, rho_a_vals, lambda_vals, mu_vals
  integer :: iz

  ! Build a vector with evenly spaced altitudes from 0 to 12000m
  do iz=1, nz
     zvals(iz)=real(iz-1, kind=iprec)
  enddo

  ! Calculate 

  call usstd_profiles(nz, zvals, pvals, tvals)

  do iz=1, nz
     P = pvals(iz)
     T = tvals(iz)
     rho_a_vals(iz)  = rho(P, T)
     lambda_vals(iz) = mfp(P, T)
     mu_vals(iz)     = mu(P,T)
  enddo
  
  print*,'==============================================='
  print*,'SAMPLE VALUES'
  call testvalues(-7., -3.)
  print*,'==============================================='
  print*,'RANGE -7. TO -6.'
  call calcspeed(-7., -6., 100000000, nz, lambda_vals, rho_a_vals, mu_vals)
  print*,'==============================================='
  print*,'RANGE -6. TO -5.'
  call calcspeed(-6., -5. , 100000000, nz, lambda_vals, rho_a_vals, mu_vals)
  print*,'==============================================='
  print*,'RANGE -5. TO -4.'
  call calcspeed(-5., -4. , 100000000, nz, lambda_vals, rho_a_vals, mu_vals)
  print*,'==============================================='
  print*,'RANGE -4. TO -3.'
  call calcspeed(-4., -3. , 100000000, nz, lambda_vals, rho_a_vals, mu_vals)
  print*,'==============================================='

end program testspeed

subroutine testvalues(min_exp, max_exp)
  use FreeFall, only : iprec, vinfty_iter, vinfty_bisection, vinfty_aersett
  use thermo, only : mfp, rho, mu
  
  implicit none
  real(kind=iprec), intent(IN)  :: min_exp, max_exp
  real(kind=iprec) :: D, P, T, lambda_a, mu_a, rho_a, rho_p
  integer          :: i
  P = 101350.
  T=290.
  mu_a = mu(P,T)
  lambda_a = mfp(P,T)
  rho_a = rho(P,T)
  rho_p = 2650.
  print*,mu_a, lambda_a, rho_a
  do i=0,100
   D = 10** (min_exp+(i/100.)*(max_exp-min_exp))
   print*,i, D, vinfty_iter(D, rho_p, lambda_a, rho_a, mu_a), &
        & vinfty_bisection(D, rho_p, lambda_a, rho_a, mu_a), &
        & vinfty_aersett(D, rho_p, lambda_a, rho_a, mu_a)
   enddo

end subroutine testvalues

subroutine calcspeed(min_exp, max_exp, testsize, nz, lambda_vals, rho_a_vals, mu_vals)
  use FreeFall, only : iprec, vinfty_iter, vinfty_bisection, vinfty_aersett

  implicit none
  
  integer, intent(IN) :: testsize, nz
  real(kind=iprec), intent(IN)  :: min_exp, max_exp
  real(kind=iprec), intent(IN), dimension(nz)  :: lambda_vals, rho_a_vals, mu_vals
  real(kind=iprec), dimension(testsize) :: testvec, lambdavec, rho_a_vec, muvec
  integer, dimension(testsize) :: izvec
  real(kind=iprec)    :: dumsum, dum, avgiter, xrand
  real(kind=iprec), parameter    :: rho_p = 2650.
  
  integer      :: i, j, niter
  integer(kind=8) :: nitersum
  integer :: count_ini, count, count_max
  real(kind=8) :: count_rate
  
  ! Generate a large random vector of R between 1 and 1.e5

  ! Generate a list of diameters between 10**min_exp and 10**max_exp
  call random_number(testvec)
  testvec=min_exp+(max_exp-min_exp)*testvec
  testvec=10**testvec

  ! Generate a list of random altitudes
  do i=1, testsize
     call random_number(xrand)
     if(xrand==1.)call random_number(xrand)
     izvec(i) = 1+floor(xrand*nz)
  enddo
  do i=1, testsize
     lambdavec(i) = lambda_vals(izvec(i))
     muvec(i) = mu_vals(izvec(i))
     rho_a_vec(i) = rho_a_vals(izvec(i))
!     print*,izvec(i), lambdavec(i), muvec(i), rho_a_vec(i)
  enddo
  
  dumsum=0.
  CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
  ! test speed of AerSett
  do i=1, testsize
     dumsum=dumsum+vinfty_aersett(testvec(i), rho_p, lambdavec(i), rho_a_vec(i), muvec(i))
  enddo
  CALL SYSTEM_CLOCK(count, count_rate, count_max)
  print*, 'Aersett execution time (ns/call):', (count-count_ini)/count_rate*1.e+9/testsize
  print*,'dumsum',dumsum


  ! Count iterations of iter
  dumsum=0.
  nitersum=0
  do i=1, testsize
     dum=vinfty_iter(testvec(i), rho_p, lambdavec(i), rho_a_vec(i), muvec(i), niter)
     nitersum=nitersum+niter
  enddo
  avgiter=(1.*nitersum)/(1.*testsize)

  ! test speed of iter
  dumsum=0.
  CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
  do i=1, testsize
     dumsum=dumsum+vinfty_iter(testvec(i), rho_p, lambdavec(i), rho_a_vec(i), muvec(i))
  enddo
  CALL SYSTEM_CLOCK(count, count_rate, count_max)
  print*, 'Iter execution time (ns/call):', (count-count_ini)/count_rate*1.e+9/testsize
  print*,'dumsum',dumsum
  print*,'Number of iterations', avgiter

  ! Count iterations of bisection
  dumsum=0.
  nitersum=0
  do i=1, testsize
     dum=vinfty_bisection(testvec(i), rho_p, lambdavec(i), rho_a_vec(i), muvec(i), niter)
     nitersum=nitersum+niter
  enddo
  avgiter=(1.*nitersum)/(1.*testsize)

  ! test speed of bisection
  dumsum=0.
  CALL SYSTEM_CLOCK(count_ini, count_rate, count_max)
  do i=1, testsize
     dumsum=dumsum+vinfty_bisection(testvec(i), rho_p, lambdavec(i), rho_a_vec(i), muvec(i))
  enddo
  CALL SYSTEM_CLOCK(count, count_rate, count_max)
  print*, 'Bisection execution time (ns/call):', (count-count_ini)/count_rate*1.e+9/testsize
  print*,'dumsum',dumsum
  print*,'Number of iterations', avgiter

end subroutine calcspeed
