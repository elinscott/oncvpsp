subroutine poisson_solver(rr, mmax, rho, zion, eps, vo)

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   !Input vaiables
   integer, intent(in)  :: mmax
   real(dp), intent(in) :: rho(mmax), rr(mmax)
   real(dp), intent(in) :: zion
   real(dp), intent(in) :: eps(mmax)

   !Output variable
   real(dp), intent(out) :: vo(mmax)

   !Local variables
   integer ii
   real(dp) :: al, tv
   real(dp), allocatable :: rvp(:), rv(:)
   real(dp), parameter :: rho0 = 0.00035d0 ! atomic units
   real(dp), parameter :: beta = 1.3d0 ! dimensionless

   !Function
   real(dp) :: aii

   allocate (rvp(mmax), rv(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))

! integration for electrostatic potential
   do ii = 1, mmax
      rvp(ii) = rho(ii)*al*rr(ii)**3
   end do

   rv(mmax) = zion
   rv(mmax - 1) = zion
   rv(mmax - 2) = zion

   do ii = mmax - 2, 2, -1
      rv(ii - 1) = rv(ii) + aii(rvp, ii)
   end do

   ! Divide by the dielectric permittivity
   do ii = 1, mmax
      rv(ii) = rv(ii)/eps(ii)
   end do

   do ii = 1, mmax
      rvp(ii) = rho(ii)*al*rr(ii)**2
   end do

   tv = 0.0d0
   do ii = mmax - 2, 2, -1
      tv = tv + aii(rvp, ii)
      rv(ii - 1) = rv(ii - 1) - rr(ii - 1)*tv
   end do

   do ii = 1, mmax
      vo(ii) = rv(ii)/rr(ii)
   end do

   deallocate (rvp, rv)

end subroutine poisson_solver
