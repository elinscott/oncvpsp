subroutine screen_potential_via_gaussian_charge(mmax, rr, zz, eps, v)

   implicit none
   integer, parameter :: dp = kind(1.0d0)
   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp

   integer, intent(in) :: mmax
   real(dp), intent(in) :: rr(mmax), eps(mmax)
   real(dp), intent(in) :: zz
   real(dp), intent(inout) :: v(mmax)

   integer :: ii
   real(dp), allocatable :: rho_gaussian(:), eps_vacuum(:), v_gaussian_vacuum(:), v_gaussian_eps(:)
   real(dp) :: al, amesh
   real(dp), parameter :: sigma = 0.1d0 ! The smearing of the Gaussian

   allocate(rho_gaussian(mmax), eps_vacuum(mmax), v_gaussian_eps(mmax), v_gaussian_vacuum(mmax))
   eps_vacuum(:) = 1.0d0
   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   ! The potential due to a Gaussian charge distribution in a vacuum has an analytical solution
   do ii = 1, mmax
      v_gaussian_vacuum(ii) =  -1.0d0 * zz / rr(ii) * erf(rr(ii) / sqrt(2.0d0) / sigma)
   end do

   ! The potential due to a Gaussian charge distribution in a non-homogenous dielectric does not have an analytical solution
   ! We solve it numerically by solving the Poisson equation

   do ii = 1, mmax
      rho_gaussian(ii) = -1.0_dp * zz / (sqrt(pi / 2.0_dp) * sigma**3) * dexp(-rr(ii)**2 / (2.0d0*sigma**2))
   end do
   call poisson_solver(rr, mmax, rho_gaussian, -1.0d0 * zz, eps, v_gaussian_eps)

   ! Update the potential
   v(:) = v(:) + v_gaussian_eps(:) - v_gaussian_vacuum(:)

   ! Printing for debugging
   ! do ii = 1, mmax
   !    write(80, '(7f25.8)') rr(ii), rho_gaussian(ii), v(ii), -1.0d0 * zz/rr(ii), v_gaussian_vacuum(ii), v_gaussian_eps(ii), 
   ! end do
   ! write(80, *) ''
   ! stop

   deallocate(rho_gaussian, eps_vacuum, v_gaussian_vacuum, v_gaussian_eps)
end subroutine screen_potential_via_gaussian_charge