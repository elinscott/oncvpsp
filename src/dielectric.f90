function dielectric(rho, eps_infinity)

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   !Input vaiables
   real(dp), allocatable :: rho(:)
   real(dp) :: eps_infinity

   !Output variable
   real(dp), allocatable :: dielectric(:)

   !Local variables
   integer :: ii
   integer :: mmax
   real(dp), parameter :: rho0 = 0.00035d0 ! atomic units
   real(dp), parameter :: beta = 1.3d0 ! dimensionless

   mmax = size(rho)
   allocate(dielectric(mmax))

   do ii = 1, mmax
      dielectric(ii) = 1.0d0 + (eps_infinity - 1)/2*(1 + (1 - (rho(ii)/rho0)**(2*beta))/(1 + (rho(ii)/rho0)**(2*beta)))
   end do

end function dielectric
