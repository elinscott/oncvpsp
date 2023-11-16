! Function that defines an external confining potential that can be deployed in order to reduce
! the spread of weakly-bound (or even initially unbound) states

! The confining potential is set to zero as r -> infinity

function vconfining(rho, rho0, depth, beta)

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Arguments
   real(dp), allocatable :: rho(:) ! The total density of the system
   real(dp) :: rho0 ! The density below which the confining potential is applied
   real(dp) :: depth ! The depth of the confining potential
   real(dp) :: beta ! A parameter that affects the smoothness of the confining potential

   ! Returns
   real(dp), allocatable :: vconfining(:) ! The confining potential

   ! Local variables
   integer :: ii
   integer :: mmax
   real(dp) :: rho_max

   mmax = size(rho)

   allocate(vconfining(mmax))

   rho_max = rho(mmax)
   do ii = mmax, 1, -1
      ! We define the confining potential as a function of the maximum density for radii greater than the current value of r. This ensures that it is monotonically increasing.
      if (rho(ii) > rho_max) rho_max = rho(ii)
      vconfining(ii) = 0.5_dp * depth * ((1 - (rho_max/rho0)**(2*beta))/(1 + (rho_max/rho0)**(2*beta)) - 1)
   end do



end function vconfining
