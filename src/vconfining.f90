! Function that defines an external confining potential that can be deployed in order to reduce
! the spread of weakly-bound (or even initially unbound) states

function vconfining(rho, height)

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Arguments
   real(dp), allocatable :: rho(:) ! The total density of the system
   real(dp) :: height ! The height of the confining potential

   ! Returns
   real(dp), allocatable :: vconfining(:) ! The confining potential

   ! Local variables
   integer :: ii
   integer :: mmax
   real(dp), parameter :: rho_cutoff = 0.00035_dp ! The density below which the confining potential is applied
   real(dp), parameter :: beta = 1.3_dp ! A parameter that affects the smoothness of the confining potential

   mmax = size(rho)

   allocate(vconfining(mmax))

   do ii = 1, mmax
      vconfining(ii) = 0.5_dp * height * (1 + (1 - (rho(ii)/rho_cutoff)**(2*beta))/(1 + (rho(ii)/rho_cutoff)**(2*beta)))
   end do

end function vconfining
