! Function that defines an external tanh potential that can be deployed in order to reduce
! the spread of weakly-bound (or even initially unbound) states

! The tanh potential is set to zero as r -> infinity

function vtanh(rr, rc, depth, beta)

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Arguments
   real(dp), allocatable :: rr(:) ! The radius at which the tanh potential is applied
   real(dp) :: rc
   real(dp) :: depth ! The depth of the tanh potential
   real(dp) :: beta ! A parameter that affects the smoothness of the tanh potential

   ! Returns
   real(dp), allocatable :: vtanh(:) ! The tanh potential

   ! Local variables
   integer :: ii, mmax

   mmax = size(rr)

   allocate(vtanh(mmax))

   do ii = 1, mmax
      vtanh(ii) = depth / 2.0_dp * (tanh(beta * (rr(ii) - rc)) - 1)
   end do

end function vtanh
