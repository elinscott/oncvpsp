! Function that defines an external tanh potential that can be deployed in order to reduce
! the spread of weakly-bound (or even initially unbound) states

! The tanh potential is set to zero as r -> infinity

function vtanhlog(rr, rc, depth, beta)

   implicit none

   integer, parameter :: dp = kind(1.0d0)

   ! Arguments
   real(dp), allocatable :: rr(:) ! The radius at which the tanh potential is applied
   real(dp) :: rc
   real(dp) :: depth ! The depth of the tanh potential
   real(dp) :: beta ! A parameter that affects the smoothness of the tanh potential

   ! Returns
   real(dp), allocatable :: vtanhlog(:) ! The tanh potential

   ! Local variables
   integer :: ii, mmax

   mmax = size(rr)

   allocate(vtanhlog(mmax))

   do ii = 1, mmax
      vtanhlog(ii) = -1.0_dp * depth / ( 1.0_dp + (rr(ii)/rc) ** (2.0_dp * beta))
   end do

end function vtanhlog
