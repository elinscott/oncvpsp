function spherically_integrate(mmax, rr, obj)
   ! Integrates a radial function on a logarithmic grid

   implicit none

   ! Input
   integer, parameter :: dp = kind(1.0d0)
   integer :: mmax
   real(dp) :: rr(mmax), obj(mmax)

   ! Output
   real(dp) :: spherically_integrate

   ! Local variables
   integer :: ii
   real(dp) :: al
   real(dp), allocatable :: objr3(:)

   allocate(objr3(mmax))
   objr3(:) = 0.0d0

   al = 0.01d0*dlog(rr(101)/rr(1))

   do ii = 1, mmax
      objr3(ii) = obj(ii)*rr(ii)**3
   end do

   spherically_integrate = (9.0d0*objr3(1) + 28.0d0*objr3(2) + 23.0d0*objr3(3))/24.0d0
   do ii = 4, mmax
      spherically_integrate = spherically_integrate + objr3(ii)
   end do
   spherically_integrate = al*spherically_integrate + objr3(1)/3.0d0

   deallocate(objr3)

end function
