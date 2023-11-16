subroutine lschvkbb_integrate_inward(ll, nin, mmax, mch, rr, vloc, cf, ee, uu, up, inflexion)

   implicit none
   integer, parameter :: dp = kind(1.0d0)

   ! Argument variables
   integer, intent(in) :: ll, nin, mmax
   integer, intent(inout) :: mch
   real(dp), intent(in) :: rr(mmax), vloc(mmax), cf(mmax), ee
   real(dp), intent(inout) :: uu(mmax), up(mmax)
   logical, intent(out) :: inflexion

   ! Local variables
   integer :: ii, it, jj
   real(dp) :: al, sls, xkap, max_uu
   real(dp), allocatable :: upp(:)

   ! Functions
   interface
      function aei(yy, ii)
         implicit none
         integer, parameter :: dp = kind(1.0d0)
         integer :: ii
         real(dp) :: yy(*)
         real(dp) :: aei
      end function aei
   end interface

   interface
      function aii(yy, ii)
         implicit none
         integer, parameter :: dp = kind(1.0d0)
         integer :: ii
         real(dp) :: yy(*)
         real(dp) :: aii
      end function aii
   end interface

   al = 0.01d0*dlog(rr(101)/rr(1))
   sls = ll*(ll + 1)

   allocate(upp(mmax))
   upp(:) = 0.0d0

   inflexion = .false.

   ! start inward integration at 10*classical turning
   ! point with simple exponential
   xkap = dsqrt(sls/rr(nin)**2 + 2.0d0*(vloc(nin) - ee))
   
   do ii = nin, nin + 4
      uu(ii) = exp(-xkap*(rr(ii) - rr(nin)))
      up(ii) = -rr(ii)*al*xkap*uu(ii)
      upp(ii) = al*up(ii) + cf(ii)*uu(ii)
   end do
         
   ! integrate inward
   do ii = nin, mch + 1, -1
      uu(ii - 1) = uu(ii) + aei(up, ii)
      up(ii - 1) = up(ii) + aei(upp, ii)
      do it = 1, 2
         upp(ii - 1) = al*up(ii - 1) + cf(ii - 1)*uu(ii - 1)
         up(ii - 1) = up(ii) + aii(upp, ii)
         uu(ii - 1) = uu(ii) + aii(up, ii)
      end do

      ! Keep the exponential growth under control
      max_uu = maxval(uu(ii - 1:nin))
      if (max_uu > 1.0d10) then
         do jj = nin, ii - 1, -1
            uu(jj) = uu(jj) / max_uu
            up(jj) = up(jj) / max_uu
            upp(jj) = upp(jj) / max_uu
         end do
      end if

      ! if (upp(ii) < 0.0d0) then
      !    mch = min(ii + 10, mmax - 10)
      !    inflexion = .true.
      !    exit
      ! end if
   end do
   open(100, file='lschvkbb_integrate_inward.dat', status='unknown')
   do ii = nin, mch, -1
      write(100,'(4e16.8)') rr(ii), uu(ii), up(ii), upp(ii)
   end do
   close(100)

   deallocate(upp)

end subroutine lschvkbb_integrate_inward
