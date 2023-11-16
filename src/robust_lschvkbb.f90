!
! Copyright (c) 1989-2019 by D. R. Hamann, Mat-Sim Research LLC and Rutgers
! University
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
subroutine robust_lschvkbb(nn, ll, nvkb, ierr, ee, &
                           & rr, vloc, vkb, evkb, uu, up, mmax, mch, istep, nextra_nodes)

! Finds bound states of a pseudopotential with
! Vanderbilt-Kleinman-Bylander non-local projectors

!nn  principal quantum number
!ll  angular-momentum quantum number
!nvkb  = number of VKB projectors to be used
!ierr  non-zero return if error
!ee  bound-state energy, input guess and output calculated value
!emin  externally generaated estimate of lower bound for ee
!emax  externally generaated estimate of upper bound for ee
!rr  log radial mesh
!vloc  local part of psp
!vkb  VKB projectors
!evkb coefficients of BKB projectors
!uu  output radial wave function (*rr)
!up  d(uu)/dr
!mmax  size of log grid
!mch matching mesh point for inward-outward integrations

   implicit none
   integer, parameter :: dp = kind(1.0d0)

!Input Variables
   real(dp) :: rr(mmax), vloc(mmax), vkb(mmax, nvkb), evkb(nvkb)
   integer :: nn, ll, nvkb, mmax

!Output variables
   real(dp) :: uu(mmax), up(mmax)
   real(dp) :: ee  !in/out, really - needs starting guess
   integer :: ierr, mch
   integer :: istep
   integer :: nextra_nodes

!Local variables
   real(dp) :: cn
   real(dp) :: de
   real(dp) :: eps, ro, sc
   real(dp) :: sls, sn, uout, upin, upout
   real(dp) :: amesh, al, als
   real(dp) :: rc
   integer :: ii, nin, nint, node
   real(dp) :: e_lower, e_upper
   logical :: e_lower_found, inflexion

   real(dp), allocatable :: cf(:)

   interface
      subroutine lschvkbb_integrate_inward(ll, nin, mmax, mch, rr, vloc, cf, ee, uu, up, inflexion)
         implicit none
         integer, parameter :: dp = kind(1.0d0)
         integer, intent(in) :: ll, nin, mmax
         integer, intent(inout) :: mch
         real(dp), intent(in) :: rr(mmax), vloc(mmax), cf(mmax), ee
         real(dp), intent(inout) :: uu(mmax), up(mmax)
         logical, intent(out) :: inflexion
      end subroutine lschvkbb_integrate_inward
   end interface

   allocate (cf(mmax))

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)
   node = 0

! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
!eps=1.0d-10
   eps = 1.0d-8
   ierr = 100

   sls = ll*(ll + 1)

! null arrays to remove leftover garbage
   uu(:) = 0.0d0
   up(:) = 0.0d0

   als = al**2

   e_lower_found = .false.
   e_upper = 0.0_dp
   e_lower = -1.0_dp ! dummy value, will be overwritten

   write(*, *) nn, ll
   ! if (nn == 2 .and. ll == 0) write(*,*) 'GRID SEARCH'

! return point for bound state convergence
   do nint = 1, 1000

      ! if (nn == 2 .and. ll == 0) then
      !    ee = -15.0d0 + 10.0d0 * (float(nint) - 1) / 1000.0d0
      ! end if

! coefficient array for u in differential eq.
      do ii = 1, mmax
         cf(ii) = als*sls + 2.0d0*als*(vloc(ii) - ee)*rr(ii)**2
      end do

! find classical turning point for matching
      mch = 1
      do ii = mmax, 2, -1
         if (cf(ii - 1) <= 0.d0 .and. cf(ii) > 0.d0) then
            mch = ii
            exit
         end if
      end do

      if (mch == 1 .and. nvkb == 0) then
         write (6, '(/a)') 'lschvkbb: ERROR no classical turning point'
         write (6, '(a,i2,a,i6)') 'lschvkbb: ERROR nvkb=', nvkb, 'mch=', mch
         write (6, '(a,i2,a,f8.4,a)') 'lschvkbb: ERROR l=', ll, '  e=', ee
!   stop
      end if

      if (nvkb > 0) then
! find cutoff radius for projectors
         rc = 0.0d0
         do ii = mmax, 1, -1
            if (abs(vkb(ii, 1)) > 0.0d0) then
               rc = rr(ii)
               exit
            end if
         end do

         if (mch == 1) rc = 1.25d0*rc

! adjust matching radius if necessary
         if (rr(mch) < rc) then
            do ii = mch, mmax
               if (rr(ii) > rc) then
                  mch = ii
                  exit
               end if
            end do
         end if
      end if

      ! ! Manually set mch
      ! if (nn == 2 .and. ll == 0) then
      !    if (rr(mch) < 6.0_dp) then
      !       do ii = mch, mmax
      !          if (rr(ii) > 6.0_dp) then
      !             mch = ii
      !             exit
      !          end if
      !       end do
      !    end if
      ! end if

      ! Inward integration to make sure that we won't get any nodes for r > r(mch)
      ! call lschvkbb_integrate_inward(ll, mmax - 4, mmax, mch, rr, vloc, cf, ee, uu, up, inflexion)
      ! if (inflexion) then
      !    write (*, *) 'inflexion found; updating rr(mch) to ', rr(mch)
      ! end if

      write (*, '(i4,3f18.12)', advance='no') nint, ee, e_lower, e_upper

      ! Make sure mch encounters the step potential
      ! mch = max(mch, istep + int(log(1.25_dp) / al))

! outward integration

      call vkboutwf(ll, nvkb, ee, vkb, evkb, rr, vloc, uu, up, node, mmax, mch)

      uout = uu(mch)
      upout = up(mch)

      ! integrate inward
      nin = mch + 2.3d0/al
      if (nin + 4 > mmax) nin = mmax - 4
      call lschvkbb_integrate_inward(ll, nin, mmax, mch, rr, vloc, cf, ee, uu, up, inflexion)

      ! if (inflexion) then
      !    write (*, *) 'INFLEXION POINT FOUND; SHOULD NOT HAPPEN'
      !    stop
      ! end if

      ! scale outside wf for continuity
      sc = uout/uu(mch)

      do ii = mch, nin
         up(ii) = sc*up(ii)
         uu(ii) = sc*uu(ii)
      end do

      ! do ii = mch, nin
      !    write(100 + nint, '(2f18.12)') rr(ii), uu(ii)
      ! end do

      upin = up(mch)

! perform normalization sum

      ro = rr(1)/dsqrt(amesh)
      sn = ro**(2*ll + 3)/dfloat(2*ll + 3)

      do ii = 1, nin - 3
         sn = sn + al*rr(ii)*uu(ii)**2
      end do

      sn = sn + al*(23.0d0*rr(nin - 2)*uu(nin - 2)**2 &
 &              + 28.0d0*rr(nin - 1)*uu(nin - 1)**2 &
 &              + 9.0d0*rr(nin)*uu(nin)**2)/24.0d0

! normalize u

      cn = 1.0d0/dsqrt(sn)
      uout = cn*uout
      upout = cn*upout
      upin = cn*upin

      do ii = 1, nin
         up(ii) = cn*up(ii)
         uu(ii) = cn*uu(ii)
      end do
      do ii = nin + 1, mmax
         uu(ii) = 0.0d0
      end do

      open (unit=1000 + nint, status='unknown')
      close (1000 + nint, status='delete')
      do ii = 1, mmax
         write (1000 + nint, '(2f18.12)') rr(ii), uu(ii)
      end do

      ! Count the number of nodes
      node = 0
      do ii = 1, mch - 1
         if (uu(ii) * uu(ii + 1) < 0.0_dp) then
            node = node + 1
         end if
      end do

      if (node - nn + ll + 1 - nextra_nodes == 0) then

         ! perturbation theory for energy shift
         de = 0.5d0*uout*(upout - upin)/(al*rr(mch))

         ! convergence test and possible exit
         if (dabs(de) < dmax1(dabs(ee), 0.2d0)*eps) then
            ierr = 0
            write (*, '(a)') ' converged'
            exit
         end if

         if (de > 0.0d0) then
            ! ee is too low
            write (*, '(a)') ' 2 ee is too low (de > 0); increasing by de'
            e_lower = ee
            e_lower_found = .true.
         else
            ! ee is too high
            write (*, '(a)') ' 3 ee is too high (de < 0); decreasing by de'
            e_upper = ee
         end if

         ! Update ee (capping de at 0.1)
         ee = ee + sign(min(dabs(de), 0.1), de)

         ! Make sure ee stays within the bounds
         if (ee > e_upper) then
            if (e_lower_found) then
               ee = 0.5d0*(e_upper + e_lower)
            else
               ! Cautiously inch upwards
               ee = 0.5_dp*(ee - de + e_upper)
            end if
         else if (e_lower_found) then
            if (ee < e_lower) ee = 0.5d0*(e_upper + e_lower)
         end if

      else if (node - nn + ll + 1 - nextra_nodes < 0) then
         ! too few nodes; ee is too low
         write (*, '(a,i1,a)') ' 1 ee is too low (', node, ' nodes); increasing'
         e_lower_found = .true.
         e_lower = ee
         ee = 0.5_dp*(ee + e_upper)
      else
         ! too many nodes; ee is too high
         write (*, '(a,i1,a)') ' 4 ee is too high (', node, ' nodes); decreasing'
         e_upper = ee
         if (e_lower_found) then
            ee = 0.5_dp*(e_upper + e_lower)
         else
            ! Cautiously inch downwards
            ee = ee - 0.1_dp
         end if
      end if

   end do
!fix sign to be positive at rr->oo
   if (uu(mch) < 0.0d0) then
      uu(:) = -uu(:)
      up(:) = -up(:)
   end if

   deallocate (cf)

   return

end subroutine robust_lschvkbb
