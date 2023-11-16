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
! self-consistent pseudoatom calculation

subroutine psatom(na, la, ea, fat, nv, it, rhoc, rho, &
                  & rr, rcmax, mmax, mxprj, iexc, etot, nproj, vpuns, lloc, vkb, evkb, plot, robust, label, istep, ierr)

!na  principal quantum number array, dimension nv
!la  angular-momenta
!ea  eigenvalues (input starting guess, output)
!fa  occupancies
!rpk  radius of outermost peak of wave function
!nv  number of valence states
!it  number of iterations (output)
!rr  log radial mesh
!rcmax  maximum core radius for psp
!mmax  size of log grid
!mxprj  dimension of number of projectors
!iexc  exchange-correlation function to be used
!etot  pseudoatom total energy (output)
!nproj  number of VKB projectora to use for each l
!vpuns  unscreened semi-local pseudopotentials (plus differenv vloc if lloc==4)
!lloc  index-1 of local potential
!vkb   Vanderbilt-Kleinman-Bylander projectors
!evkb VKB projector coefficients
!uua Array of pseudo-atomic orbitals (output)

   implicit none
   integer, parameter :: dp = kind(1.0d0)

!Input variables

   integer :: mmax, mxprj, iexc, nv, lloc, okb
   integer :: na(nv), la(nv), nproj(5)
   real(dp) :: rcmax
   real(dp) :: fat(30, 2), rr(mmax)
   real(dp) :: vpuns(mmax, 5), vkb(mmax, mxprj, 4), evkb(mxprj, 4)
   logical :: plot
   character(len=*) :: label
   integer :: istep
   logical :: robust

!Output variables
   integer :: it
   real(dp) :: etot
   real(dp) :: ea(nv)
   real(dp) :: rho(mmax), rhoc(mmax), vi(mmax)

!Local variables
   integer :: nin, mch
   real(dp) :: amesh, al
   real(dp) :: dr, eeel, eexc, et, emin, emax, rl, rl1, sd, sn, sls, eeig
   real(dp) :: thl, vn, zval, dfa
   real(dp) :: fa(30)
   integer :: ii, jj, l1, ierr, icx, nprj, kk, lextra, nextra_nodes
   logical :: convg

   real(dp), allocatable :: uua(:, :)
   real(dp), allocatable :: up(:)
   real(dp), allocatable :: vo(:), vi1(:), vo1(:), vxc(:), vp(:)
   real(dp), allocatable :: vtot(:)

! blend parameter for Anderson iterative potential mixing
   real(dp), parameter ::  bl = 0.5d0

   allocate (uua(mmax, nv), up(mmax))
   allocate (vo(mmax), vi1(mmax), vo1(mmax), vxc(mmax), vp(mmax))
   allocate (vtot(mmax))

! this seems necessary in sratom, so I might as well do it  here too
! don't ask why!
   uua(:, :) = 0.d0; up(:) = 0.d0
   vo(:) = 0.d0; vi1(:) = 0.d0; vo1(:) = 0.d0; vxc(:) = 0.d0; vp(:) = 0.d0
   vtot(:) = 0.d0
   dr = 0.d0; eeel = 0.d0; eexc = 0.d0; et = 0.d0; emin = 0.d0; emax = 0.d0; 
   rl = 0.d0; rl1 = 0.d0; sd = 0.d0; sn = 0.d0; sls = 0.d0; eeig = 0.d0; 
   thl = 0.d0; vn = 0.d0; zval = 0.d0; 
   nin = 0; mch = 0

   al = 0.01d0*dlog(rr(101)/rr(1))
   amesh = dexp(al)

   nprj = 0
   do l1 = 1, 4
      nprj = max(nprj, nproj(l1))
   end do

! total valence charge, initial config
   zval = 0.0d0
   do ii = 1, nv
      fa(ii) = fat(ii, 1)
      zval = zval + fa(ii)
   end do

! test if new states are added
   icx = 50
!do ii=1,nv
!  if(fat(ii,2)>fat(ii,1)) icx=100
!end do

! screening potential for pseudocharge
! input rho is assumed to be valence rho from all-electron calculation
! which shouldn''t be a bad approximation for a starting screened
! potential

   call vout(1, rho, rhoc, vi, vxc, zval, eeel, eexc, rr, mmax, iexc)

! big self  self-consietency loop
   do it = 1, 100
      write(*, *) 'psatom iteration step ', it

! convergence is only to be considered if we are fully in the target
! configurtion

      if (it > 2) then
         convg = .true.
      else
         convg = .false.
      end if

! solve for bound states in turn

      eeig = 0.0d0

      vtot(:) = vpuns(:, lloc + 1) + vi(:)

      rho(:) = 0.0d0

      do ii = 1, nv

         if (fa(ii) == 0.0d0) then
            cycle
         end if
         et = ea(ii)
         ierr = 0
         l1 = la(ii) + 1

         emin = 1.10d0*ea(ii)
         emax = 0.90d0*ea(ii)
         if (robust) then
            do nextra_nodes = 0, 3
               write(*, '(a,i1,a)') 'Looking for a solution with ', na(ii) - la(ii) - 1 + nextra_nodes, ' nodes'
               call robust_lschvkbb(na(ii), la(ii), nproj(l1), ierr, et, &
    &                       rr, vtot, vkb(1, 1, l1), evkb(1, l1), uua(:, ii), up, mmax, mch, istep, nextra_nodes)
               if (ierr == 0) then
                  exit
               end if
            end do
         else
            call lschvkbb(na(ii), la(ii), nproj(l1), ierr, et, emin, emax, &
    &                    rr, vtot, vkb(1, 1, l1), evkb(1, l1), uua(:, ii), up, mmax, mch)
         end if

         if (ierr .ne. 0) then
            write (6, '(/a,3i4)') 'psatom: WARNING lschvkbb convergence error n,l,iter=', &
     &       na(ii), la(ii), it
            write (6, '(a,1p,4e14.6)') 'ea(ii),et,emin,emax', ea(ii), et, emin, emax
            exit
         end if

! overall convergence criterion based on eps within lschf*
         if (ea(ii) /= et) convg = .false.
         ea(ii) = et

! accumulate charge and eigenvalues
         eeig = eeig + fa(ii)*ea(ii)
         rho(:) = rho(:) + fa(ii)*(uua(:, ii)/rr(:))**2

      end do !ii

      if (ierr /= 0) exit

! total valence charge
      zval = 0.0d0
      do ii = 1, nv
         zval = zval + fa(ii)
      end do

! output potential
      call vout(1, rho, rhoc, vo, vxc, zval, eeel, eexc, rr, mmax, iexc)

! generate next iteration using d. g. anderson''s
! method
      thl = 0.0d0
      if (it > icx + 1) then
         sn = 0.0d0
         sd = 0.0d0
         do ii = 1, mmax
            rl = vo(ii) - vi(ii)
            rl1 = vo1(ii) - vi1(ii)
            dr = rl - rl1
            sn = sn + rl*dr*rr(ii)**2
            sd = sd + dr*dr*rr(ii)**2
         end do
         thl = sn/sd
      end if

      do ii = 1, mmax
         vn = (1.0d0 - bl)*((1.0d0 - thl)*vi(ii) + thl*vi1(ii)) &
      &   + bl*((1.0d0 - thl)*vo(ii) + thl*vo1(ii))
         vi1(ii) = vi(ii)
         vo1(ii) = vo(ii)
         vi(ii) = vn
      end do

      if (convg) exit

      if (it == 400 .and. .not. convg) then
         write (6, '(/a)') 'psatom: potential failed to converge'
         ierr = 101
      end if

! switch from reference to new configuration
! new or increased occupancy added in second stage
      do ii = 1, nv
         dfa = 0.02d0*(fat(ii, 2) - fat(ii, 1))
         if (it <= 50) then
            fa(ii) = 0.02d0*(50 - it)*fat(ii, 1) + 0.02d0*it*fat(ii, 2)
         else
            fa(ii) = fat(ii, 2)
         end if
      end do

   end do !it

! total energy output

! output potential for e-e interactions

   call vout(1, rho, rhoc, vo, vxc, zval, eeel, eexc, &
  &          rr, mmax, iexc)

   etot = eeig + eexc - 0.5d0*eeel

   ! Plotting
   if (plot) then
      open(20, file='paos' // label // '.dat', status='unknown')
      ! Header
      write(20, '(2i)') mmax, nv
      do kk = 1, nv
         write(20, '(i2)', advance='no') la(kk)
      end do

      ! Header
      ! write(20, '(a6,a10,a16)', advance='no') '# PAOs', 'r', 'rho'
      ! do kk = 1, nv
      !    write(20, '(a15,i1)', advance='no') 'psi', kk
      ! end do
      ! do kk = 1, 5
      !    write(20, '(a14,i1,a)', advance='no') 'v(l=', kk, ')'
      ! end do
      write(20, '()')

      ! Content
      do ii = 1, mmax
         write(20, '(2f20.15)', advance='no') log(rr(ii)), rr(ii)
         do kk = 1, nv
            write(20, '(f20.15)', advance='no') uua(ii, kk)
         end do
         ! do kk = 1, 5
         !    write(20, '(f16.8)', advance='no') vpuns(ii, kk)
         ! end do
         write(20, '()')
      end do
      close(20)
   end if

   deallocate (uua, up)
   deallocate (vo, vi1, vo1, vxc, vp)
   deallocate (vtot)
   return

end subroutine psatom
