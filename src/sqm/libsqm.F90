#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
!
!        SQM  stand-alone quantum program
!        Shared library version by Xiaoliang Pan in 02/2022
!

module sqm

   use file_io_dat, only : MAX_FN_LEN

   implicit none

   _REAL_ escf
   _REAL_ born_radii(1000), one_born_radii(1000)
   _REAL_ intdiel, extdiel, Arad
   integer ntpr
   character(len=MAX_FN_LEN) mdin, mdout
   ! external charge
   _REAL_ excharge(40000)
   integer chgatnum(10000)
   integer :: igb, maxcyc
   _REAL_  :: grms_tol
   _REAL_  :: total_energy
   logical :: master=.true.

contains

subroutine sqm_init( natom, atnum, ncharge )

   use qmmm_module, only : qmmm_struct, qmmm_mpi, qm2_struct
   use sqm_qmmm_read_and_alloc, only : read_qmmm_nm_and_alloc

   implicit none

   integer, intent(inout) :: natom
   integer, intent(in)  :: atnum(*)
   integer, intent(in) :: ncharge

   integer ier

   ! ==== Initialise first_call flags for QMMM ====
   qmmm_struct%qm_mm_first_call = .true.
   qmmm_struct%fock_first_call = .true.
   qmmm_struct%fock2_2atm_first_call = .true.
   ! qmmm_struct%qm2_deriv_qm_analyt_first_call = .true.
   qmmm_struct%qm2_allocate_e_repul_first_call = .true.
   ! qmmm_struct%qm2_rotate_qmqm_first_call = .true.
   qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
   qmmm_struct%qm2_scf_first_call = .true.
   qmmm_struct%zero_link_charges_first_call = .true.
   qmmm_struct%adj_mm_link_pair_crd_first_call = .true.

   !     --- default file names ---

   mdin   = 'mdin'
   mdout  = 'mdout'
   igb = 0
   call amopen(5,mdin,'O','F','R')
   call amopen(6,mdout,'R','F','W')

   write(6,*) '           --------------------------------------------------------'
   write(6,*) '                            AMBER SQM VERSION 19'
   write(6,*) ''
   write(6,*) '                                    By'
   write(6,*) '             Ross C. Walker, Michael F. Crowley, Scott Brozell,'
   write(6,*) '                        Tim Giese, Andreas W. Goetz,'
   write(6,*) '                       Tai-Sung Lee and David A. Case'
   write(6,*) ''
   write(6,*) '           --------------------------------------------------------'
   write(6,*) ''

   call read_qmmm_nm_and_alloc(natom,igb,atnum,maxcyc,grms_tol,ntpr, &
                               ncharge,excharge,chgatnum )
   close(5)
   call qm_assign_atom_types

   ! Set default QMMM MPI parameters - for single cpu operation.
   ! These will get overwritten by qmmm_mpi_setup if MPI is on.
   ! qmmm_mpi%master = master
   qmmm_mpi%commqmmm_master = master
   qmmm_mpi%numthreads = 1
   qmmm_mpi%mytaskid = 0
   qmmm_mpi%natom_start = 1
   qmmm_mpi%natom_end = natom
   qmmm_mpi%nquant_nlink_start = 1
   qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink
   call allocate_qmgb(qmmm_struct%nquant_nlink)

   allocate( qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat = ier )
   REQUIRE(ier == 0)

   allocate ( qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat = ier )
   REQUIRE(ier == 0)

end subroutine sqm_init



subroutine sqm_calc( coords, extcharge )
   use qmmm_module, only : qmmm_nml, qmmm_struct, qm_gb, qmmm_mpi, qm2_struct, qmmm_scratch
   use qm2_dftb_module, only : ks_struct
   use constants, only : EV_TO_KCAL, KCAL_TO_EV, zero, one, alpb_alpha
   use qm2_pm6_hof_module, only : cct, nsp2, print, strlen

   use UtilitiesModule, only : print

   implicit none

   _REAL_, intent(in) :: coords(qmmm_struct%nquant_nlink*3)
   _REAL_, intent(in) :: extcharge(qmmm_struct%qm_mm_pairs*4)


!Locals
   _REAL_ :: alpb_beta

   integer :: ier=0
   integer i, i3, j, natom, iqmp

   character(len=strlen) :: string

   natom = qmmm_struct%nquant_nlink + qmmm_struct%qm_mm_pairs
   qmmm_struct%qm_xcrd = 0.0D0
   j=0
   do i=1,qmmm_struct%qm_mm_pairs
      qmmm_struct%qm_xcrd(1,i) = extcharge(j+1)
      qmmm_struct%qm_xcrd(2,i) = extcharge(j+2)
      qmmm_struct%qm_xcrd(3,i) = extcharge(j+3)
      qmmm_struct%qm_xcrd(4,i) = extcharge(j+4)
      j=j+4
   end do

!=============================================================================
!                   START OF QMMM SETUP: allocate list memory
!=============================================================================

!  If this is the first call to the routine, do some initial allocation
!  that has not been done elsewhere.
   if (qmmm_struct%qm_mm_first_call) then

     allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant_nlink), stat=ier )
                !Stores the REAL and link atom qm coordinates
     REQUIRE(ier == 0)

     !Allocation for QM_GB (qmgb==2)
     if (qmmm_nml%qmgb == 2) then
       !Calculate dielectric factor
       if (qm_gb%alpb_on) then
         alpb_beta=alpb_alpha*(intdiel/extdiel)
         qm_gb%intdieli = one/(intdiel*(one + alpb_beta))
         qm_gb%extdieli = one/(extdiel*(one + alpb_beta))
         qm_gb%one_Arad_beta = alpb_beta/Arad
       else
         qm_gb%intdieli = 1.0d0/intdiel
         qm_gb%extdieli = 1.0d0/extdiel
       end if
       qm_gb%mmcut2 = 999.d0
     end if
   end if ! ---- first call endif ----------

   ! call qm_extract_coords(coords)
   i3 = 0
   !do i=1,natom
   do i=1,qmmm_struct%nquant_nlink
      qmmm_struct%qm_coords(1,i) = coords(i3+1)
      qmmm_struct%qm_coords(2,i) = coords(i3+2)
      qmmm_struct%qm_coords(3,i) = coords(i3+3)
      i3 = i3 + 3
   end do

!=============================================================================
!                   START OF REST OF QMMM SETUP
!=============================================================================
   if(qmmm_struct%qm_mm_first_call) then
       if (qmmm_mpi%commqmmm_master) then
          write(6,'(/80(1H-)/''  QM CALCULATION INFO'',/80(1H-))')
       end if

       call qm2_load_params_and_allocate(.false.) !Load the parameters
             !Also does a lot of memory allocation and pre-calculates all
             !the STO-6G orbital expansions.

       if (qmmm_mpi%commqmmm_master) then
          ! call qm_print_dyn_mem(natom,qmmm_struct%qm_mm_pairs)
          call qm_print_coords(0,.true.)
          !Finally print the result header that was skipped in sander.
          write(6,'(/80(1H-)/''  RESULTS'',/80(1H-)/)')
       end if
   end if !if (qmmm_struct%qm_mm_first_call)

!======================END OF QMMM SETUP ======================================

   !Calculate RIJ and many related equations here. Necessary memory allocation
   !is done inside the routine.
!Parallel
   call qm2_calc_rij_and_eqns(qmmm_struct%qm_coords, qmmm_struct%nquant_nlink, &
          qmmm_struct%qm_xcrd, natom, qmmm_struct%qm_mm_pairs)
                                !and store them in memory to save time later.

   !============================
   ! Calculate SCF Energy
   !============================
   call qm2_energy(escf, qm2_struct%scf_mchg, natom, born_radii, one_born_radii)

   !=============================
   !   Print Mulliken Charges
   !=============================

   if (qmmm_nml%printcharges .and. qmmm_mpi%commqmmm_master) then
     call qm2_print_charges(1,qmmm_nml%dftb_chg,qmmm_struct%nquant_nlink, &
                            qm2_struct%scf_mchg,qmmm_struct%iqm_atomic_numbers)
   end if
   !=============================
   !   Print Dipole Charges
   !=============================

   select case (qmmm_nml%printdipole)
      case (1)
         call qm2_calc_dipole(coords)
      case (2)
         write (6,'("QMMM: Not MM part; please check your selection")')
      case default
   end select

   ! Print some extra informatiom about energy contributions
   ! (This is really only required in sander since we print the energies anyways
   !  but kept here for historical reasons)

   if (qmmm_mpi%commqmmm_master) then
      call qm2_print_energy(qmmm_nml%verbosity, qmmm_nml%qmtheory, escf, qmmm_struct)
   end if

   qmmm_struct%qm_mm_first_call = .false.

   ! --------------------
   ! print MO eigenvalues
   ! --------------------
   if ( (qmmm_nml%print_eigenvalues > 0) .and. qmmm_mpi%commqmmm_master) then
      write (6,*) ''
      if (qmmm_nml%qmtheory%DFTB) then
         call print(' Final MO eigenvalues (au)', ks_struct%ev(1:ks_struct%ind(qmmm_struct%nquant_nlink+1)))
      else
         call print(' Final MO eigenvalues (eV)', qm2_struct%eigen_values)
      end if
   end if

   ! ----------------
   ! print SCF energy
   ! ----------------
   string = 'CC triple bond correction (unpublished)'
   call print(cct, string)
   string = 'Nitrogen pyramidalization correction'
   call print(nsp2, string)
   ! at present sqm does only pure QM, need to update this for QM/MM
   total_energy = qmmm_struct%elec_eng +  qmmm_struct%enuclr_qmqm + qmmm_struct%enuclr_qmmm
   write(6,'(/,a,f20.8,a,f18.8,a)') ' Heat of formation   =', &
        escf, ' kcal/mol  (', escf*KCAL_TO_EV, ' eV)'
   write(6,'(/a,f20.8,a,f18.8,a)')   ' Total SCF energy    =', &
        total_energy*EV_TO_KCAL, ' kcal/mol  (', total_energy, ' eV)'
   write(6,'(a,f20.8,a,f18.8,a)')   ' Electronic energy   =', &
        qmmm_struct%elec_eng*EV_TO_KCAL, ' kcal/mol  (', qmmm_struct%elec_eng, ' eV)'
   write(6,'(a,f20.8,a,f18.8,a)')   ' Core-core repulsion =', &
        (qmmm_struct%enuclr_qmqm+qmmm_struct%enuclr_qmmm)*EV_TO_KCAL, ' kcal/mol  (', &
        (qmmm_struct%enuclr_qmqm+qmmm_struct%enuclr_qmmm), ' eV)'
    if (qmmm_nml%qmtheory%DISPERSION .or. qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
       write(6,'(/a,f20.8,a,f18.8,a)')  ' Dispersion energy   =', &
            qmmm_struct%dCorrection, ' kcal/mol  (', qmmm_struct%dCorrection*KCAL_TO_EV, ' eV)'
    end if
    if (qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
       write(6,'(a,f20.8,a,f18.8,a)')   ' H-bond energy       =', &
            qmmm_struct%hCorrection, ' kcal/mol  (', qmmm_struct%hCorrection*KCAL_TO_EV, ' eV)'
    end if

   write(6,*) ''

   if ( qmmm_nml%printbondorders ) then
      write(6,*) ''
      write(6,*) 'Bond Orders'
      call qm2_print_bondorders()
   end if

  qmmm_struct%dxyzqm=zero
   if (qmmm_nml%qmtheory%DFTB) then
     call qm2_dftb_get_qm_forces(qmmm_struct%dxyzqm)
   else
     !standard semi-empirical
     call qm2_get_qm_forces(qmmm_struct%dxyzqm)
   end if

   !Calculate forces between QM and MM atoms
   if (qmmm_nml%qmmm_int > 0 .and. (qmmm_nml%qmmm_int /= 5) ) then
     qmmm_struct%dxyzcl=zero
     qmmm_struct%mm_esp=zero
     if (qmmm_nml%qmtheory%DFTB) then
       iqmp = qmmm_struct%qm_mm_pairs
       call qm2_dftb_get_qmmm_forces(qmmm_struct%dxyzcl, &
              qmmm_struct%dxyzqm, qmmm_scratch%qm_real_scratch, &
              qmmm_scratch%qm_real_scratch(natom+1:natom+iqmp), &
              qmmm_scratch%qm_real_scratch(2*natom+1:2*natom+iqmp), &
              qmmm_scratch%qm_real_scratch(3*natom+1:3*natom+iqmp))
     else
       call qm2_get_qmmm_forces(qmmm_struct%dxyzqm, qmmm_struct%qm_xcrd, &
                                qmmm_struct%dxyzcl, qm2_struct%scf_mchg, qmmm_struct%mm_esp)
     end if
   end if

   write(6,*)

   write(6,*) '          --------- Calculation Completed ----------'
   write(6,*)
   close(6)

end subroutine sqm_calc


subroutine sqm_clean
   use qmmm_module, only : qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params, deallocate_qmmm
   implicit none
   call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
end subroutine sqm_clean

end module sqm

!======================END OF QM_MM ======================================


subroutine qm_gradient(output)
   use qmmm_module, only : qmmm_struct
   implicit none
   _REAL_, intent(out) :: output(*)
   integer :: i, i3
   i3 = 0
   do i=1,qmmm_struct%nquant_nlink
      output(i3+1) = qmmm_struct%dxyzqm(1,i)
      output(i3+2) = qmmm_struct%dxyzqm(2,i)
      output(i3+3) = qmmm_struct%dxyzqm(3,i)
      i3 = i3 + 3
   end do
end subroutine qm_gradient


subroutine mm_esp(output)
   use qmmm_module, only : qmmm_struct
   implicit none
   _REAL_, intent(out) :: output(*)
   integer :: i, i4
   i4 = 0
   do i=1,qmmm_struct%qm_mm_pairs
      output(i4+1) = qmmm_struct%mm_esp(1,i)
      output(i4+2) = qmmm_struct%mm_esp(2,i)
      output(i4+3) = qmmm_struct%mm_esp(3,i)
      output(i4+4) = qmmm_struct%mm_esp(4,i)
      i4 = i4 + 4
   end do
end subroutine mm_esp


subroutine qm_mulliken(output)
   use qmmm_module, only : qmmm_struct, qm2_struct
   implicit none
   _REAL_, intent(out) :: output(*)
   integer :: i
   do i=1,qmmm_struct%nquant_nlink
      output(i) = qm2_struct%scf_mchg(i)
   end do
end subroutine qm_mulliken


subroutine den_matrix(output)
   use qmmm_module, only : qm2_struct
   implicit none
   _REAL_, intent(out) :: output(*)
   integer :: i
   do i=1,qm2_struct%matsize
      output(i) = qm2_struct%den_matrix(i)
   end do
end subroutine den_matrix
