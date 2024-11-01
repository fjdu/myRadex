! This file is part of our code to solve the statistical equilibrium problem of
! the energy level population of a uniform volume of gas under the effect of
! background radiation and collisional excitation.
!
! The code works similar to radex
! (see https://sron.rug.nl/~vdtak/radex/index.shtml)
!
! One difference from radex is that we solve the statistical equilibrium
! equation using an ode solver (ODEPACK, see netlib.org/odepack/).
!
! Written by Fujun Du, fujun.du@gmail.com, fdu@umich.edu
!
! 2014-01-02 Thu 02:43:39
! 2014-01-02 Thu 22:58:37
!


module statistic_equilibrium

use trivials
use phy_const
implicit none

integer, parameter, private :: const_len_energy_level = 12
integer, parameter, private :: const_len_molecule = 12


type :: type_energy_level
  character(len=const_len_energy_level) :: name_energy
  integer id
  double precision :: energy
  double precision :: weight
end type type_energy_level


type :: type_rad_transition
  double precision Eup, Elow, dE, freq, lambda
  double precision Aul, Bul, Blu, beta, J_ave, J_cont_bg, cooling_rate
  double precision tau
  double precision, dimension(:), allocatable :: critical_densities
  integer iup, ilow
end type type_rad_transition


type :: type_collisional_transition
  character(len=const_len_molecule) :: name_partner
  double precision dens_partner
  integer n_transition, n_T
  integer, dimension(:), allocatable :: iup, ilow
  double precision, dimension(:), allocatable :: T_coll
  double precision, dimension(:,:), allocatable :: Cul
end type type_collisional_transition


type :: type_rad_set
  integer n_transition
  type(type_rad_transition), dimension(:), allocatable :: list
end type type_rad_set


type :: type_colli_set
  integer n_partner
  type(type_collisional_transition), dimension(:), allocatable :: list
end type type_colli_set


type :: type_molecule_energy_set
  character(len=const_len_molecule) name_molecule
  character(len=16) :: geotype = ''
  integer iSpe, iType
  double precision Tkin, density_mol, dv, length_scale, cooling_rate_total
  integer n_level
  type(type_energy_level), dimension(:), allocatable :: level_list
  double precision, dimension(:), allocatable :: f_occupation
  type(type_rad_set), allocatable :: rad_data
  type(type_colli_set), allocatable :: colli_data
  double precision :: abundance_factor = 1D0
end type type_molecule_energy_set


type :: type_statistic_equil_params
  integer nitem
  double precision :: RTOL = 1D-4, ATOL = 1D-20
  double precision :: t_max = 1D6, dt_first_step = 1D-4, ratio_tstep = 1.1D0
  real :: max_runtime_allowed = 10.0
  integer n_record
  integer :: &
        NERR, &
        NEQ, &
        ITOL = 1, &
        ITASK = 1, &
        ISTATE = 1, &
        IOPT = 1, &
        LIW, &
        LRW, &
        MF = 21
  double precision, dimension(:), allocatable :: RWORK
  integer, dimension(:), allocatable :: IWORK
  logical is_good
end type type_statistic_equil_params


type :: type_continuum_lut
  integer :: n=0
  double precision, dimension(:), allocatable :: lam, alpha, J
end type type_continuum_lut

type(type_molecule_energy_set), pointer :: a_mol_using => null()

type(type_statistic_equil_params) statistic_equil_params

type(type_continuum_lut) cont_lut_bg, cont_lut_in, cont_lut_out


contains


subroutine reset_statistic_equil_params
  statistic_equil_params%is_good = .true.
  statistic_equil_params%NERR = 0
  statistic_equil_params%ITASK = 1
  statistic_equil_params%ISTATE = 1
  statistic_equil_params%IOPT = 1
  !
  statistic_equil_params%RWORK = 0D0
  statistic_equil_params%IWORK = 0
  statistic_equil_params%IWORK(6) = 5000
end subroutine reset_statistic_equil_params


subroutine calc_cooling_rate
  integer i
  a_mol_using%cooling_rate_total = 0D0
  do i=1, a_mol_using%rad_data%n_transition
    associate( &
          n_mol => a_mol_using%density_mol, &
          beta  => a_mol_using%rad_data%list(i)%beta, &
          nu    => a_mol_using%rad_data%list(i)%freq, &
          Aul   => a_mol_using%rad_data%list(i)%Aul, &
          Bul   => a_mol_using%rad_data%list(i)%Bul, &
          Blu   => a_mol_using%rad_data%list(i)%Blu, &
          iup   => a_mol_using%rad_data%list(i)%iup, &
          ilow  => a_mol_using%rad_data%list(i)%ilow, &
          J_ave => a_mol_using%rad_data%list(i)%J_ave)
      a_mol_using%rad_data%list(i)%cooling_rate = &
        beta * phy_hPlanck_CGS * nu * n_mol * &
        ((Aul + Bul * J_ave) * a_mol_using%f_occupation(iup) - &
          Blu * J_ave * a_mol_using%f_occupation(ilow))
      a_mol_using%cooling_rate_total = a_mol_using%cooling_rate_total + &
        a_mol_using%rad_data%list(i)%cooling_rate
    end associate
  end do
end subroutine calc_cooling_rate


subroutine load_moldata_LAMBDA(filename, recalculateFreqWithEupElow)
  character(len=*) filename
  logical, intent(in), optional :: recalculateFreqWithEupElow
  logical recalcFreq
  character(len=512) strtmp
  integer, parameter :: nstr_split = 64
  character(len=32), dimension(nstr_split) :: str_split
  integer i, j, k, fU, nout
  character(len=8), parameter :: strfmt_row = '(A512)'
  character(len=8), parameter :: strfmt_float = '(F16.3)'
  character(len=8), parameter :: strfmt_int = '(I6)'
  character(len=8), parameter :: strfmt_int_long = '(I64)'
  double precision, parameter :: GHz_to_Hz = 1D9
  !
  integer iup, ilow
  !
  integer n_T_, n_transition_
  !
  recalcFreq = .false.
  if (present(recalculateFreqWithEupElow)) then
    recalcFreq = recalculateFreqWithEupElow
  end if
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a free file unit.  In load_moldata_LAMBDA.'
    stop
  end if
  call openFileSequentialRead(fU, filename, 99999)
  ! Get molecule name.
  read(fU,'(A1)') strtmp
  read(fU, strfmt_row) strtmp
  call split_str_by_space(strtmp, str_split, nstr_split, nout)
  a_mol_using%name_molecule = trim(str_split(1))
  ! Get energy level list
  read(fU,'(A1)') strtmp
  read(fU,'(A1)') strtmp
  read(fU,'(A1)') strtmp
  read(fU, strfmt_int_long) a_mol_using%n_level
  read(fU,'(A1)') strtmp
  write(*, '(A,I6)') 'Number of energy levels:', a_mol_using%n_level
  allocate(a_mol_using%level_list(a_mol_using%n_level), &
           a_mol_using%f_occupation(a_mol_using%n_level))
  do i=1, a_mol_using%n_level
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    read(str_split(2), strfmt_float) a_mol_using%level_list(i)%energy
    read(str_split(3), strfmt_float) a_mol_using%level_list(i)%weight
  end do
  !
  ! Get radiative transitions
  allocate(a_mol_using%rad_data)
  read(fU,'(A1)') strtmp
  read(fU, strfmt_int_long) a_mol_using%rad_data%n_transition
  write(*,'(A,I6)') 'Number of radiative transitions:', a_mol_using%rad_data%n_transition
  read(fU,'(A1)') strtmp
  allocate(a_mol_using%rad_data%list(a_mol_using%rad_data%n_transition))
  do i=1, a_mol_using%rad_data%n_transition
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    read(str_split(2), strfmt_int) a_mol_using%rad_data%list(i)%iup
    read(str_split(3), strfmt_int) a_mol_using%rad_data%list(i)%ilow
    read(str_split(4), strfmt_float) a_mol_using%rad_data%list(i)%Aul
    read(str_split(5), strfmt_float) a_mol_using%rad_data%list(i)%freq
    !read(str_split(6), strfmt_float) a_mol_using%rad_data%list(i)%Eup
    !
    iup  = a_mol_using%rad_data%list(i)%iup
    ilow = a_mol_using%rad_data%list(i)%ilow
    !
    if (recalcFreq) then
      ! The frequencies in the LAMDA data file may be inconsistent with the energy levels,
      ! so here I recompute them using the energy differences.  The results are in Hz.
      ! However, note that the energy levels in the LAMDA data file may not be accurate,
      ! so by default this recalculation is NOT done.
      a_mol_using%rad_data%list(i)%freq = phy_SpeedOfLight_CGS * &
        (a_mol_using%level_list(iup)%energy - &
         a_mol_using%level_list(ilow)%energy)
    else
      ! The frequency unit in the LAMDA file is in GHz. Now convert to Hz.
      a_mol_using%rad_data%list(i)%freq = a_mol_using%rad_data%list(i)%freq * GHz_to_Hz
    end if
    a_mol_using%rad_data%list(i)%Eup  = a_mol_using%level_list(iup)%energy * phy_cm_1_2K
    a_mol_using%rad_data%list(i)%Elow = a_mol_using%level_list(ilow)%energy * phy_cm_1_2K
  end do
  !
  ! Convert the energy unit into Kelvin from cm-1
  a_mol_using%level_list%energy = a_mol_using%level_list%energy * phy_cm_1_2K
  !
  ! Lambda in micron
  a_mol_using%rad_data%list%lambda = phy_SpeedOfLight_SI/a_mol_using%rad_data%list%freq*1D6
  !
  a_mol_using%rad_data%list%Bul = a_mol_using%rad_data%list%Aul / &
    ((2D0*phy_hPlanck_CGS/phy_SpeedOfLight_CGS**2) * &
     (a_mol_using%rad_data%list%freq)**3)
  do i=1, a_mol_using%rad_data%n_transition
    j = a_mol_using%rad_data%list(i)%iup
    k = a_mol_using%rad_data%list(i)%ilow
    a_mol_using%rad_data%list(i)%Blu = a_mol_using%rad_data%list(i)%Bul * &
        a_mol_using%level_list(j)%weight / a_mol_using%level_list(k)%weight
  end do
  !
  ! Collisional transitions
  allocate(a_mol_using%colli_data)
  ! Get the number of collisional partners
  read(fU,'(A1)') strtmp
  read(fU, strfmt_int_long) a_mol_using%colli_data%n_partner
  write(*,'(A,I6)') 'Number of collisional partners:', a_mol_using%colli_data%n_partner
  allocate(a_mol_using%colli_data%list(a_mol_using%colli_data%n_partner))
  !
  do i=1, a_mol_using%rad_data%n_transition
    allocate(a_mol_using%rad_data%list(i)%critical_densities(a_mol_using%colli_data%n_partner))
  end do
  !
  do i=1, a_mol_using%colli_data%n_partner
    ! Get the name of partner
    read(fU,'(A1)') strtmp
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    a_mol_using%colli_data%list(i)%name_partner = trim(str_split(4))
    if (a_mol_using%colli_data%list(i)%name_partner .eq. 'electron') then
      a_mol_using%colli_data%list(i)%name_partner = 'e'
    else if (a_mol_using%colli_data%list(i)%name_partner .eq. 'with') then
      a_mol_using%colli_data%list(i)%name_partner = str_split(5)
    end if
    write(*, '(A, I6, " ", A)') 'Collisional partner:', i, a_mol_using%colli_data%list(i)%name_partner
    ! Get the number of transitions and temperatures
    read(fU,'(A1)') strtmp
    read(fU, strfmt_int_long) a_mol_using%colli_data%list(i)%n_transition
    write(*,'(A,2I6)') 'Number of collisional transitions:', i, a_mol_using%colli_data%list(i)%n_transition
    read(fU,'(A1)') strtmp
    read(fU, strfmt_int_long) a_mol_using%colli_data%list(i)%n_T
    write(*,'(A,2I6)') 'Number of temperatures:', i, a_mol_using%colli_data%list(i)%n_T
    !
    ! Name too long...
    n_transition_ = a_mol_using%colli_data%list(i)%n_transition
    n_T_ = a_mol_using%colli_data%list(i)%n_T
    if ((n_T_+3) .gt. nstr_split) then
      write(*,'(A)') 'The number of different temperatures is too large!'
      write(*,'(A, I6)') 'nstr_split = ', nstr_split
      write(*,'(A)') 'Change nstr_split of the source code to a higher value.'
      stop
    end if
    !
    allocate(a_mol_using%colli_data%list(i)%iup(n_transition_), &
             a_mol_using%colli_data%list(i)%ilow(n_transition_), &
             a_mol_using%colli_data%list(i)%T_coll(n_T_), &
             a_mol_using%colli_data%list(i)%Cul(n_T_, n_transition_))
    !
    ! Get the list of temperatures
    read(fU,'(A1)') strtmp
    read(fU, strfmt_row) strtmp
    call split_str_by_space(strtmp, str_split, nstr_split, nout)
    do j=1, n_T_
      read(str_split(j), strfmt_float) a_mol_using%colli_data%list(i)%T_coll(j)
    end do
    ! Get the collision coefficients
    read(fU,'(A1)') strtmp
    do j=1, n_transition_
      read(fU, strfmt_row) strtmp
      call split_str_by_space(strtmp, str_split, nstr_split, nout)
      read(str_split(2), strfmt_int) a_mol_using%colli_data%list(i)%iup(j)
      read(str_split(3), strfmt_int) a_mol_using%colli_data%list(i)%ilow(j)
      do k=1, n_T_
        read(str_split(3+k), strfmt_float) a_mol_using%colli_data%list(i)%Cul(k, j)
        !iup_ = a_mol_using%colli_data%list(i)%iup(j)
        !ilow_ = a_mol_using%colli_data%list(i)%ilow(j)
      end do
    end do
  end do
  !
  close(fU)
  ! Test the results
  ! write(*,*) a_mol_using%name_molecule
  ! do i=1, a_mol_using%n_level
  !   write(*,*) i, a_mol_using%level_list(i)%energy, a_mol_using%level_list(i)%weight
  ! end do
  ! do i=1, a_mol_using%rad_data%n_transition
  !   write(*,*) i, a_mol_using%rad_data%list(i)%iup, &
  !     a_mol_using%rad_data%list(i)%ilow, &
  !     a_mol_using%rad_data%list(i)%Aul, &
  !     a_mol_using%rad_data%list(i)%Bul, &
  !     a_mol_using%rad_data%list(i)%Blu, &
  !     a_mol_using%rad_data%list(i)%freq, &
  !     a_mol_using%rad_data%list(i)%Eup
  ! end do
  ! do i=1, a_mol_using%colli_data%n_partner
  !   write(*,*) i, a_mol_using%colli_data%list(i)%name_partner
  !   write(*,*) a_mol_using%colli_data%list(i)%T_coll
  !   do j=1, a_mol_using%colli_data%list(i)%n_transition
  !     write(*,*) i, j, a_mol_using%colli_data%list(i)%iup(j), &
  !       a_mol_using%colli_data%list(i)%ilow(j), &
  !       a_mol_using%colli_data%list(i)%Cul(:, j)
  !   end do
  ! end do
end subroutine load_moldata_LAMBDA



subroutine deallocate_a_mol_using()
  if (associated(a_mol_using)) then
    deallocate(a_mol_using%rad_data%list)
    deallocate(a_mol_using%colli_data%list)
    deallocate(a_mol_using%level_list, &
               a_mol_using%f_occupation, &
               a_mol_using%rad_data, &
               a_mol_using%colli_data)
    deallocate(a_mol_using)
    nullify(a_mol_using)
  end if
end subroutine deallocate_a_mol_using


subroutine statistic_equil_solve
  use my_timer
  external stat_equili_ode_f, stat_equili_ode_jac
  integer i
  double precision t, tout, t_step, tscal, edot
  type(atimer) timer
  real time_thisstep, runtime_thisstep, time_laststep, runtime_laststep
  !
  call reset_statistic_equil_params
  !
  t = 0D0
  tout = statistic_equil_params%dt_first_step
  t_step = statistic_equil_params%dt_first_step
  !
  call timer%init('Stati_equil')
  time_laststep = timer%elapsed_time()
  runtime_laststep = huge(0.0)
  !
  statistic_equil_params%n_record = ceiling( &
    log(statistic_equil_params%t_max / statistic_equil_params%dt_first_step * &
        (statistic_equil_params%ratio_tstep - 1D0) + 1D0) &
    / log(statistic_equil_params%ratio_tstep))
  !
  do i=2, statistic_equil_params%n_record
    !write (*, '(A, 25X, "Solving... ", I5, " (", F5.1, "%)", "  t = ", ES9.2, "  tStep = ", ES9.2)') &
    !  CHAR(27)//'[A', i, real(i*100)/real(statistic_equil_params%n_record), t, t_step
    !
    call DLSODE( &
         stat_equili_ode_f, &
         !
         a_mol_using%n_level, &
         a_mol_using%f_occupation, &
         !
         t, &
         tout, &
         !
         statistic_equil_params%ITOL, &
         statistic_equil_params%RTOL, &
         statistic_equil_params%ATOL, &
         statistic_equil_params%ITASK, &
         statistic_equil_params%ISTATE, &
         statistic_equil_params%IOPT, &
         statistic_equil_params%RWORK, &
         statistic_equil_params%LRW, &
         statistic_equil_params%IWORK, &
         statistic_equil_params%LIW, &
         !
         stat_equili_ode_jac, &
         !
         statistic_equil_params%MF)
    !
    time_thisstep = timer%elapsed_time()
    runtime_thisstep = time_thisstep - time_laststep
    if ((runtime_thisstep .gt. &
         max(100D0*runtime_laststep, &
             2D0*statistic_equil_params%max_runtime_allowed)) &
        .or. &
        (time_thisstep .gt. statistic_equil_params%max_runtime_allowed)) then
      write(*, '(A, ES9.2/)') 'Premature finish: t = ', t
      statistic_equil_params%is_good = .false.
      exit
    end if
    time_laststep = time_thisstep
    runtime_laststep = runtime_thisstep
    !
    if (statistic_equil_params%ISTATE .LT. 0) then
      statistic_equil_params%NERR = statistic_equil_params%NERR + 1
      !write(*, '(A, I3/)') 'Error: ', statistic_equil_params%ISTATE
      statistic_equil_params%ISTATE = 3
    end if
    t_step = t_step * statistic_equil_params%ratio_tstep
    tout = t + t_step
  end do
  !
  if (3*statistic_equil_params%NERR .gt. statistic_equil_params%n_record) then
    write(*, '(/A)') 'Errors occurred:'
    write(*, '(2ES12.4, 2I6, /)') t, statistic_equil_params%t_max, &
      statistic_equil_params%NERR, statistic_equil_params%n_record
    statistic_equil_params%is_good = .false.
  end if
  !
  a_mol_using%f_occupation = a_mol_using%f_occupation / sum(a_mol_using%f_occupation)
  !
  call stat_equili_ode_f(a_mol_using%n_level, &
    t, a_mol_using%f_occupation, &
    statistic_equil_params%RWORK(1:a_mol_using%n_level))
  !
  do i=1, a_mol_using%n_level
    if (a_mol_using%f_occupation(i) .lt. -1D1*statistic_equil_params%ATOL) then
      write(*, '(/A)') 'Negative occupation number occurred.'
      write(*, *) i, a_mol_using%f_occupation(i)
      statistic_equil_params%is_good = .false.
    end if
    tscal = a_mol_using%f_occupation(i)/(abs(statistic_equil_params%RWORK(i))+1D-50)
    if ((a_mol_using%f_occupation(i) .ge. statistic_equil_params%ATOL) .and. &
        (tscal .le. 1D-4*statistic_equil_params%t_max)) then
      write(*, '(A, I4, 3ES12.2)') 'Not yet equilibrium:', i, a_mol_using%f_occupation(i), statistic_equil_params%RWORK(i), tscal
    end if
    if (a_mol_using%f_occupation(i) .lt. 0D0) then
      a_mol_using%f_occupation(i) = 0D0
    end if
  end do
  !
  call calc_energy_derivative(a_mol_using%n_level, a_mol_using%f_occupation, edot)
  write(*,*) 'd/dt \sum \dot{x_i}^2:', edot
end subroutine statistic_equil_solve


subroutine statistic_equil_solve_Newton
  external stat_equili_fcn, stat_equili_jac
  double precision, dimension(a_mol_using%n_level) :: XSCAL
  double precision :: t=0D0, tscal
  integer IERR, i
  integer, dimension(50) :: IOPT
  !
  XSCAL = 1D-10
  !
  IOPT = 0
  statistic_equil_params%IWORK = 0
  statistic_equil_params%RWORK = 0D0
  !
  IOPT(3) = 1 ! JACGEN
  IOPT(11) = 3 ! MPRERR
  IOPT(13) = 1 ! MPRMON
  IOPT(31) = 4 ! NONLIN
  IOPT(32) = 0 ! 1:Broyden
  IOPT(33) = 0 ! 0:Damping
  statistic_equil_params%IWORK(31) = 5000 ! NITER
  statistic_equil_params%RWORK(22) = 1D-20 ! FCMIN
  !
#ifdef USE_EQ_SOLVER_NLEQ2
  call NLEQ2&
#else
  call NLEQ1&
#endif
    (a_mol_using%n_level, &
     stat_equili_fcn, &
     stat_equili_jac, &
     a_mol_using%f_occupation, &
     XSCAL, &
     statistic_equil_params%RTOL, &
     IOPT, &
     IERR, &
     statistic_equil_params%LIW, &
     statistic_equil_params%IWORK, &
     statistic_equil_params%LRW, &
     statistic_equil_params%RWORK)
  if (IERR .eq. 0) then
    write(*, '(A)') 'Iteration converged!'
  else if (IERR .eq. -1) then
    write(*, '(A)') 'More iterations are needed!'
  else if (IERR .eq. 10) then
    write(*, '(A)') 'LIW too small!'
  else
    write(*, '("Error code: ", I4)') IERR
    write(*,*) IOPT
  end if
  !
  a_mol_using%f_occupation = a_mol_using%f_occupation / sum(a_mol_using%f_occupation)
  !
  call stat_equili_ode_f(a_mol_using%n_level, &
    t, a_mol_using%f_occupation, &
    statistic_equil_params%RWORK(1:a_mol_using%n_level))
  !
  do i=1, a_mol_using%n_level
    if (a_mol_using%f_occupation(i) .lt. -1D1*statistic_equil_params%ATOL) then
      write(*, '(/A)') 'Negative occupation number occurred.'
      write(*, *) i, a_mol_using%f_occupation(i)
      statistic_equil_params%is_good = .false.
    end if
    tscal = a_mol_using%f_occupation(i)/(abs(statistic_equil_params%RWORK(i))+1D-50)
        
    if ((a_mol_using%f_occupation(i) .ge. statistic_equil_params%ATOL) .and. &
        (tscal .le. 1D-4*statistic_equil_params%t_max)) then
      write(*, '(A, I4, 3ES12.2)') 'Not yet equilibrium:', i, a_mol_using%f_occupation(i), statistic_equil_params%RWORK(i), tscal
    end if
  end do
end subroutine statistic_equil_solve_Newton


subroutine get_cont_alpha(cont_lut, lam, alp, J)
  type(type_continuum_lut), intent(in) :: cont_lut
  double precision, intent(in) :: lam
  double precision, intent(out) :: alp, J
  integer i, imin, imax, imid, k
  integer, parameter :: ITH = 5
  if (cont_lut%n .le. 0) then
    alp = 0D0
    J = 0D0
    return
  end if
  if (lam .lt. cont_lut%lam(1)) then
    alp = cont_lut%alpha(1)
    J = cont_lut%J(1)
  else if (lam .gt. cont_lut%lam(cont_lut%n)) then
    alp = cont_lut%alpha(cont_lut%n)
    J = cont_lut%J(cont_lut%n)
  else
    imin = 1
    imax = cont_lut%n
    do i=1, cont_lut%n
      if (imin .ge. imax-ITH) then
        do k=imin, imax-1
          if ((cont_lut%lam(k) .le. lam) .and. &
              (cont_lut%lam(k+1) .gt. lam)) then
            alp = cont_lut%alpha(k)
            J = cont_lut%J(k)
            return
          end if
        end do
        exit
      else
        imid = (imin + imax) / 2
        if (cont_lut%lam(imid) .le. lam) then
          imin = imid
        else
          imax = imid
        end if
      end if
    end do
  end if
end subroutine get_cont_alpha


subroutine calc_beta(tau, geotype, beta, dbeta_dtau)
  double precision, intent(in) :: tau
  character(len=*), intent(in) :: geotype
  double precision, intent(out) :: beta, dbeta_dtau
  double precision tmp, t1, t2, A, tt
  double precision, parameter :: const_small_tau = 1D-6
  double precision, parameter :: const_mid_tau = 5D-1
  double precision, parameter :: const_large_tau = 30D0
  double precision, parameter :: LVG_c = 2.34D0 * 0.5D0
  !
  tt = abs(tau)
  select case(geotype)
    case ('spherical', 'Spherical', 'SPHERICAL', 'sphere', 'Sphere', 'SPHERE')
      if (tt .le. const_small_tau) then
        ! Error < const_small_tau
        beta = 1D0
        dbeta_dtau = -3D0/8D0
      else if ((tt .gt. const_small_tau) .and. (tt .le. const_mid_tau)) then
        ! Error < const_mid_tau**7/1e4
        t1 = tau
        beta = 1D0 + &
          t1 * (-3D0/8D0 + t1 * (0.1D0 + t1 * (-1D0/48D0 + t1 * ( &
            1D0/280D0 + t1 * (-1D0/1920D0 + t1 * 1D0/15120D0)))))
        dbeta_dtau = -3D0/8 + &
          t1 * (1D0/5D0 + t1 * (-1D0/16D0 + t1 * (1D0/70D0 + t1 * &
            (-1D0/384D0 + t1 * (1D0/2520D0 - t1 * 1D0/19200D0)))))
        !Series[3/2/t*(1-2/t/t+2*(1/t+1/t/t)*exp(-t)), {t, 0, 0.1}]
        !Series[Diff[3/2/t*(1-2/t/t+2*(1/t+1/t/t)*exp(-t)), t], {t, 0,1}]
      else if ((tt .gt. const_mid_tau) .and. (tt .le. const_large_tau)) then
        ! Exact formula
        ! Error < const_mid_tau**3 * 1e-15 (not sure)
        tmp = exp(-tau)
        t1 = 1D0 / tau
        t2 = t1 * t1
        A = 1D0 - 2D0 * t2 + &
            2D0 * tmp * (t1 + t2)
        beta = 1.5D0 * t1 * A
        dbeta_dtau = -1.5D0 * t2 * A + &
            1.5D0 * t1 * (4D0 * t1 * t2 - &
                          2D0 * tmp * (t1 + 2D0 * t2 + 2D0 * t1 * t2))
      else
        ! Error < exp(-const_large_tau)/tau/tau
        t1 = 1D0 / tau
        t2 = t1 * t1
        beta = 1.5D0 * t1 * (1D0 - 2D0 * t2)
        dbeta_dtau = t2 * (-1.5D0 + 9D0 * t2)
      end if
    case ('lvg', 'LVG', 'Sobolev')
      ! 1980A&A, 91, 68, De Jong et al.: Hydrostatic models of molecular clouds
      ! Their formula (B-7) is modified according to radex to make it (roughly)
      ! continuous.
      if (tt .le. const_small_tau) then
        beta = 1D0
        dbeta_dtau = -LVG_c * 0.5D0
      else if ((tt .gt. const_small_tau) .and. (tt .le. 14D0)) then
        ! In radex the threshold is set to 7.
        ! However, 14 seems to be better for continuity.
        ! This discrepancy with the original De Jong paper seems to stem from
        ! tau is redefined to tau/2 here (as in radex).
        A = LVG_c * tau
        tmp = exp(-A)
        beta = (1D0 - tmp) / A
        dbeta_dtau = LVG_c * (tmp * (1D0 + A) - 1D0) / (A*A)
      else
        A = log(tau * (0.5D0/sqrt(phy_Pi)))
        t2 = sqrt(A)
        t1 = tau * t2
        beta = 1D0 / t1
        dbeta_dtau = -1D0 / (tau * t1) - 0.5D0/(t1*t1*t2)
      end if
    case ('slab', 'Slab', 'SLAB')
      if (tt .le. const_small_tau) then
        beta = 1D0
        dbeta_dtau = -1.5D0
      else
        tmp = exp(-3D0 * tau)
        beta = (1D0 - tmp) / (3D0 * tau)
        dbeta_dtau = tmp * (1D0/tau + 1D0/(3D0*tau*tau)) - 1D0/(3D0*tau*tau)
      end if
    case default
      ! write(*,*) 'Using default geotype: beta=(1-exp(-t))/t'
      if (tt .le. const_small_tau) then
        beta = 1D0
        dbeta_dtau = -0.5D0
      else
        tmp = exp(-tau)
        beta = (1D0 - tmp) / tau
        dbeta_dtau = ((1D0 + tau) * tmp - 1D0) / (tau*tau)
      end if
  end select
end subroutine calc_beta

end module statistic_equilibrium



subroutine stat_equili_fcn(N, X, F, IFAIL)
  integer, intent(in) :: N
  double precision, intent(in), dimension(N) :: X
  double precision, intent(out), dimension(N) :: F
  integer, intent(inout) :: IFAIL
  integer i
  double precision tmp
  call stat_equili_ode_f(N, 0D0, X, F)
  tmp = 0D0
  do i=1, N
    !if (X(i) .lt. 0D0) then
    !  IFAIL = 1
    !  write(*,*) 'X(i)<0: ', 'x(', i, ') = ', x(i)
    !end if
    tmp = tmp + X(i)
  end do
  F(N) = tmp - 1D0
end subroutine stat_equili_fcn



subroutine stat_equili_jac(N, LDJAC, X, DFDX, IFAIL)
  integer, intent(in) :: N, LDJAC
  double precision, intent(in), dimension(N) :: X
  double precision, intent(out), dimension(LDJAC, N) :: DFDX
  integer, intent(inout) :: IFAIL
  integer i
  call stat_equili_ode_jac(N, 0D0, X, 0, 0, DFDX, LDJAC)
  do i=1, N
    !if (X(i) .lt. 0D0) then
    !  IFAIL = 1
    !  exit
    !end if
  end do
  DFDX(LDJAC, 1:N) = 1D0
  !write(*,'(I6,I6,/)') LDJAC,N
end subroutine stat_equili_jac



subroutine stat_equili_ode_f(NEQ, t, y, ydot)
  use statistic_equilibrium
  use phy_const
  implicit none
  integer NEQ
  double precision t, y(NEQ), ydot(NEQ)
  integer i, j, itmp, iup, ilow, iL, iR
  double precision nu, J_ave, rtmp, Tkin, Cul, Clu, TL, TR, deltaE, &
    del_nu, alpha, tau, beta, dbeta_dtau
  double precision lambda, cont_alpha, cont_J_bg, cont_J_in, cont_J_out, tmp
  double precision jnu, knu
  double precision t1
  ydot = 0D0
  Tkin = a_mol_using%Tkin
  do i=1, a_mol_using%rad_data%n_transition
    iup = a_mol_using%rad_data%list(i)%iup
    ilow = a_mol_using%rad_data%list(i)%ilow
    nu = a_mol_using%rad_data%list(i)%freq
    lambda = a_mol_using%rad_data%list(i)%lambda
    del_nu = nu * a_mol_using%dv / phy_SpeedOfLight_CGS * phy_GaussFWHM_c
    call get_cont_alpha(cont_lut_bg, lambda, tmp, cont_J_bg)
    call get_cont_alpha(cont_lut_in, lambda, cont_alpha, cont_J_in)
    call get_cont_alpha(cont_lut_out, lambda, tmp, cont_J_out)
    !
    t1 = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_mol_using%density_mol / del_nu
    jnu = y(iup) *  a_mol_using%rad_data%list(i)%Aul
    knu = y(ilow) * a_mol_using%rad_data%list(i)%Blu - &
          y(iup)  * a_mol_using%rad_data%list(i)%Bul
    alpha = t1 * knu + cont_alpha
    tau = alpha * a_mol_using%length_scale
    !
    call calc_beta(tau, a_mol_using%geotype, beta, dbeta_dtau)
    !
    if ((knu .gt. 1D-30) .or. (knu .lt. -1D-30)) then
      J_ave = jnu / knu
    else
      J_ave = jnu * a_mol_using%length_scale * t1
    end if
    !
    J_ave = J_ave * (1D0 - beta) + (cont_J_bg + cont_J_out) * beta + cont_J_in
    !
    a_mol_using%rad_data%list(i)%tau = tau
    a_mol_using%rad_data%list(i)%beta = beta
    a_mol_using%rad_data%list(i)%J_ave = J_ave
    a_mol_using%rad_data%list(i)%J_cont_bg = cont_J_bg
    !
    rtmp = a_mol_using%rad_data%list(i)%Aul * y(iup) + &
           a_mol_using%rad_data%list(i)%Bul * J_ave * y(iup) - &
           a_mol_using%rad_data%list(i)%Blu * J_ave * y(ilow)
    ydot(iup) = ydot(iup)   - rtmp
    ydot(ilow) = ydot(ilow) + rtmp
    !write(*, *) rtmp, beta, alpha, t1, knu, y(ilow), y(iup), iup, ilow, i
    !if (isnan(knu)) stop
  end do
  do i=1, a_mol_using%colli_data%n_partner
    ! Find the T interval
    itmp = a_mol_using%colli_data%list(i)%n_T
    if (Tkin .le. a_mol_using%colli_data%list(i)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. a_mol_using%colli_data%list(i)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do j=2, a_mol_using%colli_data%list(i)%n_T
        if ((Tkin .ge. a_mol_using%colli_data%list(i)%T_coll(j-1)) .and. &
            (Tkin .le. a_mol_using%colli_data%list(i)%T_coll(j))) then
          iL = j-1
          iR = j
          exit
        end if
      end do
    end if
    do j=1, a_mol_using%colli_data%list(i)%n_transition
      iup = a_mol_using%colli_data%list(i)%iup(j)
      ilow = a_mol_using%colli_data%list(i)%ilow(j)
      deltaE = a_mol_using%level_list(iup)%energy - a_mol_using%level_list(ilow)%energy
      if (iL .eq. iR) then
        Cul = a_mol_using%colli_data%list(i)%Cul(iL, j)
      else
        TL = a_mol_using%colli_data%list(i)%T_coll(iL)
        TR = a_mol_using%colli_data%list(i)%T_coll(iR)
        Cul = (a_mol_using%colli_data%list(i)%Cul(iL, j) * (TR - Tkin) + &
                a_mol_using%colli_data%list(i)%Cul(iR, j) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             a_mol_using%level_list(iup)%weight / &
             a_mol_using%level_list(ilow)%weight
      rtmp = (Cul * y(iup) - Clu * y(ilow)) * a_mol_using%colli_data%list(i)%dens_partner
      ydot(iup) = ydot(iup)   - rtmp
      ydot(ilow) = ydot(ilow) + rtmp
    end do
  end do
  !do i=1, NEQ
  !  write(*, '(ES12.4, I4, ES12.4, ES12.4)') t, i, y(i), ydot(i)
  !end do
end subroutine stat_equili_ode_f




subroutine stat_equili_ode_jac(NEQ, t, y, ML, MU, PD, NROWPD)
  ! PD(i,j) = \partial f_i / \partial x_j
  use statistic_equilibrium
  use phy_const
  implicit none
  double precision t
  integer ML, MU, NROWPD
  double precision, dimension(NEQ) :: y
  double precision, dimension(NROWPD, *) :: PD
  integer NEQ
  integer i, j, itmp, iup, ilow, iL, iR
  double precision nu, J_ave, &
        Tkin, Cul, Clu, TL, TR, deltaE, del_nu, alpha, tau, beta
  double precision lambda, cont_alpha, cont_J_bg, cont_J_in, cont_J_out, tmp
  double precision S, dbeta_dtau, dtau_dy_up, dtau_dy_low, &
    dJ_ave_dy_up, dJ_ave_dy_low, drtmp_dy_up, drtmp_dy_low, &
    dS_dy_up, dS_dy_low
  double precision jnu, knu
  double precision t1
  do j=1,NEQ
    PD(:,j) = 0D0
  end do
  Tkin = a_mol_using%Tkin
  do i=1, a_mol_using%rad_data%n_transition
    iup = a_mol_using%rad_data%list(i)%iup
    ilow = a_mol_using%rad_data%list(i)%ilow
    nu = a_mol_using%rad_data%list(i)%freq
    lambda = a_mol_using%rad_data%list(i)%lambda
    del_nu = nu * a_mol_using%dv / phy_SpeedOfLight_CGS * phy_GaussFWHM_c
    call get_cont_alpha(cont_lut_bg, lambda, tmp, cont_J_bg)
    call get_cont_alpha(cont_lut_in, lambda, cont_alpha, cont_J_in)
    call get_cont_alpha(cont_lut_out, lambda, tmp, cont_J_out)
    !
    t1 = phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_mol_using%density_mol / del_nu
    jnu = y(iup) *  a_mol_using%rad_data%list(i)%Aul
    knu = y(ilow) * a_mol_using%rad_data%list(i)%Blu - &
          y(iup)  * a_mol_using%rad_data%list(i)%Bul
    alpha = t1 * knu + cont_alpha
    tau = alpha * a_mol_using%length_scale
    !
    call calc_beta(tau, a_mol_using%geotype, beta, dbeta_dtau)
    !
    if ((knu .gt. 1D-30) .or. (knu .lt. -1D-30)) then
      S = jnu / knu
      dS_dy_up = (a_mol_using%rad_data%list(i)%Aul + &
                  S * a_mol_using%rad_data%list(i)%Bul) / knu
      dS_dy_low = -S * a_mol_using%rad_data%list(i)%Blu / knu
    else
      S = jnu * a_mol_using%length_scale * t1
      dS_dy_up = a_mol_using%rad_data%list(i)%Aul * a_mol_using%length_scale * t1
      dS_dy_low = 0D0
    end if
    !
    J_ave = S * (1D0 - beta) + (cont_J_bg + cont_J_out) * beta + cont_J_in
    !
    dtau_dy_up = a_mol_using%length_scale * &
                 phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_mol_using%density_mol * &
                 (-a_mol_using%rad_data%list(i)%Bul) / del_nu
    dtau_dy_low = a_mol_using%length_scale * &
                  phy_hPlanck_CGS * nu / (4D0*phy_Pi) * a_mol_using%density_mol * &
                  (a_mol_using%rad_data%list(i)%Blu) / del_nu
    !
    dJ_ave_dy_up  = -S * dbeta_dtau * dtau_dy_up + dS_dy_up * (1D0 - beta)
    dJ_ave_dy_low = -S * dbeta_dtau * dtau_dy_low + dS_dy_low * (1D0 - beta)
    !
    drtmp_dy_up = a_mol_using%rad_data%list(i)%Aul + &
             a_mol_using%rad_data%list(i)%Bul * J_ave + &
             (a_mol_using%rad_data%list(i)%Bul * y(iup) - &
              a_mol_using%rad_data%list(i)%Blu * y(ilow)) * dJ_ave_dy_up
    drtmp_dy_low = -a_mol_using%rad_data%list(i)%Blu * J_ave + &
             (a_mol_using%rad_data%list(i)%Bul * y(iup) - &
              a_mol_using%rad_data%list(i)%Blu * y(ilow)) * dJ_ave_dy_low
    PD(iup,  iup)  = PD(iup,  iup)  - drtmp_dy_up
    PD(ilow, iup)  = PD(ilow, iup)  + drtmp_dy_up
    PD(iup,  ilow) = PD(iup,  ilow) - drtmp_dy_low
    PD(ilow, ilow) = PD(ilow, ilow) + drtmp_dy_low
  end do
  do i=1, a_mol_using%colli_data%n_partner
    ! Find the T interval
    itmp = a_mol_using%colli_data%list(i)%n_T
    if (Tkin .le. a_mol_using%colli_data%list(i)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. a_mol_using%colli_data%list(i)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do j=2, a_mol_using%colli_data%list(i)%n_T
        if ((Tkin .ge. a_mol_using%colli_data%list(i)%T_coll(j-1)) .and. &
            (Tkin .le. a_mol_using%colli_data%list(i)%T_coll(j))) then
          iL = j-1
          iR = j
          exit
        end if
      end do
    end if
    do j=1, a_mol_using%colli_data%list(i)%n_transition
      iup = a_mol_using%colli_data%list(i)%iup(j)
      ilow = a_mol_using%colli_data%list(i)%ilow(j)
      deltaE = a_mol_using%level_list(iup)%energy - &
               a_mol_using%level_list(ilow)%energy
      if (iL .eq. iR) then
        Cul = a_mol_using%colli_data%list(i)%Cul(iL, j)
      else
        TL = a_mol_using%colli_data%list(i)%T_coll(iL)
        TR = a_mol_using%colli_data%list(i)%T_coll(iR)
        Cul = (a_mol_using%colli_data%list(i)%Cul(iL, j) * (TR - Tkin) + &
               a_mol_using%colli_data%list(i)%Cul(iR, j) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             a_mol_using%level_list(iup)%weight / &
             a_mol_using%level_list(ilow)%weight
      drtmp_dy_up  = Cul  * a_mol_using%colli_data%list(i)%dens_partner
      drtmp_dy_low = -Clu * a_mol_using%colli_data%list(i)%dens_partner
      PD(iup,  iup)  = PD(iup,  iup) - drtmp_dy_up
      PD(ilow, iup)  = PD(ilow, iup) + drtmp_dy_up
      PD(iup,  ilow) = PD(iup,  ilow) - drtmp_dy_low
      PD(ilow, ilow) = PD(ilow, ilow) + drtmp_dy_low
      !write(*,*) iup, ilow, y(iup), y(ilow), drtmp_dy_up, drtmp_dy_low
      !if (isnan(drtmp_dy_up) .or. isnan(drtmp_dy_low)) stop
    end do
  end do
end subroutine stat_equili_ode_jac



subroutine calc_critical_density_for_one_transition(iTran, tau)
  use statistic_equilibrium
  use phy_const
  implicit none
  integer, intent(in) :: iTran
  double precision, intent(in) :: tau
  !
  integer i, icp, ict, itmp, iup_rad, ilow_rad, iup_col, ilow_col, iL, iR
  double precision nu, J_ave, Tkin, Cul, Clu, TL, TR, deltaE, &
    del_nu, alpha, beta, dbeta_dtau
  double precision lambda, cont_alpha, cont_J_bg, cont_J_in, cont_J_out, tmp
  double precision jnu, knu
  double precision radiative_transfer_rate_from_iup_to_ilow
  double precision collisional_transfer_rate_from_iup_to_others
  Tkin = a_mol_using%Tkin
  iup_rad = a_mol_using%rad_data%list(iTran)%iup
  ilow_rad = a_mol_using%rad_data%list(iTran)%ilow
  nu = a_mol_using%rad_data%list(iTran)%freq
  lambda = a_mol_using%rad_data%list(iTran)%lambda
  !
  call calc_beta(tau, a_mol_using%geotype, beta, dbeta_dtau)
  !
  radiative_transfer_rate_from_iup_to_ilow = a_mol_using%rad_data%list(iTran)%Aul * beta
  !
  do icp=1, a_mol_using%colli_data%n_partner
    ! Find the T interval for interpolation
    itmp = a_mol_using%colli_data%list(icp)%n_T
    if (Tkin .le. a_mol_using%colli_data%list(icp)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. a_mol_using%colli_data%list(icp)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do i=2, a_mol_using%colli_data%list(icp)%n_T
        if ((Tkin .ge. a_mol_using%colli_data%list(icp)%T_coll(i-1)) .and. &
            (Tkin .le. a_mol_using%colli_data%list(icp)%T_coll(i))) then
          iL = i-1
          iR = i
          exit
        end if
      end do
    end if
    !
    collisional_transfer_rate_from_iup_to_others = 0D0
    do ict=1, a_mol_using%colli_data%list(icp)%n_transition
      iup_col = a_mol_using%colli_data%list(icp)%iup(ict)
      ilow_col = a_mol_using%colli_data%list(icp)%ilow(ict)
      if ((iup_col .ne. iup_rad) .and. (ilow_col .ne. iup_rad)) then
        cycle
      end if
      deltaE = a_mol_using%level_list(iup_col)%energy - a_mol_using%level_list(ilow_col)%energy
      ! Do the interpolation
      if (iL .eq. iR) then
        Cul = a_mol_using%colli_data%list(icp)%Cul(iL, ict)
      else
        TL = a_mol_using%colli_data%list(icp)%T_coll(iL)
        TR = a_mol_using%colli_data%list(icp)%T_coll(iR)
        Cul = (a_mol_using%colli_data%list(icp)%Cul(iL, ict) * (TR - Tkin) + &
                a_mol_using%colli_data%list(icp)%Cul(iR, ict) * (Tkin - TL)) / (TR - TL)
      end if
      !
      Clu = Cul * exp(-deltaE/Tkin) * &
             a_mol_using%level_list(iup_col)%weight / &
             a_mol_using%level_list(ilow_col)%weight
      if (iup_col .eq. iup_rad) then
        collisional_transfer_rate_from_iup_to_others = &
          collisional_transfer_rate_from_iup_to_others + Cul
      else if (ilow_col .eq. iup_rad) then
        collisional_transfer_rate_from_iup_to_others = &
          collisional_transfer_rate_from_iup_to_others + Clu
      else
        write(*,*) "Something is wrong in calc_critical_density_for_one_transition!"
      end if
    end do
    a_mol_using%rad_data%list(iTran)%critical_densities(icp) = &
      radiative_transfer_rate_from_iup_to_ilow / &
      (collisional_transfer_rate_from_iup_to_others + 1d-200)
  end do
end subroutine calc_critical_density_for_one_transition





subroutine calc_critical_density_old_def_for_one_transition(iTran, tau)
  use statistic_equilibrium
  use phy_const
  implicit none
  integer, intent(in) :: iTran
  double precision, intent(in) :: tau
  integer i, icp, ict, itmp, iup_rad, ilow_rad, iup_col, ilow_col, iL, iR
  double precision nu, J_ave, Tkin, Cul, Clu, TL, TR, deltaE, &
    del_nu, alpha, beta, dbeta_dtau
  double precision lambda, cont_alpha, cont_J_bg, cont_J_in, cont_J_out, tmp
  double precision jnu, knu
  double precision radiative_transfer_rate_from_iup_to_ilow
  double precision collisional_transfer_rate_from_iup_to_others
  Tkin = a_mol_using%Tkin
  iup_rad = a_mol_using%rad_data%list(iTran)%iup
  ilow_rad = a_mol_using%rad_data%list(iTran)%ilow
  nu = a_mol_using%rad_data%list(iTran)%freq
  lambda = a_mol_using%rad_data%list(iTran)%lambda
  !
  call calc_beta(tau, a_mol_using%geotype, beta, dbeta_dtau)
  !
  radiative_transfer_rate_from_iup_to_ilow = a_mol_using%rad_data%list(iTran)%Aul * beta
  !
  do icp=1, a_mol_using%colli_data%n_partner
    ! Find the T interval for interpolation
    itmp = a_mol_using%colli_data%list(icp)%n_T
    if (Tkin .le. a_mol_using%colli_data%list(icp)%T_coll(1)) then
      iL = 1
      iR = 1
    else if (Tkin .ge. a_mol_using%colli_data%list(icp)%T_coll(itmp)) then
      iL = itmp
      iR = itmp
    else
      do i=2, a_mol_using%colli_data%list(icp)%n_T
        if ((Tkin .ge. a_mol_using%colli_data%list(icp)%T_coll(i-1)) .and. &
            (Tkin .le. a_mol_using%colli_data%list(icp)%T_coll(i))) then
          iL = i-1
          iR = i
          exit
        end if
      end do
    end if
    collisional_transfer_rate_from_iup_to_others = 0D0
    do ict=1, a_mol_using%colli_data%list(icp)%n_transition
      iup_col = a_mol_using%colli_data%list(icp)%iup(ict)
      ilow_col = a_mol_using%colli_data%list(icp)%ilow(ict)
      if ((iup_col .ne. iup_rad) .or. (ilow_col .ne. ilow_rad)) then
        ! This makes the difference with calc_critical_density_f
        cycle
      end if
      deltaE = a_mol_using%level_list(iup_col)%energy - a_mol_using%level_list(ilow_col)%energy
      if (iL .eq. iR) then
        Cul = a_mol_using%colli_data%list(icp)%Cul(iL, ict)
      else
        TL = a_mol_using%colli_data%list(icp)%T_coll(iL)
        TR = a_mol_using%colli_data%list(icp)%T_coll(iR)
        Cul = (a_mol_using%colli_data%list(icp)%Cul(iL, ict) * (TR - Tkin) + &
                a_mol_using%colli_data%list(icp)%Cul(iR, ict) * (Tkin - TL)) / (TR - TL)
      end if
      Clu = Cul * exp(-deltaE/Tkin) * &
             a_mol_using%level_list(iup_col)%weight / &
             a_mol_using%level_list(ilow_col)%weight
      if (iup_col .eq. iup_rad) then
        collisional_transfer_rate_from_iup_to_others = &
          collisional_transfer_rate_from_iup_to_others + Cul
      else if (ilow_col .eq. iup_rad) then
        collisional_transfer_rate_from_iup_to_others = &
          collisional_transfer_rate_from_iup_to_others + Clu
      else
        write(*,*) "Something is wrong in calc_critical_density_old_def_for_one_transition!"
      end if
    end do
    a_mol_using%rad_data%list(iTran)%critical_densities(icp) = &
      radiative_transfer_rate_from_iup_to_ilow / &
      (collisional_transfer_rate_from_iup_to_others + 1d-200)
  end do
end subroutine calc_critical_density_old_def_for_one_transition



subroutine calc_critical_density_f(tau)
  use statistic_equilibrium
  use phy_const
  implicit none
  double precision, intent(in) :: tau
  integer irt
  do irt=1, a_mol_using%rad_data%n_transition
    call calc_critical_density_for_one_transition(irt, tau)
  end do
end subroutine calc_critical_density_f



subroutine calc_critical_density_old_def_f(tau)
  use statistic_equilibrium
  use phy_const
  implicit none
  double precision, intent(in) :: tau
  integer irt
  do irt=1, a_mol_using%rad_data%n_transition
    call calc_critical_density_old_def_for_one_transition(irt, tau)
  end do
end subroutine calc_critical_density_old_def_f



subroutine calc_energy_derivative(n, x, edot)
  external stat_equili_ode_f, stat_equili_ode_jac
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x
  double precision, intent(out) :: edot
  double precision :: t = 0D0
  double precision, dimension(n) :: f
  double precision, dimension(n,n) :: jac
  integer i,j
  call stat_equili_ode_f(n, t, x, f)
  call stat_equili_ode_jac(n, t, x, 0, 0, jac, n)
  edot = dot_product(f,matmul(jac,f))
end subroutine calc_energy_derivative
