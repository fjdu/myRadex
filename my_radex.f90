! This file is part of our code to solve the statistical equilibrium problem of
! the energy level population of a uniform volume of gas under the effect of
! background radiation and collisional excitation.
!
! The code works similar to radex
! (see http://home.strw.leidenuniv.nl/~moldata/radex.html)
!
! One difference from radex is that we solve the statistical equilibrium
! equation using an ode solver (ODEPACK, see netlib.org/odepack/).
!
! Written by Fujun Du, fujun.du@gmail.com, fdu@umich.edu
!
! 2014-01-02 Thu 02:43:39
! 2014-01-02 Thu 22:58:10
!


module my_radex

use phy_const
use trivials
use statistic_equilibrium

implicit none

! Max number of elements for vectorized config params
integer, parameter :: ndim_cfg_vec = 100
integer, parameter :: ndim_Tbg = 8

type :: type_rdxx_cfg
  character(len=128) :: dir_transition_rates = ''
  character(len=128) :: dir_save = ''
  character(len=128) :: filename_molecule = ''
  character(len=128) :: filename_save = ''
  logical :: recalculateFreqWithEupElow = .false.
  logical :: iLevel_subtract_one = .false.
  logical :: verbose = .true.
  !
  double precision :: freqmin=0D0, freqmax=1D99
  !
  double precision :: max_evol_time = 1D16
  double precision :: max_code_run_time = 10.0
  double precision :: rtol = 1D-4, atol = 1D-20
  !
  character(len=16)  :: geotype = ''
  !
  integer :: nTbg = 1
  double precision, dimension(ndim_Tbg) :: Tbg = 0D0
  double precision, dimension(ndim_Tbg) :: Tbgscaling = 1D0
  logical :: provide_bgfile = .false.
  character(len=128) :: dir_bg = ''
  character(len=128) :: filename_bg = ''
  !
  integer :: nTin = 0
  double precision, dimension(ndim_Tbg) :: Tin = 0D0
  double precision, dimension(ndim_Tbg) :: Tinscaling = 1D0
  logical :: provide_infile = .false.
  character(len=128) :: dir_in = ''
  character(len=128) :: filename_in = ''
  !
  integer :: nTout = 0
  double precision, dimension(ndim_Tbg) :: Tout = 0D0
  double precision, dimension(ndim_Tbg) :: Toutscaling = 1D0
  logical :: provide_outfile = .false.
  character(len=128) :: dir_out = ''
  character(len=128) :: filename_out = ''
  !
  character(len=32) :: solve_method = ''
  character(len=32) :: f_occupation_init_method = ''
  !
  logical :: provideLength = .false.
  double precision length_scale
  !
  double precision :: beam_FWHM_in_arcsec = 0.4D0
  !
  double precision, dimension(ndim_cfg_vec) :: &
    Tkin, dv, Ncol_x, n_x, &
    n_H2, n_HI, n_oH2, n_pH2, n_Hplus, n_E, n_He   
  integer :: nTkin=1, ndv=1, nn_x=1, nNcol_x=1, ndens=1 ! Vector sizes
  integer iTkin, idv, in_x, iNcol_x, idens ! Loop indices
  double precision ::  opH2_ratio = 3D0
  logical :: opH2eq3 = .false.
  integer fU
  integer :: collisioPartnerCrit = 1
end type type_rdxx_cfg

type(type_rdxx_cfg) :: rdxx_cfg

namelist /rdxx_configure/ rdxx_cfg


contains


subroutine do_my_radex(do_init)
  logical, intent(in), optional :: do_init
  logical do_ini
  integer i, itot, ntot, ilowSav, iupSav
  integer iTkin, idv, in_x, iNcol_x, idens
  double precision fup, flow, gup, glow, Tex, Tr, flux_CGS, flux_K_km_s, flux_Jy
  double precision Inu_t, tau, t1, t2
  double precision beam_area
  double precision critical_density, critical_density_old
  integer flag_good
  character(len=9) beam_FWHM_str
  !
  if (present(do_init)) then
    do_ini = do_init
  else
    do_ini = .true.
  end if
  !
  write(*, '(A)') 'Code running...'
  if (.not. rdxx_cfg%verbose) then
    write(*, '(A)') 'Runtime message disabled.'
  end if
  !
  if (do_ini) then
    ! Load the molecular data, etc.
    call my_radex_prepare
  end if
  !
  write(beam_FWHM_str, '("(", F4.1, """)")') rdxx_cfg%beam_FWHM_in_arcsec
  !
  call openFileSequentialWrite(rdxx_cfg%fU, &
    combine_dir_filename(rdxx_cfg%dir_save, &
    rdxx_cfg%filename_save), 999, 1)
  write(rdxx_cfg%fU, '(2A5, A12, 2A15, 10A12, 2A7, &
            &5A12, A2)') &
    '! iup', 'ilow', 'Eup', 'freq', 'lam', 'Tex', 'tau', 'Tr', &
    'fup', 'flow', 'flux_Kkms', 'flux_int', 'flux_Jy', 'beta', &
    'Jnu', 'gup', 'glow', 'Aul', 'Bul', 'Blu', 'n_crit', 'n_crit_old', 'q'
  write(rdxx_cfg%fU, '(2A5, A12, 2A15, 10A12, 2A7, &
            &5A12, A2)') &
    '!    ', '  ', 'K', 'Hz', 'micron', 'K', '', 'K', &
    '   ', '    ', 'K km/s', 'erg/cm2/s', 'Jy'//trim(beam_FWHM_str), '    ', &
    '...', '   ', '    ', 's-1', '...', '...', 'cm-3', 'cm-3', ''
  !
  write(*,*) 'Critical density will be included for collisional partner: ',  rdxx_cfg%collisioPartnerCrit
  !
  itot = 0
  ntot = rdxx_cfg%nTkin * rdxx_cfg%ndv * rdxx_cfg%nn_x * &
         rdxx_cfg%nNcol_x * rdxx_cfg%ndens
  !
  ! Big loop starts here
  !
  do iTkin=1, rdxx_cfg%nTkin
  do idv=1, rdxx_cfg%ndv
  do in_x=1, rdxx_cfg%nn_x
  do iNcol_x=1, rdxx_cfg%nNcol_x
  do idens=1, rdxx_cfg%ndens
    !
    itot = itot + 1
    !
    rdxx_cfg%iTkin   = iTkin
    rdxx_cfg%idv     = idv
    rdxx_cfg%in_x    = in_x
    rdxx_cfg%iNcol_x = iNcol_x
    rdxx_cfg%idens   = idens
    !
    write(rdxx_cfg%fU, '("!", I6, "/", I6, A, 5I5, " / ", 5I5, 2X, A)') &
      itot, ntot, ': ', &
      rdxx_cfg%iTkin, rdxx_cfg%idv, rdxx_cfg%in_x, &
      rdxx_cfg%iNcol_x, rdxx_cfg%idens, &
      rdxx_cfg%nTkin, rdxx_cfg%ndv, rdxx_cfg%nn_x, &
      rdxx_cfg%nNcol_x, rdxx_cfg%ndens, &
      'Loop order: Tkin, dv, n_x, Ncol_x, dens'
    write(rdxx_cfg%fU, '(3("!", ES12.4E2, " =", A8/), ("!", ES12.4E2, " =", A8))') &
      rdxx_cfg%Tkin(rdxx_cfg%iTkin), 'Tkin', &
      rdxx_cfg%dv(rdxx_cfg%idv), 'dv', &
      rdxx_cfg%n_x(rdxx_cfg%in_x), 'n_x', &
      rdxx_cfg%Ncol_x(rdxx_cfg%iNcol_x), 'Ncol_x'
    !
    if (rdxx_cfg%verbose) then
      write(*, '("!", I6, "/", I6, A, 5I5, " / ", 5I5, 2X, A)') &
        itot, ntot, ': ', &
        rdxx_cfg%iTkin, rdxx_cfg%idv, rdxx_cfg%in_x, &
        rdxx_cfg%iNcol_x, rdxx_cfg%idens, &
        rdxx_cfg%nTkin, rdxx_cfg%ndv, rdxx_cfg%nn_x, &
        rdxx_cfg%nNcol_x, rdxx_cfg%ndens, &
        'Loop order: Tkin, dv, n_x, Ncol_x, dens'
      write(*, '(3("!", ES12.4E2, " =", A8/), ("!", ES12.4E2, " =", A8))') &
        rdxx_cfg%Tkin(rdxx_cfg%iTkin), 'Tkin', &
        rdxx_cfg%dv(rdxx_cfg%idv), 'dv', &
        rdxx_cfg%n_x(rdxx_cfg%in_x), 'n_x', &
        rdxx_cfg%Ncol_x(rdxx_cfg%iNcol_x), 'Ncol_x'
    end if
    !
    call my_radex_prepare_molecule
    select case(rdxx_cfg%solve_method)
      case ('ODE', 'ode')
        call statistic_equil_solve
      case ('Newton', 'NEWTON')
        call statistic_equil_solve_Newton
      case default
        write(*,*) 'Unknown method: "', trim(rdxx_cfg%solve_method), '"'
    end select
    !
    if (statistic_equil_params%is_good) then
      flag_good = 1
    else
      flag_good = 0
    end if
    !
    do i=1, a_mol_using%rad_data%n_transition
      associate(r => a_mol_using%rad_data%list(i))
        if ((r%freq .lt. rdxx_cfg%freqmin) .or. &
            (r%freq .gt. rdxx_cfg%freqmax)) then
          cycle
        end if
        !
        fup  = a_mol_using%f_occupation(r%iup)
        flow = a_mol_using%f_occupation(r%ilow)
        gup  = a_mol_using%level_list(r%iup)%weight
        glow = a_mol_using%level_list(r%ilow)%weight
        Tex = -(r%Eup - r%Elow) / log(fup*glow / (flow*gup))
        !
        tau = r%tau
        t1 = exp(-tau)
        if (abs(tau) .lt. 1D-6) then
          t2 = tau
        else
          t2 = 1D0 - t1
        end if
        Inu_t = planck_B_nu(Tex, r%freq) * t2 + t1 * r%J_cont_bg
        !
        Tr = (Inu_t - r%J_cont_bg) * phy_SpeedOfLight_CGS**2 / &
          (2D0 * r%freq**2 * phy_kBoltzmann_CGS)
        flux_K_km_s = Tr * a_mol_using%dv / 1D5 * phy_GaussFWHM_c
        flux_CGS = (Inu_t - r%J_cont_bg) * &
          a_mol_using%dv * r%freq / phy_SpeedOfLight_CGS * (4D0 * phy_Pi * phy_GaussFWHM_c)
        beam_area = FWHM_to_area(rdxx_cfg%beam_FWHM_in_arcsec)
        flux_Jy = (Inu_t - r%J_cont_bg) * beam_area / phy_jansky2CGS
        ! flux_CGS = \int I_\nu d\nu d\Omega
        !
        ilowSav = r%ilow
        iupSav = r%iup
        if (rdxx_cfg%iLevel_subtract_one) then
          ! Some people may want the level numbers be subtracted by one.
          ilowSav = r%ilow - 1
          iupSav = r%iup - 1
        end if
        !
        call calc_critical_density_for_one_transition(i, tau)
        critical_density = a_mol_using%rad_data%list(i)%critical_densities(rdxx_cfg%collisioPartnerCrit)
        call calc_critical_density_old_def_for_one_transition(i, tau)
        critical_density_old = a_mol_using%rad_data%list(i)%critical_densities(rdxx_cfg%collisioPartnerCrit)
        !
        write(rdxx_cfg%fU, '(2I5, F12.4, 2ES15.7E2, 10ES12.3E3, 2F7.1, &
                  &5ES12.3E3, I2)') &
          iupSav, ilowSav, r%Eup, r%freq, r%lambda, Tex, r%tau, Tr, &
          fup, flow, flux_K_km_s, flux_CGS, flux_Jy, r%beta, &
          r%J_ave, gup, glow, r%Aul, r%Bul, r%Blu, critical_density, critical_density_old, flag_good
      end associate
    end do
    flush(rdxx_cfg%fU)
  end do
  end do
  end do
  end do
  end do
  !
  ! Big loop ends here
  !
  close(rdxx_cfg%fU)
  nullify(a_mol_using)
end subroutine do_my_radex


pure function FWHM_to_area(FWHM) result(res)
  ! Given f(r) = exp(-r^2/a):
  !     Area = Int[exp(-r^2/a) 2pi rdr] = a pi
  ! exp(-r^2/a) = 1/2
  !   r = sqrt(a ln2)
  !   FWHM = 2sqrt(a ln2)
  ! Hence
  ! Area = pi FWHM^2 / (4ln2)
  double precision, intent(in) :: FWHM
  double precision :: res
  res = phy_Pi / log(16D0) * (FWHM / 3.6D3 * phy_Deg2Rad)**2
end function FWHM_to_area


subroutine my_radex_prepare_molecule
  integer i
  !
  if (rdxx_cfg%verbose) then
    write(*,*) 'Using geotype:"', rdxx_cfg%geotype, '"'
  end if
  a_mol_using%geotype = rdxx_cfg%geotype ! Geometric type
  !
  a_mol_using%Tkin = rdxx_cfg%Tkin(rdxx_cfg%iTkin) ! K
  a_mol_using%dv = rdxx_cfg%dv(rdxx_cfg%idv) ! cm s-1
  !
  ! When the continuum opacity is zero, the density of the molecule being
  ! studied does not really enter the calculation.
  !
  a_mol_using%density_mol = rdxx_cfg%n_x(rdxx_cfg%in_x)
  if (a_mol_using%density_mol .le. 1D-20) then
    ! If not set, set it to a non-harmful value.
    a_mol_using%density_mol = 1D0
  end if
  !
  if (rdxx_cfg%provideLength) then
    a_mol_using%length_scale = rdxx_cfg%length_scale ! cm
  else
    a_mol_using%length_scale = rdxx_cfg%Ncol_x(rdxx_cfg%iNcol_x) / &
                               a_mol_using%density_mol
  end if
  !
  ! Set the initial occupation
  select case (rdxx_cfg%f_occupation_init_method)
    case ('Boltzmann', 'BOLTZMANN', 'boltzmann')
      a_mol_using%f_occupation = a_mol_using%level_list%weight * &
        exp(-a_mol_using%level_list%energy / a_mol_using%Tkin)
    case default
      call random_number(a_mol_using%f_occupation)
  end select
  ! Normalize
  a_mol_using%f_occupation = a_mol_using%f_occupation / &
                             sum(a_mol_using%f_occupation)
  !
  ! Ortho/para ratio of H2.
  if (rdxx_cfg%opH2eq3) then
    rdxx_cfg%opH2_ratio = 3D0
  else
    rdxx_cfg%opH2_ratio = calc_ortho_para_H2_ratio(a_mol_using%Tkin)
  end if
  !
  ! Set the density of the collisional partners
  do i=1, a_mol_using%colli_data%n_partner
    select case (a_mol_using%colli_data%list(i)%name_partner)
      case ('H2', 'h2')
      !
        a_mol_using%colli_data%list(i)%dens_partner = &
          rdxx_cfg%n_H2(rdxx_cfg%idens)
      !
      case ('o-H2', 'oH2', 'o_H2')
      !
        if (rdxx_cfg%n_oH2(rdxx_cfg%idens) .le. 1D-20) then
          a_mol_using%colli_data%list(i)%dens_partner = &
            rdxx_cfg%n_H2(rdxx_cfg%idens) * &
            rdxx_cfg%opH2_ratio / (1D0 + rdxx_cfg%opH2_ratio)
        else
          a_mol_using%colli_data%list(i)%dens_partner = &
            rdxx_cfg%n_oH2(rdxx_cfg%idens)
        end if
      !
      case ('p-H2', 'pH2', 'p_H2')
      !
        if (rdxx_cfg%n_pH2(rdxx_cfg%idens) .le. 1D-20) then
          a_mol_using%colli_data%list(i)%dens_partner = &
            rdxx_cfg%n_H2(rdxx_cfg%idens) * &
            1D0 / (1D0 + rdxx_cfg%opH2_ratio)
        else
          a_mol_using%colli_data%list(i)%dens_partner = &
            rdxx_cfg%n_pH2(rdxx_cfg%idens)
        end if
      !
      case ('H', 'h')
      !
        a_mol_using%colli_data%list(i)%dens_partner = &
          rdxx_cfg%n_HI(rdxx_cfg%idens)
      !
      case ('H+', 'h+')
      !
        a_mol_using%colli_data%list(i)%dens_partner = &
          rdxx_cfg%n_Hplus(rdxx_cfg%idens)
      !
      case ('E', 'e', 'E-', 'e-')
      !
        a_mol_using%colli_data%list(i)%dens_partner = &
          rdxx_cfg%n_E(rdxx_cfg%idens)
      !
      case ('He', 'HE')
      !
        a_mol_using%colli_data%list(i)%dens_partner = &
          rdxx_cfg%n_He(rdxx_cfg%idens)
      !
      case default
      !
        write(*, '(/A, A)') 'Unknown collisional partner: ', &
          a_mol_using%colli_data%list(i)%name_partner
        a_mol_using%colli_data%list(i)%dens_partner = 0D0
      !
    end select
    !
    write(rdxx_cfg%fU, '("!", ES12.4E2, " =", A8)') &
      a_mol_using%colli_data%list(i)%dens_partner, &
      trim(a_mol_using%colli_data%list(i)%name_partner)
    if (rdxx_cfg%verbose) then
      write(*, '("!", ES12.4E2, " =", A8)') &
        a_mol_using%colli_data%list(i)%dens_partner, &
        trim(a_mol_using%colli_data%list(i)%name_partner)
    end if
  end do
  !
end subroutine my_radex_prepare_molecule



subroutine my_radex_prepare(loadMoleculeData, makeLUT)
  logical, intent(in), optional :: loadMoleculeData, makeLUT
  logical loadMoleculeData_, makeLUT_
  double precision lam_min, lam_max
  integer, parameter :: n_cont_lam = 200
  !
  loadMoleculeData_ = .true.
  if (present(loadMoleculeData)) then
    loadMoleculeData_ = loadMoleculeData
  end if
  makeLUT_ = .true.
  if (present(makeLUT)) then
    makeLUT_ = makeLUT
  end if
  !
  ! Load the molecular data
  if (loadMoleculeData_) then
    if (associated(a_mol_using)) then
      call deallocate_a_mol_using()
    end if
    allocate(a_mol_using)
    !
    call load_moldata_LAMBDA(&
      combine_dir_filename(rdxx_cfg%dir_transition_rates, &
      rdxx_cfg%filename_molecule), rdxx_cfg%recalculateFreqWithEupElow)
  end if
  !
  statistic_equil_params%max_runtime_allowed = rdxx_cfg%max_code_run_time
  statistic_equil_params%rtol = rdxx_cfg%rtol
  statistic_equil_params%atol = rdxx_cfg%atol
  !
  ! Evolution time for the differential equation
  statistic_equil_params%t_max = rdxx_cfg%max_evol_time
  !
  ! Prepare for the storage
  statistic_equil_params%NEQ = a_mol_using%n_level
  statistic_equil_params%LIW = 50 + statistic_equil_params%NEQ*2
  statistic_equil_params%LRW = 61 + 13*statistic_equil_params%NEQ*2 + &
                               statistic_equil_params%NEQ*statistic_equil_params%NEQ
  if (statistic_equil_params%NEQ .gt. a_mol_using%n_level) then
    if (allocated(statistic_equil_params%IWORK)) then
      deallocate(statistic_equil_params%IWORK, statistic_equil_params%RWORK)
    end if
  end if
  if (.not. allocated(statistic_equil_params%IWORK)) then
    allocate(statistic_equil_params%IWORK(statistic_equil_params%LIW), &
             statistic_equil_params%RWORK(statistic_equil_params%LRW))
  end if
  !
  if (makeLUT_) then
    lam_min = minval(a_mol_using%rad_data%list%lambda) ! micron
    lam_max = maxval(a_mol_using%rad_data%list%lambda) ! micron
    lam_min = lam_min * (1D0 - 10D0 * 1D7/phy_SpeedOfLight_CGS)
    lam_max = lam_max * (1D0 + 10D0 * 1D7/phy_SpeedOfLight_CGS)
    !
    call make_local_cont_lut( &
      trim(combine_dir_filename(rdxx_cfg%dir_bg, rdxx_cfg%filename_bg)), &
      rdxx_cfg%provide_bgfile, &
      rdxx_cfg%Tbg, rdxx_cfg%nTbg, rdxx_cfg%Tbgscaling, &
      cont_lut_bg, lam_min, lam_max, n_cont_lam)
    !
    call make_local_cont_lut( &
      trim(combine_dir_filename(rdxx_cfg%dir_in, rdxx_cfg%filename_in)), &
      rdxx_cfg%provide_infile, &
      rdxx_cfg%Tin, rdxx_cfg%nTin, rdxx_cfg%Tinscaling, &
      cont_lut_in, lam_min, lam_max, n_cont_lam)
    !
    call make_local_cont_lut( &
      trim(combine_dir_filename(rdxx_cfg%dir_out, rdxx_cfg%filename_out)), &
      rdxx_cfg%provide_outfile, &
      rdxx_cfg%Tout, rdxx_cfg%nTout, rdxx_cfg%Toutscaling, &
      cont_lut_out, lam_min, lam_max, n_cont_lam)
  end if
  !
end subroutine my_radex_prepare


subroutine make_local_cont_lut(filename, usefile, Ts, nTs, scaling, &
  cont_lut, lam_min, lam_max, n)
  ! Prepare for the continuum radiation field (usually just cmb).
  character(len=*), intent(in) :: filename
  logical, intent(in) :: usefile
  integer, intent(in) :: nTs
  double precision, dimension(:), intent(in) :: Ts, scaling
  type(type_continuum_lut), intent(out) :: cont_lut
  double precision, intent(in) :: lam_min, lam_max
  integer, intent(in) :: n
  integer i, j
  double precision dlam, dlam0, lam_ratio, lam, freq
  integer nrow, ios, fUnit, idx
  double precision, dimension(:), allocatable :: bglam, bgval, bgalpha
  double precision frac, tmp, tmp1
  character(len=128) str
  character, parameter :: commentstr = '!'
  !
  if ((.not. usefile) .and. (nTs .lt. 1)) then
    cont_lut%n = 0
    return
  end if
  !
  if (.not. allocated(cont_lut%lam)) then
    cont_lut%n = n
    allocate(cont_lut%lam(cont_lut%n), &
             cont_lut%alpha(cont_lut%n), &
             cont_lut%J(cont_lut%n))
  end if
  !
  if (usefile) then
    nrow = GetFileLen_comment_blank(filename, commentstr)
    allocate(bglam(nrow), bgval(nrow), bgalpha(nrow))
    bglam = 0D0
    bgval = 0D0
    bgalpha = 0D0
    call openFileSequentialRead(fUnit, filename, 999, 1)
    i = 0
    ios = 0
    do
      i = i + 1
      read(fUnit, '(A128)', IOSTAT=ios) str
      if (ios .ne. 0) then
        exit
      end if
      if ((str(1:1) .ne. commentstr) .and. (str(1:1) .ne. ' ')) then
        read(str, '(3F16.6)', IOSTAT=ios) bglam(i), bgval(i), bgalpha(i)
        if (ios .ne. 0) then
          exit
        end if
      end if
    end do
    close(fUnit)
    if ((i-1) .ne. nrow) then
      write(*, '(A)') 'Error loading background file.'
      write(*, '(A)') 'Maybe something not correct in the file format.'
      write(*, '("nrow,i,ios=", I6, I6, I6)') nrow, i, ios
    end if
  end if
  !
  dlam = (lam_max - lam_min) / dble(n-1)
  dlam0 = dlam / (1D0 + 5D0 * dlam / lam_min)
  lam_ratio = get_ratio_of_interval_log(lam_min, lam_max, dlam0, n)
  !
  dlam = dlam0
  cont_lut%lam(1) = lam_min
  !
  do i=1, cont_lut%n
    if (i .gt. 1) then
      cont_lut%lam(i) = cont_lut%lam(i-1) + dlam
      dlam = dlam * lam_ratio
    end if
    !
    cont_lut%alpha(i) = 0D0
    !
    lam = cont_lut%lam(i) + dlam * 0.5D0
    freq = phy_SpeedOfLight_CGS / (lam * phy_micron2cm)
    !
    ! Energy per unit area per unit frequency per second per sqradian
    cont_lut%J(i) = 0D0
    do j=1, nTs
      cont_lut%J(i) = cont_lut%J(i) + planck_B_nu(Ts(j), freq) * scaling(j)
    end do
    !
    if (usefile) then
      idx = binary_search(bglam, nrow, lam, 1)
      if ((idx .ge. 1) .and. (idx .le. nrow)) then
        if (idx .lt. nrow) then
          frac = (lam - bglam(idx)) / (bglam(idx+1) - bglam(idx))
          tmp = bgval(idx) * (1D0 - frac) + bgval(idx+1) * frac
          tmp1 = bgalpha(idx) * (1D0 - frac) + bgalpha(idx+1) * frac
        else
          frac = (lam - bglam(idx-1)) / (bglam(idx) - bglam(idx-1))
          tmp = bgval(idx) * frac + bgval(idx-1) * (1D0 - frac)
          tmp1 = bgalpha(idx) * frac + bgalpha(idx-1) * (1D0 - frac)
        end if
        cont_lut%J(i) = cont_lut%J(i) + tmp
        cont_lut%alpha(i) = cont_lut%alpha(i) + tmp1
      end if
    end if
    !write(*,*) i, cont_lut%lam(i), cont_lut%J(i), lam, lam_ratio, dlam, dlam0
  end do
  if (allocated(bglam)) then
    deallocate(bglam, bgval, bgalpha)
  end if
end subroutine make_local_cont_lut



function planck_B_lambda(T, lambda_CGS)
  double precision planck_B_lambda
  double precision, intent(in) :: T, lambda_CGS
  double precision tmp
  double precision, parameter :: THTINY = 1D-6
  tmp = (phy_hPlanck_CGS * phy_SpeedOfLight_CGS) / &
        (lambda_CGS * phy_kBoltzmann_CGS * T)
  if (abs(tmp) .gt. THTINY) then
    tmp = exp(tmp) - 1D0
  end if
  planck_B_lambda = &
    2D0*phy_hPlanck_CGS * phy_SpeedOfLight_CGS**2 / lambda_CGS**5 / tmp
end function planck_B_lambda



function planck_B_nu(T, nu)
  double precision planck_B_nu
  double precision, intent(in) :: T, nu
  double precision tmp
  double precision, parameter :: THTINY = 1D-6, THBIG = 100D0
  tmp = (phy_hPlanck_CGS*nu) / (phy_kBoltzmann_CGS*T)
  if (abs(tmp) .lt. THTINY) then
    planck_B_nu = 2D0*(nu/phy_SpeedOfLight_CGS)**2 * (phy_kBoltzmann_CGS*T)
  else if (tmp .gt. THBIG) then
    planck_B_nu = 2D0*phy_hPlanck_CGS * nu**3 / &
                  (phy_SpeedOfLight_CGS**2) * exp(-tmp)
  else
    planck_B_nu = 2D0*phy_hPlanck_CGS * nu**3 / &
                  (phy_SpeedOfLight_CGS**2 * (exp(tmp) - 1D0))
  end if
end function planck_B_nu


function calc_ortho_para_H2_ratio(T)
  ! Takahashi, J. 2001, ApJ, 561, 254
  double precision calc_ortho_para_H2_ratio
  double precision, intent(in) :: T
  double precision, parameter :: rotB = 87.6D0 ! K
  double precision, parameter :: thres = 20D0
  integer i, j1, j2
  double precision s1, s2, tmp1, tmp2, tt
  !
  tt = rotB / T
  !
  s1 = 0D0
  s2 = 0D0
  do i=0, 100
    j1 = 2*i + 1
    j2 = 2*i
    tmp1 = tt * dble(j1*(j1+1))
    tmp2 = tt * dble(j2*(j2+1))
    s1 = s1 + dble(2*j1+1) * exp(-tmp1)
    s2 = s2 + dble(2*j2+1) * exp(-tmp2)
    if ((tmp1 .gt. thres) .and. (tmp1 .gt. thres)) then
      exit
    end if
  end do
  calc_ortho_para_H2_ratio = 3D0 * s1 / s2
end function calc_ortho_para_H2_ratio


end module my_radex
