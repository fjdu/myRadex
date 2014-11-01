module myradex_wrapper

implicit none


logical :: flag_good

integer, parameter :: n_item_column = 19

character(len=64) :: &
    About = 'Author: Fujun Du (fujun.du@gmail.com, fdu@umich.edu)'

character(len=128) :: column_names = &
    'iup' //' '// 'ilow' //' '// 'Eup' //' '// 'freq' //' '// 'lam' //' ' &
    // 'Tex' //' '// 'tau' //' '// 'Tr' //' '// &
    'fup' //' '// 'flow' //' '// 'flux_K' //' '// 'flux' //' '// 'beta' //' ' &
    // 'Jnu' //' '// 'gup' //' '// 'glow' //' '// 'Aul' //' '// 'Bul' //' '// 'Blu'

character(len=32) :: molecule_name = ''

contains


subroutine config_basic(dir_transition_rates, filename_molecule, verbose, &
    n_levels, n_item, n_transitions)
  use my_radex
  character(len=*), intent(in) :: dir_transition_rates, filename_molecule
  logical, intent(in) :: verbose
  integer, intent(out) :: n_levels, n_item, n_transitions
  !
  rdxx_cfg%dir_transition_rates = dir_transition_rates
  rdxx_cfg%filename_molecule = filename_molecule
  rdxx_cfg%verbose = verbose
  !
  if (verbose) then
    write(*, '(A, A)') 'About this tool: ', About
    write(*, '(A, A)') 'Column names of the output: ', column_names
  end if
  !
  call my_radex_prepare
  !
  n_levels = a_mol_using%n_level
  n_item = n_item_column
  n_transitions = a_mol_using%rad_data%n_transition
  molecule_name = a_mol_using%name_molecule
  !
  if (verbose) then
    write(*, '(A, A)') 'Molecule name: ', molecule_name
  end if
end subroutine config_basic


subroutine run_one_params( &
    Tkin, dv_CGS, &
    dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_density_CGS, &
    HII_density_CGS, Electron_density_CGS, &
    n_levels, n_item, n_transitions, &
    energies, f_occupations, data_transitions, cooling_rate)
  !
  use my_radex
  use statistic_equilibrium
  !
  double precision, intent(in) :: Tkin, dv_CGS, dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_density_CGS, &
    HII_density_CGS, Electron_density_CGS
  !
  integer, intent(in) :: n_levels, n_item, n_transitions
  double precision, dimension(n_levels), intent(out) :: energies, f_occupations
  double precision, dimension(n_item, n_transitions), intent(out) :: data_transitions
  double precision, intent(out) :: cooling_rate
  !
  double precision fup, flow, gup, glow, Tex, Tr, flux_CGS, flux_K_km_s
  double precision Inu_t, tau, t1, t2
  integer i
  !
  rdxx_cfg%nTkin   = 1
  rdxx_cfg%ndv     = 1
  rdxx_cfg%nn_x    = 1
  rdxx_cfg%nNcol_x = 1
  rdxx_cfg%ndens   = 1
  !
  rdxx_cfg%iTkin   = 1
  rdxx_cfg%idv     = 1
  rdxx_cfg%in_x    = 1
  rdxx_cfg%iNcol_x = 1
  rdxx_cfg%idens   = 1
  !
  rdxx_cfg%Tkin(1) = Tkin
  rdxx_cfg%dv(1) = dv_CGS
  rdxx_cfg%n_x(1) = dens_X_CGS
  rdxx_cfg%Ncol_x(1) = Ncol_X_CGS
  !
  rdxx_cfg%n_H2(1) = H2_density_CGS
  rdxx_cfg%n_HI(1) = HI_density_CGS
  rdxx_cfg%n_oH2(1) = oH2_density_CGS
  rdxx_cfg%n_pH2(1) = pH2_density_CGS
  rdxx_cfg%n_Hplus(1) = HII_density_CGS
  rdxx_cfg%n_E(1) = Electron_density_CGS
  !
  call my_radex_prepare_molecule
  call statistic_equil_solve
  call calc_cooling_rate
  !
  flag_good = statistic_equil_params%is_good
  !
  cooling_rate = a_mol_using%cooling_rate_total
  !
  energies = a_mol_using%level_list%energy
  f_occupations = a_mol_using%f_occupation
  !
  rdxx_cfg%freqmin = 1D-99
  rdxx_cfg%freqmax = 1D99
  !
  do i=1, n_transitions
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
        a_mol_using%dv * r%freq / phy_SpeedOfLight_CGS
      !
      data_transitions(:, i) = &
        (/ &
        dble(r%iup-1), dble(r%ilow-1), r%Eup, r%freq, r%lambda, Tex, r%tau, Tr, &
        fup, flow, flux_K_km_s, flux_CGS, r%beta, &
        r%J_ave, gup, glow, r%Aul, r%Bul, r%Blu /)
    end associate
  end do
  !
end subroutine run_one_params


end module myradex_wrapper
