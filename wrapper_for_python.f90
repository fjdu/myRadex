module myradex_wrapper

implicit none


logical :: flag_good

integer, parameter :: n_item_column = 21

character(len=64) :: &
    About = 'Author: Fujun Du (fjdu@pmo.ac.cn, fujun.du@gmail.com)'

character(len=128) :: column_names = &
    'iup' //' '// 'ilow' //' '// 'Eup' //' '// 'freq' //' '// 'lam' //' ' &
    // 'Tex' //' '// 'tau' //' '// 'Tr' //' '// &
    'fup' //' '// 'flow' //' '// 'flux_K' //' '// 'flux_int' // ' ' // 'flux_Jy' //' '// 'beta' //' ' &
    // 'Jnu' //' '// 'gup' //' '// 'glow' //' '// 'Aul' //' '// 'Bul' //' '// 'Blu' // ' ' // 'n_crit' // ' ' // 'n_crit_old'

character(len=32) :: molecule_name = ''

contains


subroutine config_basic(dir_transition_rates, filename_molecule, &
    Tbg, &
    verbose, &
    recalculateFreqWithEupElow, iLevel_subtract_one, &
    max_code_run_time, &
    max_evol_time, &
    rtol, atol, &
    solve_method, &
    n_levels, n_item, n_transitions, n_partners)
  use my_radex
  character(len=*), intent(in) :: dir_transition_rates, filename_molecule
  double precision, intent(in) :: Tbg
  logical, intent(in), optional :: verbose
  logical, intent(in), optional :: recalculateFreqWithEupElow, iLevel_subtract_one
  double precision, intent(in), optional :: max_code_run_time, max_evol_time, rtol, atol
  character(len=*), intent(in), optional :: solve_method
  integer, intent(out) :: n_levels, n_item, n_transitions, n_partners
  logical verbs
  !
  if (present(recalculateFreqWithEupElow)) then
    rdxx_cfg%recalculateFreqWithEupElow = recalculateFreqWithEupElow
  end if
  if (present(iLevel_subtract_one)) then
    rdxx_cfg%iLevel_subtract_one = iLevel_subtract_one
  end if
  !
  if (present(max_code_run_time)) then
    rdxx_cfg%max_code_run_time = max_code_run_time
  end if
  if (present(max_evol_time)) then
    rdxx_cfg%max_evol_time = max_evol_time
  end if
  if (present(rtol)) then
    rdxx_cfg%rtol = rtol
  end if
  if (present(atol)) then
    rdxx_cfg%atol = atol
  end if
  if (present(solve_method)) then
    rdxx_cfg%solve_method = solve_method
  end if
  !
  rdxx_cfg%dir_transition_rates = dir_transition_rates
  rdxx_cfg%filename_molecule = filename_molecule
  !
  rdxx_cfg%nTbg = 1
  rdxx_cfg%Tbg(1) = Tbg
  !
  verbs = .true.
  !
  if (present(verbose)) then
    verbs = verbose
  end if
  !
  rdxx_cfg%verbose = verbs
  !
  if (verbs) then
    write(*, '(A, A)') 'Column names of the output: ', column_names
  end if
  !
  call my_radex_prepare
  !
  n_levels = a_mol_using%n_level
  n_item = n_item_column
  n_transitions = a_mol_using%rad_data%n_transition
  molecule_name = a_mol_using%name_molecule
  n_partners = a_mol_using%colli_data%n_partner
  !
  if (verbs) then
    write(*, '(A, A)') 'Molecule name: ', molecule_name
    write(*, '(A, I10)') 'Number of levels: ', n_levels
    write(*, '(A, I10)') 'Number of columns: ', n_item
    write(*, '(A, I10)') 'Number of transitions: ', n_transitions
  end if
end subroutine config_basic


subroutine run_one_params( &
    Tkin, dv_CGS, &
    dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS, &
    n_levels, n_item, n_transitions, &
    donotsolve, &
    collisioPartnerCrit, &
    energies, f_occupations, data_transitions, cooling_rate)
  !
  use my_radex
  use statistic_equilibrium
  !
  double precision, intent(in) :: Tkin, dv_CGS, dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS
  !
  integer, intent(in) :: n_levels, n_item, n_transitions
  logical, intent(in), optional :: donotsolve
  integer, intent(in), optional :: collisioPartnerCrit
  double precision, dimension(n_levels), intent(out) :: energies, f_occupations
  double precision, dimension(n_item, n_transitions), intent(out) :: data_transitions
  double precision, intent(out) :: cooling_rate
  !
  type(type_rad_transition) r
  double precision fup, flow, gup, glow, Tex, Tr, flux_CGS, flux_K_km_s
  double precision Inu_t, tau, t1, t2
  double precision critical_density, critical_density_old
  integer i, iupSav, ilowSav
  logical donotsolve_
  integer iCollPartner
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
  rdxx_cfg%n_pH2(1) = pH2_densty_CGS
  rdxx_cfg%n_Hplus(1) = HII_density_CGS
  rdxx_cfg%n_E(1) = Electron_density_CGS

  if (present(donotsolve)) then
    donotsolve_ = donotsolve
  else
    donotsolve_ = .false.
  end if
  !
  iCollPartner = 1
  if (present(collisioPartnerCrit)) then
    iCollPartner = collisioPartnerCrit
  end if
  !
  call my_radex_prepare_molecule
  if (.not. donotsolve_) then
    if (trim(rdxx_cfg%solve_method) .eq. 'ODE') then
      call statistic_equil_solve
    else if (trim(rdxx_cfg%solve_method) .eq. 'Newton') then
      call statistic_equil_solve_Newton
    else
      write(*,*) 'Unknown solving method: ', trim(rdxx_cfg%solve_method)
    end if
  end if
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
      r = a_mol_using%rad_data%list(i)
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
      iupSav = r%iup
      ilowSav = r%ilow
      if (rdxx_cfg%iLevel_subtract_one) then
        iupSav = r%iup - 1
        ilowSav = r%ilow - 1
      end if
      !
      call calc_critical_density_for_one_transition(i, tau)
      critical_density = a_mol_using%rad_data%list(i)%critical_densities(iCollPartner)
      call calc_critical_density_old_def_for_one_transition(i, tau)
      critical_density_old = a_mol_using%rad_data%list(i)%critical_densities(iCollPartner)
      !
      data_transitions(:, i) = &
        (/ &
        dble(iupSav), dble(ilowSav), r%Eup, r%freq, r%lambda, Tex, r%tau, Tr, &
        fup, flow, flux_K_km_s, flux_CGS, r%beta, &
        r%J_ave, gup, glow, r%Aul, r%Bul, r%Blu, critical_density, critical_density_old /)
  end do
end subroutine run_one_params


subroutine run_one_params_geometry( &
    Tkin, dv_CGS, &
    dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS, &
    n_levels, n_item, n_transitions, &
    geotype, &
    donotsolve, &
    collisioPartnerCrit, &
    energies, f_occupations, data_transitions, cooling_rate)
  !
  use my_radex
  use statistic_equilibrium
  !
  double precision, intent(in) :: Tkin, dv_CGS, dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS
  !
  integer, intent(in) :: n_levels, n_item, n_transitions
  character(len=16), intent(in) :: geotype
  logical, intent(in), optional :: donotsolve
  integer, intent(in), optional :: collisioPartnerCrit
  !
  double precision, dimension(n_levels), intent(out) :: energies, f_occupations
  double precision, dimension(n_item, n_transitions), intent(out) :: data_transitions
  double precision, intent(out) :: cooling_rate
  logical donotsolve_
  integer iCollPartner
  !
  if (present(donotsolve)) then
    donotsolve_ = donotsolve
  else
    donotsolve_ = .false.
  end if
  !
  iCollPartner = 1
  if (present(collisioPartnerCrit)) then
    iCollPartner = collisioPartnerCrit
  end if
  !
  rdxx_cfg%geotype = adjustl(geotype)
  if (rdxx_cfg%verbose) then
    write(*,*) 'Using geotype:', rdxx_cfg%geotype
  end if
  !
  call run_one_params( &
    Tkin, dv_CGS, &
    dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS, &
    n_levels, n_item, n_transitions, donotsolve_, iCollPartner, &
    energies, f_occupations, data_transitions, cooling_rate)
end subroutine run_one_params_geometry


subroutine calc_critical_density(tau, n_partner, n_transitions, &
    critical_densities, iup, ilow)
  use statistic_equilibrium
  double precision, intent(in) :: tau
  integer, intent(in) :: n_partner, n_transitions
  double precision, dimension(n_partner,n_transitions), intent(out) :: critical_densities
  integer, dimension(n_transitions), intent(out) :: iup, ilow
  integer ipt, itr
  call calc_critical_density_f(tau)
  do ipt=1, n_partner
    do itr=1, n_transitions
      critical_densities(ipt, itr) = a_mol_using%rad_data%list(itr)%critical_densities(ipt)
    end do
  end do
  do itr=1, n_transitions
    iup(itr) = a_mol_using%rad_data%list(itr)%iup
    ilow(itr) = a_mol_using%rad_data%list(itr)%ilow
  end do
end subroutine calc_critical_density



subroutine calc_critical_density_old_def(tau, n_partner, n_transitions, &
    critical_densities, iup, ilow)
  use statistic_equilibrium
  double precision, intent(in) :: tau
  integer, intent(in) :: n_partner, n_transitions
  double precision, dimension(n_partner,n_transitions), intent(out) :: critical_densities
  integer, dimension(n_transitions), intent(out) :: iup, ilow
  integer ipt, itr
  call calc_critical_density_old_def_f(tau)
  do ipt=1, n_partner
    do itr=1, n_transitions
      critical_densities(ipt, itr) = a_mol_using%rad_data%list(itr)%critical_densities(ipt)
    end do
  end do
  do itr=1, n_transitions
    iup(itr) = a_mol_using%rad_data%list(itr)%iup
    ilow(itr) = a_mol_using%rad_data%list(itr)%ilow
  end do
end subroutine calc_critical_density_old_def

end module myradex_wrapper
