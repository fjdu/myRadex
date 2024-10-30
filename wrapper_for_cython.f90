module wrapper_for_cython

USE, INTRINSIC :: ISO_C_BINDING
implicit none

contains

subroutine copy_char_arr_to_str(src, dest, len)
  character(kind=c_char), dimension(*), intent(in) :: src
  character(len=*), intent(out) :: dest
  integer, intent(in) :: len
  integer i
  dest = ''
  do i=1,len
    dest(i:i) = src(i)
  end do
end subroutine copy_char_arr_to_str


subroutine copy_str_to_char_arr(src, dest, len)
  character(len=*), intent(in) :: src
  character(kind=c_char), dimension(len), intent(out) :: dest
  integer, intent(in) :: len
  integer i
  dest = ''
  do i=1,len
    dest(i) = src(i:i)
  end do
end subroutine copy_str_to_char_arr

subroutine f_config(&
    dir_transition_rates, &
    filename_molecule, &
    solve_method, &
    f_occupation_init_method, &
    o_column_names, o_molecule_name, &
    len1,len2,len3,len4, &
    Tbg, &
    beam_FWHM_in_arcsec,&
    verbose, &
    recalculateFreqWithEupElow, iLevel_subtract_one, &
    max_code_run_time, &
    max_evol_time, &
    rtol, atol, &
    n_levels, n_item, n_transitions, n_partners,len5,len6) bind(C, name="c_config")

  use myradex_wrapper, only: config_basic, column_names, &
    molecule_name, len_column_names, len_molecule_name
  use my_radex, only: rdxx_cfg

  character(kind=c_char), dimension(*), intent(in) :: &
      dir_transition_rates, filename_molecule, solve_method, f_occupation_init_method
  character(kind=c_char), dimension(128), intent(out) :: o_column_names
  character(kind=c_char), dimension(32), intent(out) :: o_molecule_name
  integer(kind=c_int), intent(in), value :: len1, len2, len3, len4
  real(kind=c_double), intent(in), value :: &
      Tbg, beam_FWHM_in_arcsec, max_code_run_time, max_evol_time, rtol, atol
  logical(kind=c_bool), intent(in), value :: &
      verbose, recalculateFreqWithEupElow, iLevel_subtract_one
  integer(kind=c_int), intent(out) :: n_levels, n_item, n_transitions, n_partners
  integer(kind=c_int), intent(out) :: len5, len6

  character(len=128) :: s1, s2, s3, s4

  call copy_char_arr_to_str(dir_transition_rates, s1, len1)
  call copy_char_arr_to_str(filename_molecule, s2, len2)
  call copy_char_arr_to_str(solve_method, s3, len3)
  call copy_char_arr_to_str(f_occupation_init_method, s4, len4)
  call config_basic(&
    s1, s2, &
    Tbg, &
    logical(verbose), &
    logical(recalculateFreqWithEupElow), logical(iLevel_subtract_one), &
    max_code_run_time, &
    max_evol_time, &
    rtol, atol, &
    s3, s4, &
    n_levels, n_item, n_transitions, n_partners)
  len5 = len(trim(column_names))
  len6 = len(trim(molecule_name))
  call copy_str_to_char_arr(column_names, o_column_names, len5)
  call copy_str_to_char_arr(molecule_name, o_molecule_name, len6)
  rdxx_cfg%beam_FWHM_in_arcsec = beam_FWHM_in_arcsec

end subroutine f_config


subroutine f_run_one_params(&
    Tkin, dv_CGS, &
    dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS, &
    n_levels, n_item, n_transitions, &
    energies, f_occupations, data_transitions, cooling_rate, &
    donotsolve, collisioPartnerCrit, &
    Tbg, beam_FWHM_in_arcsec, max_code_run_time, max_evol_time, &
    rtol, atol, solve_method, f_occupation_init_method, geotype, len3, len4, len5) bind(C, name="c_run_one_params")

  use myradex_wrapper, only: run_one_params

  real(kind=c_double), intent(in), value :: Tkin, dv_CGS, dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS
  integer(kind=c_int), intent(in), value :: n_levels, n_item, n_transitions, len3, len4, len5
  logical(kind=c_bool), intent(in), value :: donotsolve
  integer(kind=c_int), intent(in), value :: collisioPartnerCrit
  real(kind=c_double), dimension(n_levels), intent(out) :: energies, f_occupations
  real(kind=c_double), dimension(n_item, n_transitions), intent(out) :: data_transitions
  real(kind=c_double), intent(out) :: cooling_rate
  real(kind=c_double), intent(in), value :: Tbg, beam_FWHM_in_arcsec, max_code_run_time, max_evol_time, rtol, atol
  character(kind=c_char), dimension(*), intent(in) :: solve_method, f_occupation_init_method, geotype
  character(len=128) :: s3, s4, s5
  call copy_char_arr_to_str(solve_method, s3, len3)
  call copy_char_arr_to_str(f_occupation_init_method, s4, len4)
  call copy_char_arr_to_str(geotype, s5, len5)

  call run_one_params(&
    Tkin, dv_CGS, &
    dens_X_CGS, Ncol_X_CGS, &
    H2_density_CGS, HI_density_CGS, &
    oH2_density_CGS, pH2_densty_CGS, &
    HII_density_CGS, Electron_density_CGS, &
    n_levels, n_item, n_transitions, &
    energies, f_occupations, data_transitions, cooling_rate, &
    logical(donotsolve), collisioPartnerCrit, &
    Tbg, beam_FWHM_in_arcsec, max_code_run_time, max_evol_time, &
    rtol, atol, s3, s4, s5)

end subroutine f_run_one_params

end module wrapper_for_cython
