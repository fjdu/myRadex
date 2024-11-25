#include "pipe_fortran_python.hpp"
#include <iostream>

void  cc_config(
  const std::string dir_transition_rates,
  const std::string filename_molecule,
  const std::string solve_method,
  const std::string f_occupation_init_method,
  std::string& o_column_names,
  std::string& o_molecule_name,
  std::vector<std::string> *sQnum_s,
  double tbg,
  double beam_FWHM_in_arcsec,
  bool verbose,
  bool recalculatefreqwitheupelow,
  bool ilevel_subtract_one,
  double max_code_run_time,
  double max_evol_time,
  double rtol,
  double atol,
  int *n_levels,
  int *n_item,
  int *n_transitions,
  int *n_partners,
  int *len5,
  int *len6) {

  char s1[256]="", s2[128]="";

  c_config(
    dir_transition_rates.c_str(),
    filename_molecule.c_str(),
    solve_method.c_str(),
    f_occupation_init_method.c_str(),
    s1,
    s2,
    dir_transition_rates.length(),
    filename_molecule.length(),
    solve_method.length(),
    f_occupation_init_method.length(),
    tbg,
    beam_FWHM_in_arcsec,
    verbose,
    recalculatefreqwitheupelow,
    ilevel_subtract_one,
    max_code_run_time,
    max_evol_time,
    rtol,
    atol,
    n_levels,
    n_item,
    n_transitions,
    n_partners,
    len5,
    len6);
  o_column_names = std::string(s1);
  o_molecule_name = std::string(s2);

  char** s3 = new char*[*n_levels];
  for (int i=0; i<*n_levels; ++i) {
    s3[i] = nullptr;
  }

  c_get_sQum_s(s3, *n_levels);
  for (int i=0;i<*n_levels;++i) {
    (*sQnum_s).push_back(std::string(s3[i]));
  }
}


void cc_run_one_params(
    double Tkin,
    double dv_CGS,
    double dens_X_CGS,
    double Ncol_X_CGS,
    double H2_density_CGS,
    double HI_density_CGS,
    double oH2_density_CGS,
    double pH2_densty_CGS,
    double HII_density_CGS,
    double Electron_density_CGS,
    int n_levels,
    int n_item,
    int n_transitions,
    double *energies,
    double *f_occupations,
    double *data_transitions,
    double *cooling_rate,
    bool donotsolve,
    int collisioPartnerCrit,
    double Tbg,
    double beam_FWHM_in_arcsec,
    double max_code_run_time,
    double max_evol_time,
    double rtol,
    double atol,
    std::string solve_method,
    std::string f_occupation_init_method,
    std::string geotype) {
// The memory needed by energies, f_occupations, and data_transitions are assumed to have already been allocated.

  c_run_one_params(
    Tkin,
    dv_CGS,
    dens_X_CGS,
    Ncol_X_CGS,
    H2_density_CGS,
    HI_density_CGS,
    oH2_density_CGS,
    pH2_densty_CGS,
    HII_density_CGS,
    Electron_density_CGS,
    n_levels,
    n_item,
    n_transitions,
    energies,
    f_occupations,
    data_transitions,
    cooling_rate,
    donotsolve,
    collisioPartnerCrit,
    Tbg,
    beam_FWHM_in_arcsec,
    max_code_run_time,
    max_evol_time,
    rtol,
    atol,
    solve_method.c_str(),
    f_occupation_init_method.c_str(),
    geotype.c_str(),
    solve_method.length(),
    f_occupation_init_method.length(),
    geotype.length());
}


void cc_get_flag(bool *flag) {
  c_get_flag(flag);
}
