#include <string>
#include <vector>

extern "C" {

//char c_column_names[128], c_molecule_name[128];

void c_config(
  const char *dir_transition_rates,
  const char *filename_molecule,
  const char *solve_method,
  const char *f_occupation_init_method,
  char* o_column_names,
  char* o_molecule_name,
  int len1, int len2, int len3, int len4,
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
  int *len6);

void c_run_one_params(
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
    const char *solve_method,
    const char *f_occupation_init_method,
    const char *geotype,
    int len3,
    int len4,
    int len5);

void c_get_flag(bool *flag);

}


void  cc_config(
  const std::string dir_transition_rates,
  const std::string filename_molecule,
  const std::string solve_method,
  const std::string f_occupation_init_method,
  std::string& o_column_names,
  std::string& o_molecule_name,
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
  int *len6);


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
    const std::string solve_method,
    const std::string f_occupation_init_method,
    const std::string geotype);


void cc_get_flag(bool *flag);
