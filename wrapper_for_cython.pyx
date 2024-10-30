# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.set cimport set as cppset
from libcpp.utility cimport pair
from libcpp cimport bool
import numpy as np
import cython

cdef extern from 'pipe_fortran_python.hpp':
  void  cc_config(
    const string dir_transition_rates,
    const string filename_molecule,
    const string solve_method,
    const string f_occupation_init_method,
    string& column_names,
    string& molecule_name,
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
    int *len6)

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
    const string solve_method,
    const string f_occupation_init_method,
    const string geotype)

def str2cppstr(s):
    if type(s) == str:
        return s.encode('utf-8')
    return s


cdef class MyRadexModel:

  cdef public:
      int n_levels, n_item, n_transitions, n_partners, len5, len6
      double cooling_rate
      object energies, f_occupations, data_transitions
      string c_column_names, c_molecule_name
      object column_names, molecule_name
      double [:] energies_view, f_occupations_view
      double [:,:] data_transitions_view
      object Tbg, beam_FWHM_in_arcsec, max_code_run_time, max_evol_time, rtol, atol, solve_method, f_occupation_init_method

  def __init__(self,
      dir_transition_rates=None,
      filename_molecule=None,
      solve_method='ODE',
      f_occupation_init_method='Boltzmann',
      Tbg=2.725,
      beam_FWHM_in_arcsec=30.0,
      verbose=True,
      recalculateFreqWithEupElow=False,
      iLevel_subtract_one=False,
      max_code_run_time=10.0,
      max_evol_time=3.14e12,
      rtol=1e-4,
      atol=1e-14):

      if not dir_transition_rates.endswith('/'):
        dir_transition_rates += '/'

      self.Tbg = Tbg
      self.beam_FWHM_in_arcsec = beam_FWHM_in_arcsec
      self.max_code_run_time = max_code_run_time
      self.max_evol_time = max_evol_time
      self.rtol = rtol
      self.atol = atol
      self.solve_method = solve_method
      self.f_occupation_init_method = f_occupation_init_method

      cc_config(str2cppstr(dir_transition_rates),
          str2cppstr(filename_molecule),
          str2cppstr(solve_method),
          str2cppstr(f_occupation_init_method),
          self.c_column_names,
          self.c_molecule_name,
          Tbg,
          beam_FWHM_in_arcsec,
          verbose,
          recalculateFreqWithEupElow,
          iLevel_subtract_one,
          max_code_run_time,
          max_evol_time,
          rtol,
          atol,
          &self.n_levels, &self.n_item,
          &self.n_transitions, &self.n_partners,
          &self.len5, &self.len6)

      self.energies = np.zeros(self.n_levels)
      self.f_occupations = np.zeros(self.n_levels)
      self.data_transitions = np.zeros((self.n_transitions, self.n_item))

      self.energies_view = self.energies
      self.f_occupations_view = self.f_occupations
      self.data_transitions_view = self.data_transitions

      self.column_names = self.c_column_names.decode('utf-8').split()
      self.molecule_name = self.c_molecule_name.decode('utf-8')

      return


  def run_one_params(self,
      Tkin=None, dv_CGS=None,
      dens_X_CGS=1e0, Ncol_X_CGS=None,
      H2_density_CGS=0e0, HI_density_CGS=0e0,
      oH2_density_CGS=0e0, pH2_densty_CGS=0e0,
      HII_density_CGS=0e0, Electron_density_CGS=0e0,
      donotsolve=False, collisioPartnerCrit=1,
      Tbg=None, beam_FWHM_in_arcsec=None, max_code_run_time=None, max_evol_time=None,
      rtol=None, atol=None, solve_method=None, f_occupation_init_method=None, geotype=None):

      smeth = str2cppstr(self.solve_method)
      if solve_method:
        smeth = str2cppstr(solve_method)
      fimeth = str2cppstr(self.f_occupation_init_method)
      if f_occupation_init_method:
        fimeth = str2cppstr(f_occupation_init_method)
      gtp = str2cppstr('')
      if geotype:
        gtp = str2cppstr(geotype)

      cc_run_one_params(
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
        self.n_levels,
        self.n_item,
        self.n_transitions,
        &self.energies_view[0],
        &self.f_occupations_view[0],
        &self.data_transitions_view[0][0],
        &self.cooling_rate,
        donotsolve,
        collisioPartnerCrit,
        Tbg or self.Tbg,
        beam_FWHM_in_arcsec or self.beam_FWHM_in_arcsec,
        max_code_run_time or self.max_code_run_time,
        max_evol_time or self.max_evol_time,
        rtol or self.rtol,
        atol or self.atol,
        smeth,
        fimeth,
        gtp)
      return
