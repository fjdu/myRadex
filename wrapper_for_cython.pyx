# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.set cimport set as cppset
from libcpp.utility cimport pair
from libcpp cimport bool
import numpy as np
import cython
import multiprocessing as mp
import uuid

cdef extern from 'pipe_fortran_python.hpp':
  void  cc_config(
    const string dir_transition_rates,
    const string filename_molecule,
    const string solve_method,
    const string f_occupation_init_method,
    string& column_names,
    string& molecule_name,
    vector[string] *sQnum_s,
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
    double He_density_CGS,
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

  void cc_get_flag(bool *flag_good)

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
      vector[string] qnum_s
      object Tbg, beam_FWHM_in_arcsec, max_code_run_time, max_evol_time, rtol, atol, solve_method, f_occupation_init_method
      bool flag_good

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
      max_evol_time=3.15e12,
      rtol=1e-4,
      atol=1e-14):
      '''
      Initialize the code.

      The most important thing done in this step is to load the transition data file.
      Many other parameters (such as `rtol`, `atol`) can be adjusted when calling `run_one_params`.

      Parameters
      ----------
      dir_transition_rates: str
          The directory containing the transition data file.
      filename_molecule: str
          The filename of the transition data file.
      solve_method: str, optional
          The method to use for solving the statistical equilibrium problem.
          Possible values: 'ODE', 'Newton'
          Default: 'ODE'
      f_occupation_init_method: str, optional
          The method to use for initializing the occupation fractions.
          Possible values: 'Boltzmann', 'Random'
          Default: 'Boltzmann'
      Tbg: float, optional
          Background temperature.
          Default: 2.725 (= T_cmb)
      beam_FWHM_in_arcsec: float, optional
          Beam full-width-at-half-maximum, in arcsec.
          Default: 30.0
      verbose: bool, optional
          Whether to print running messages in the terminal screen.
          Default: True
      recalculateFreqWithEupElow: bool, optional
          Whether to recalculate the transition frequencies using the Eup and Elow in the data file.
          Default: False
      iLevel_subtract_one: bool, optional
          Whether to subtract the energy level numbers by one (to meet the Python convention).
          Default: False
      max_code_run_time: float, optional
          Max code run time for one model run (only for the ODE method).
          Default: 10.0 seconds
      max_evol_time: float, optional
          Max evolution time when solving with the ODE approach
          Default: 3.15e12 seconds (= 1e5 yr)
      rtol: float, optional
          Relative tolerance of the solution
          Default: 1e-4
      atol: float, optional
          Absolute tolerance of the solution
          Default: 1e-14
      '''

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
      self.flag_good = True

      cc_config(str2cppstr(dir_transition_rates),
          str2cppstr(filename_molecule),
          str2cppstr(solve_method),
          str2cppstr(f_occupation_init_method),
          self.c_column_names,
          self.c_molecule_name,
          &self.qnum_s,
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
      Tkin=None, dv_FWHM_CGS=None, dv_CGS=None,
      dens_X_CGS=1e0, Ncol_X_CGS=None, geotype='slab',
      H2_density_CGS=0e0, HI_density_CGS=0e0,
      oH2_density_CGS=0e0, pH2_densty_CGS=0e0,
      HII_density_CGS=0e0, Electron_density_CGS=0e0, He_density_CGS=0e0,
      donotsolve=False, collisioPartnerCrit=1,
      Tbg=None, beam_FWHM_in_arcsec=None, max_code_run_time=None, max_evol_time=None,
      rtol=None, atol=None, solve_method=None, f_occupation_init_method=None):
      '''
      Run a model for one set of parameters.

      Parameters
      ----------
      Tkin: float
          Kinetic temperature, in K
      dv_FWHM_CGS: float
          Velocity FWHM, in cm/s
      dv_CGS: float, obsolete
          The same as dv_FWHM_CGS.
      dens_X_CGS: float, optional
          Volumn density of the molecule under study, in cm^-3.
          This parameter is not essential.
          The code use this density and its column density to get a length scale.
          Default: 1.0
      Ncol_X_CGS: float
          Column density of the molecule under study, in cm^-2.
      geotype: str, optional
          Geometric type of the model.
          Default: 'slab'
      H2_density_CGS: float, optional
          The density of H2 molecule, in cm^-3.
          Default: 0
      HI_density_CGS: float, optional
          The density of atomic H, in cm^-3.
          Default: 0
      oH2_density_CGS: float, optional
          The density of ortho H2, in cm^-3.
          Default: 0
      pH2_density_CGS: float, optional
          The density of para H2, in cm^-3.
          Default: 0
      HII_density_CGS: float, optional
          The density of H+, in cm^-3.
          Default: 0
      Electron_density_CGS: float, optional
          The density of electron in cm^-3.
          Default: 0
      He_density_CGS: float, optional
          The density of Helium in cm^-3.
          Default: 0
      donotsolve: bool, optional
          If True, the code will not solve the problem (but will calculate the critical densities).
          Default: False
      collisioPartnerCrit: int, optional
          Select which collisional partner (when there are more than one in the file) to use for calculqting the critical densities.
          Default: 1
      Tbg: float, optional
          Background temperature.
          Default: 2.725 (= T_cmb)
      beam_FWHM_in_arcsec: float, optional
          Beam full-width-at-half-maximum, in arcsec.
          Default: 30.0
      max_code_run_time: float, optional
          Max code run time for one model run (only for the ODE method).
          Default: 10.0 seconds
      max_evol_time: float, optional
          Max evolution time when solving with the ODE approach
          Default: 3.15e12 seconds (= 1e5 yr)
      rtol: float, optional
          Relative tolerance of the solution
          Default: 1e-4
      atol: float, optional
          Absolute tolerance of the solution
          Default: 1e-14
      solve_method: str, optional
          The method to use for solving the statistical equilibrium problem.
          Possible values: 'ODE', 'Newton'
          Default: 'ODE'
      f_occupation_init_method: str, optional
          The method to use for initializing the occupation fractions.
          Possible values: 'Boltzmann', 'Random'
          Default: 'Boltzmann'
      '''

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
        dv_FWHM_CGS or dv_CGS,
        dens_X_CGS,
        Ncol_X_CGS,
        H2_density_CGS,
        HI_density_CGS,
        oH2_density_CGS,
        pH2_densty_CGS,
        HII_density_CGS,
        Electron_density_CGS,
        He_density_CGS,
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
      cc_get_flag(&self.flag_good)
      return


  def __make_qnum__(self):
    res = []
    for i in range(self.data_transitions.shape[0]):
      iup = int(self.data_transitions[i,0]) - 1
      ilow = int(self.data_transitions[i,1]) - 1
      res.append(self.qnum_s[iup].decode('utf-8') + ' -> ' + self.qnum_s[ilow].decode('utf-8'))
    return res


  def make_dataframe(self):
    '''
    Return a `pandas` dataframe from the calculated transition data for viewing.
    '''
    import pandas as pd
    df = pd.DataFrame(data=self.data_transitions, columns=self.column_names)
    df.insert(len(self.column_names), 'q_num', self.__make_qnum__())
    return df


class SolverProcess:
    def __init__(self, init_params):
        self.input_queue = mp.Queue()
        self.output_queue = mp.Queue()
        self.process = mp.Process(
            target=self._worker, 
            args=(self.input_queue, self.output_queue, init_params)
        )
        self.process.start()
    
    @staticmethod
    def _worker(input_queue, output_queue, init_params):
        solver = MyRadexModel(**init_params)

        while True:
            task = input_queue.get()
            if task is None:
                break
            
            print('Running task', task)
            task_id, data = task
            try:
                solver.run_one_params(**data)
                result = solver.make_dataframe()
                output_queue.put((task_id, 'success', result))
            except Exception as e:
                output_queue.put((task_id, 'error', str(e)))
    
    def solve(self, data):
        task_id = str(uuid.uuid4())
        self.input_queue.put((task_id, data))
        
        while True:
            result_id, status, result = self.output_queue.get()
            if result_id == task_id:
                if status == 'error':
                    raise RuntimeError(result)
                return result
    
    def cleanup(self):
        self.input_queue.put(None)
        self.process.join()


class MyRadexManager:
    '''
    Example usage:
    if __name__ == '__main__':
        manager = myRadex.MyRadexManager()
        
        solvers = [
            manager.create_solver('12CO', dir_transition_rates='./', filename_molecule='12C16O_H2.dat'),
            manager.create_solver('CI', dir_transition_rates='./', filename_molecule='catom.dat')
        ]
        
        p1 = dict(Tkin=1e2, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e4, geotype='lvg')
        p2 = dict(Tkin=1e3, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e5, geotype='lvg')
        
        result1 = manager.solve('12CO', p1)
        result2 = manager.solve('CI', p1) 
        
        manager.cleanup()
    
        print(result1.describe())
        print(result2.describe())
    '''

    def __init__(self):
        self.solvers = {}
    
    def create_solver(self, solver_id, **init_params):
        if solver_id in self.solvers:
            raise ValueError(f"Solver {solver_id} already exists")
        
        self.solvers[solver_id] = SolverProcess(init_params)
        return solver_id
    
    def solve(self, solver_id, data):
        """Solve using a specific solver instance"""
        if solver_id not in self.solvers:
            raise ValueError(f"Solver {solver_id} not found")
        
        return self.solvers[solver_id].solve(data)
    
    def cleanup(self):
        for solver in self.solvers.values():
            solver.cleanup()
        self.solvers.clear()

