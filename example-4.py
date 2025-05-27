import myRadex
import numpy as np
import pandas as pd

molecule_params = {
    'pH2CO': {
       'config': dict(dir_transition_rates='./', filename_molecule='ph2co-h2.dat', f_occupation_init_method='Boltzmann', solve_method='ODE'),
       'params': [
                 dict(Tkin=1e2, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e4, geotype='lvg'),
                 dict(Tkin=1e3, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e5, geotype='lvg')
                 ]
      },
    'CI': {
       'config': dict(dir_transition_rates='./', filename_molecule='catom.dat', f_occupation_init_method='Boltzmann', solve_method='ODE'),
       'params': [
                 dict(Tkin=1e2, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e4, geotype='lvg'),
                 dict(Tkin=1e3, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e5, geotype='lvg')
                 ]
      }
}


for mname in molecule_params:
  print('Solving', mname)
  m = myRadex.MyRadexModel(**molecule_params[mname]['config'])
  for p in molecule_params[mname]['params']:
    m.run_one_params(**p)
    print('Cooling rate: ', m.cooling_rate)
    print('Did the solving proceed correctly?', m.flag_good)
    df = m.make_dataframe()
    print(df.describe())
