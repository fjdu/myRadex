import myRadex
import numpy as np
import pandas as pd
import os
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 999)
pd.set_option('display.max_seq_items', 999)
pd.set_option('display.width', 999)

a = myRadex.MyRadexModel(
    dir_transition_rates=os.path.expanduser('~/_c/my_radex/'),
    recalculateFreqWithEupElow=False,
    filename_molecule='oh2o@rovib.dat')
    # oh2o@rovib.dat  oh2o@daniel.dat

a.run_one_params(
    Tkin=1000.0, Ncol_X_CGS=1e17,
    H2_density_CGS=1e12, dv_FWHM_CGS=1e5, dens_X_CGS=1e0,
    beam_FWHM_in_arcsec=1e-1, geotype='slab',
    solve_method='ODE', f_occupation_init_method='bottom',
    max_code_run_time=20, max_evol_time=3.15e15,
    rtol=1e-4, atol=1e-13)

print('Did the solving proceed correctly?', a.flag_good)
print('Cooling rate (erg/cm^3/s): ', a.cooling_rate)

df = a.make_dataframe()

print(df[df['lam'].between(5,30) & (df['flux_Jy']>0.1)])
print(df.describe())
