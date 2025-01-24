import myRadex
import numpy as np
import pandas as pd

a = myRadex.MyRadexModel(
    dir_transition_rates='./',
    filename_molecule='12C16O_H2.dat')

print('Quantum numbers:', a.qnum_s)

a.run_one_params(
    Tkin=167.3470613480748, dv_FWHM_CGS=1e5,
    dens_X_CGS=1e0, Ncol_X_CGS=2.1e20,
    H2_density_CGS=206913.808111479,
    geotype='lvg',
    solve_method='ODE', f_occupation_init_method='bottom')

print('Cooling rate: ', a.cooling_rate)
print('Did the solving proceed correctly?', a.flag_good)

df = a.make_dataframe()

print(df)
print(df.describe())
