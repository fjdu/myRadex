import myRadex
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200000)

a = myRadex.MyRadexModel(
    dir_transition_rates='/Users/fjdu/_o/radex/data/',
    filename_molecule='hco+.dat')

a.run_one_params(
    Tkin=20.0, dv_CGS=1e0,
    dens_X_CGS=1e0, Ncol_X_CGS=1e13,
    H2_density_CGS=1e4, solve_method='ODE', geotype='sphere')

df = pd.DataFrame(data=a.data_transitions, columns=a.column_names)
print(df)
print('Cooling rate: ', a.cooling_rate)
