import myRadex
import numpy as np
import pandas as pd

a = myRadex.MyRadexModel(
    dir_transition_rates='./',
    filename_molecule='ph2co-h2.dat')

a.run_one_params(
    Tkin=30.0, dv_CGS=1e5,
    dens_X_CGS=1e0, Ncol_X_CGS=1e13,
    oH2_density_CGS=1e5, pH2_densty_CGS=3e4, solve_method='ODE')

print(len(a.column_names))
print(a.data_transitions[:2, :])
print('Cooling rate: ', a.cooling_rate)
print('Did the solving proceed correctly?', a.flag_good)

df = pd.DataFrame(data=a.data_transitions, columns=a.column_names)

print(df)
