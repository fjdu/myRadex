import myRadex

if __name__ == '__main__':
    manager = myRadex.MyRadexManager()
    
    solvers = [
        manager.create_solver('12CO', dir_transition_rates='./', filename_molecule='12C16O_H2.dat', solve_method='Newton', f_occupation_init_method='Boltzmann'),
        manager.create_solver('CI', dir_transition_rates='./', filename_molecule='catom.dat', solve_method='Newton', f_occupation_init_method='Boltzmann')
    ]
    
    p1 = dict(Tkin=1e2, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e4, geotype='lvg')
    p2 = dict(Tkin=1e3, dv_FWHM_CGS=3e5, dens_X_CGS=1e0, Ncol_X_CGS=1e16, H2_density_CGS=1e5, geotype='lvg')
    
    result1 = manager.solve('12CO', p1)
    result2 = manager.solve('CI', p1) 
    
    manager.cleanup()

    print(result1.describe())
    print(result2.describe())
