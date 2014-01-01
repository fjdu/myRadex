module my_radex

use phy_const
!use data_struct
use trivials
use statistic_equilibrium
!use chemistry

implicit none

type :: type_my_radex_config
  character(len=128) :: dir_transition_rates = './'
  character(len=128) :: dir_save = './'
  character(len=128) :: filename_molecule = '12C16O_H2.dat'
  character(len=128) :: filename_save = 'output.dat'
  character(len=16)  :: geotype = ''
  double precision Tkin, Tcmb, dv, length_scale
  logical :: opH2eq3 = .true., provideLength = .false.
  double precision n_H, X_H2, X_oH2, X_pH2, X_HI, X_Hplus, X_E
  double precision X_x, Ncol_x
  double precision max_evol_time
end type type_my_radex_config

type(type_molecule_energy_set), target :: molecule

type(type_my_radex_config) :: my_radex_config

namelist /my_radex_configure/ &
  my_radex_config


contains


subroutine do_my_radex
  integer i, fU
  double precision fup, flow, gup, glow, Tex, Tr, flux_CGS, flux_K_km_s
  double precision Inu_t, tau, t1, t2
  !
  a_mol_using => molecule
  call my_radex_prepare_molecule
  call statistic_equil_solve
  !
  call openFileSequentialWrite(fU, &
    combine_dir_filename(my_radex_config%dir_save, &
    my_radex_config%filename_save), 999, 1)
  write(fU, '(2A5, A12, 2A15, 9A12, 2A7, &
            &3A12)') &
    '! iup', 'ilow', 'Eup', 'freq', 'lam', 'Tex', 'tau', 'Tr', &
    'fup', 'flow', 'flux_K', 'flux', 'beta', &
    'Jnu', 'gup', 'glow', 'Aul', 'Bul', 'Blu'
  write(fU, '(2A5, A12, 2A15, 9A12, 2A7, &
            &3A12)') &
    '!    ', '  ', 'K', 'Hz', 'micron', 'K', '', 'K', &
    '   ', '    ', 'K km/s', 'erg/cm2/s', '    ', &
    '...', '   ', '    ', '...', '...', '...'
  do i=1, a_mol_using%rad_data%n_transition
    associate(r => a_mol_using%rad_data%list(i))
      fup  = a_mol_using%f_occupation(r%iup)
      flow = a_mol_using%f_occupation(r%ilow)
      gup  = a_mol_using%level_list(r%iup)%weight
      glow = a_mol_using%level_list(r%ilow)%weight
      Tex = -(r%Eup - r%Elow) / log(fup*glow / (flow*gup))
      !
      tau = r%tau
      t1 = exp(-tau)
      if (abs(tau) .lt. 1D-6) then
        t2 = tau
      else
        t2 = 1D0 - t1
      end if
      Inu_t = planck_B_nu(Tex, r%freq) * t2 + t1 * r%J_cont
      !
      Tr = (Inu_t - r%J_cont) * phy_SpeedOfLight_CGS**2 / &
        (2D0 * r%freq**2 * phy_kBoltzmann_CGS)
      flux_K_km_s = Tr * a_mol_using%dv / 1D5 * sqrt(phy_pi/4D0/log(2D0))
      flux_CGS = (Inu_t - r%J_cont) * &
        a_mol_using%dv * r%freq / phy_SpeedOfLight_CGS
      write(fU, '(2I5, F12.4, 2ES15.7, 9ES12.3, 2F7.1, &
                &3ES12.3)') &
        r%iup-1, r%ilow-1, r%Eup, r%freq, r%lambda, Tex, r%tau, Tr, &
        fup, flow, flux_K_km_s, flux_CGS, r%beta, &
        r%J_ave, gup, glow, r%Aul, r%Bul, r%Blu
    end associate
  end do
  close(fU)
  nullify(a_mol_using)
end subroutine do_my_radex


subroutine my_radex_prepare
  double precision lam_min, lam_max
  ! Load the molecular data
  a_mol_using => molecule
  call load_moldata_LAMBDA(&
    combine_dir_filename(my_radex_config%dir_transition_rates, &
    my_radex_config%filename_molecule))
  !
  write(*, '(A, I5)') 'Number of levels: ', a_mol_using%n_level
  !
  ! Evolution time for the differential equation
  statistic_equil_params%t_max = my_radex_config%max_evol_time
  !
  ! Prepare for the storage
  statistic_equil_params%NEQ = a_mol_using%n_level
  statistic_equil_params%LIW = 20 + statistic_equil_params%NEQ
  statistic_equil_params%LRW = 22 + 9*statistic_equil_params%NEQ + &
                               statistic_equil_params%NEQ*statistic_equil_params%NEQ
  if (statistic_equil_params%NEQ .gt. a_mol_using%n_level) then
    if (allocated(statistic_equil_params%IWORK)) then
      deallocate(statistic_equil_params%IWORK, statistic_equil_params%RWORK)
    end if
  end if
  if (.not. allocated(statistic_equil_params%IWORK)) then
    allocate(statistic_equil_params%IWORK(statistic_equil_params%LIW), &
             statistic_equil_params%RWORK(statistic_equil_params%LRW))
  end if
  !
  lam_min = minval(a_mol_using%rad_data%list%lambda)/1.2D0 ! micron
  lam_max = maxval(a_mol_using%rad_data%list%lambda)*1.2D0 ! micron
  call make_local_cont_lut(lam_min, lam_max, 300)
  !
end subroutine my_radex_prepare


subroutine my_radex_prepare_molecule
  integer i
  !
  a_mol_using%density_mol = my_radex_config%n_H * my_radex_config%X_x
  a_mol_using%geotype = my_radex_config%geotype
  a_mol_using%Tkin = my_radex_config%Tkin ! K
  a_mol_using%dv = my_radex_config%dv ! cm s-1
  if (my_radex_config%provideLength) then
    a_mol_using%length_scale = my_radex_config%length_scale ! cm
  else
    a_mol_using%length_scale = my_radex_config%Ncol_x / a_mol_using%density_mol
  end if
  !
  ! Set the initial occupation
  a_mol_using%f_occupation = a_mol_using%level_list%weight * &
      exp(-a_mol_using%level_list%energy / &
           a_mol_using%Tkin)
  ! Normalize
  a_mol_using%f_occupation = a_mol_using%f_occupation / &
                             sum(a_mol_using%f_occupation)
  !
  if (my_radex_config%opH2eq3) then
    write(*, '(A)') 'oH2:pH2 = 3:1'
    my_radex_config%X_oH2 = 0.75D0
    my_radex_config%X_pH2 = 0.25D0
  else
    ! TODO
    my_radex_config%X_oH2 = 0.75D0
    my_radex_config%X_pH2 = 0.25D0
  end if
  ! Set the density of the collisional partners
  do i=1, a_mol_using%colli_data%n_partner
    if (a_mol_using%colli_data%list(i)%name_partner .eq. 'H2') then
    !
      a_mol_using%colli_data%list(i)%dens_partner = &
        my_radex_config%n_H * my_radex_config%X_H2
    !
    else if ((a_mol_using%colli_data%list(i)%name_partner .eq. 'o-H2') .or. &
             (a_mol_using%colli_data%list(i)%name_partner .eq. 'oH2')) then
    !
      a_mol_using%colli_data%list(i)%dens_partner = &
        my_radex_config%n_H * my_radex_config%X_oH2
    !
    else if ((a_mol_using%colli_data%list(i)%name_partner .eq. 'p-H2') .or. &
             (a_mol_using%colli_data%list(i)%name_partner .eq. 'pH2')) then
    !
      a_mol_using%colli_data%list(i)%dens_partner = &
        my_radex_config%n_H * my_radex_config%X_pH2
    !
    else if (a_mol_using%colli_data%list(i)%name_partner .eq. 'H') then
    !
      a_mol_using%colli_data%list(i)%dens_partner = &
        my_radex_config%n_H * my_radex_config%X_HI
    !
    else if (a_mol_using%colli_data%list(i)%name_partner .eq. 'H+') then
    !
      a_mol_using%colli_data%list(i)%dens_partner = &
        my_radex_config%n_H * my_radex_config%X_Hplus
    !
    else if (a_mol_using%colli_data%list(i)%name_partner .eq. 'e') then
    !
      a_mol_using%colli_data%list(i)%dens_partner = &
        my_radex_config%n_H * my_radex_config%X_E
    !
    else
    !
      a_mol_using%colli_data%list(i)%dens_partner = 0D0
    !
    end if
    !
    write(*, '("Partner", I2, 2X, A, ": ", ES12.4)') i, &
      a_mol_using%colli_data%list(i)%name_partner, &
      a_mol_using%colli_data%list(i)%dens_partner
  end do
  !
end subroutine my_radex_prepare_molecule


subroutine make_local_cont_lut(lam_min, lam_max, n)
  double precision, intent(in) :: lam_min, lam_max
  integer, intent(in) :: n
  integer i
  double precision dlam, lam, freq
  !
  if (.not. allocated(cont_lut%lam)) then
    cont_lut%n = n
    allocate(cont_lut%lam(cont_lut%n), &
             cont_lut%alpha(cont_lut%n), &
             cont_lut%J(cont_lut%n))
  end if
  !
  dlam = (lam_max - lam_min) / dble(n-1)
  do i=1, cont_lut%n
    cont_lut%lam(i) = lam_min + dlam * dble(i-1)
    cont_lut%alpha(i) = 0D0 ! absorption
    !
    lam = cont_lut%lam(i) + dlam * 0.5D0
    ! Energy per unit area per unit frequency per second per sqradian
    freq = phy_SpeedOfLight_CGS / (lam * phy_micron2cm)
    cont_lut%J(i) = planck_B_nu(my_radex_config%Tcmb, freq)
    !write(*,*) lam, freq
  end do
end subroutine make_local_cont_lut



function planck_B_lambda(T, lambda_CGS)
  double precision planck_B_lambda
  double precision, intent(in) :: T, lambda_CGS
  double precision tmp
  double precision, parameter :: TH = 1D-8
  tmp = (phy_hPlanck_CGS * phy_SpeedOfLight_CGS) / (lambda_CGS * phy_kBoltzmann_CGS * T)
  if (tmp .gt. TH) then
    tmp = exp(tmp) - 1D0
  end if
  planck_B_lambda = &
    2D0*phy_hPlanck_CGS * phy_SpeedOfLight_CGS**2 / lambda_CGS**5 / tmp
end function planck_B_lambda



function planck_B_nu(T, nu)
  double precision planck_B_nu
  double precision, intent(in) :: T, nu
  double precision tmp
  double precision, parameter :: TH = 1D-8
  tmp = (phy_hPlanck_CGS * nu) / (phy_kBoltzmann_CGS * T)
  if (tmp .gt. TH) then
    tmp = exp(tmp) - 1D0
  end if
  planck_B_nu = &
    2D0*phy_hPlanck_CGS * nu**3 / phy_SpeedOfLight_CGS**2 / tmp
end function planck_B_nu




end module my_radex
