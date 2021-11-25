! This file is part of our code to solve the statistical equilibrium problem of
! the energy level population of a uniform volume of gas under the effect of
! background radiation and collisional excitation.
!
! The code works similar to radex
! (see https://sron.rug.nl/~vdtak/radex/index.shtml)
!
! One difference from radex is that we solve the statistical equilibrium
! equation using an ode solver (ODEPACK, see netlib.org/odepack/).
!
! Written by Fujun Du, fujun.du@gmail.com, fdu@umich.edu
!
! 2014-01-02 Thu 02:43:39
!


module global_var
  implicit none
  integer, parameter :: len_filename_max = 256
  double precision g_photon_count_to_J
  double precision, dimension(:,:), allocatable :: data_inputs
  double precision, dimension(:), allocatable :: different_r
  integer len_inputs, n_different_r
  integer nSkip_inputs, ncol_inputs, n_z_step
  character(len=60), parameter :: author_name = 'Fujun Du'
  character(len=60), parameter :: author_email = 'fdu@umich.edu'
  character(len=60), parameter :: user_name = 'Whoever'
  double precision g_inf_len
  double precision dbl_NaN
end module global_var




module phy_const
  implicit none
  double precision, parameter :: phy_Pi = 3.1415926535897932384626433D0
  double precision, parameter :: phy_Pi_2 = 1.57079632679489661923132D0
  double precision, parameter :: phy_2Pi = 6.283185307179586476925D0
  double precision, parameter :: phy_sqrt2Pi = 2.5066282746310005024D0
  double precision, parameter :: phy_GaussFWHM_c = sqrt(phy_Pi/4D0/log(2D0))
  !
  double precision, parameter :: phy_elementaryCharge_SI = 1.602176487D-19
  double precision, parameter :: phy_electronClassicalRadius_SI = 2.8179403267D-15
  double precision, parameter :: phy_electronClassicalRadius_CGS = 2.8179403267D-13
  double precision, parameter :: phy_CoulombConst_SI = 8.9875517873681764D9
  double precision, parameter :: phy_mProton_SI  = 1.67262158D-27 ! kg
  double precision, parameter :: phy_mProton_CGS = 1.67262158D-24 ! g
  double precision, parameter :: phy_mElectron_SI  = 9.10938188E-31 ! kg
  double precision, parameter :: phy_mElectron_CGS = 9.10938188D-28 ! g
  double precision, parameter :: phy_kBoltzmann_SI  = 1.3806503D-23
  double precision, parameter :: phy_kBoltzmann_CGS = 1.3806503D-16
  double precision, parameter :: phy_hPlanck_SI  = 6.62606896D-34
  double precision, parameter :: phy_hPlanck_CGS = 6.62606896D-27
  double precision, parameter :: phy_hbarPlanck_SI  = 1.054571628D-34
  double precision, parameter :: phy_hbarPlanck_CGS = 1.054571628D-27
  double precision, parameter :: phy_GravitationConst_SI  = 6.67428D-11
  double precision, parameter :: phy_GravitationConst_CGS = 6.67428D-8
  double precision, parameter :: phy_SpeedOfLight_SI  = 299792458D0
  double precision, parameter :: phy_SpeedOfLight_CGS = 299792458D2
  double precision, parameter :: phy_StefanBoltzmann_SI = 5.6704D-8
  double precision, parameter :: phy_StefanBoltzmann_CGS = 5.670373D-5
  double precision, parameter :: phy_IdealGasConst_SI = 8.314472D0
  !
  double precision, parameter :: phy_Lsun_SI = 3.839D26 ! J s-1
  double precision, parameter :: phy_Lsun_CGS = 3.839D33 ! erg s-1
  double precision, parameter :: phy_Msun_SI = 1.9891D30 ! kg
  double precision, parameter :: phy_Msun_CGS = 1.9891D33 ! g
  double precision, parameter :: phy_Rsun_SI = 6.955D8 ! m
  double precision, parameter :: phy_Rsun_CGS = 6.955D10 ! cm
  !
  double precision, parameter :: phy_Rearth_CGS = 6371D5 ! cm
  double precision, parameter :: phy_Mearth_CGS = 5.97219D27 ! g
  double precision, parameter :: phy_Mmoon_CGS = 7.34767309D25 ! g
  !
  double precision, parameter :: phy_SecondsPerYear = 3600D0*24D0*365D0
  double precision, parameter :: phy_Deg2Rad = phy_Pi/180D0
  double precision, parameter :: phy_erg2joule = 1D-7
  double precision, parameter :: phy_m2cm = 1D2
  double precision, parameter :: phy_kg2g = 1D3
  double precision, parameter :: phy_eV2erg = 1.60217657D-12
  double precision, parameter :: phy_cm_1_2erg = phy_hPlanck_CGS * phy_SpeedOfLight_CGS
  double precision, parameter :: phy_cm_1_2K = phy_cm_1_2erg/phy_kBoltzmann_CGS
  double precision, parameter :: phy_AvogadroConst = 6.02214179D23
  double precision, parameter :: phy_AU2cm = 1.49597871D13
  double precision, parameter :: phy_AU2m  = 1.49597871D11
  double precision, parameter :: phy_pc2m  = 3.08567758D16
  double precision, parameter :: phy_pc2cm = 3.08567758D18
  double precision, parameter :: phy_Angstrom2micron  = 1D-4
  double precision, parameter :: phy_Angstrom2cm  = 1D-8
  double precision, parameter :: phy_micron2cm  = 1D-4
  !
  double precision, parameter :: phy_jansky2CGS = 1D-23
  double precision, parameter :: phy_jansky2SI  = 1D-26
  !
  double precision, parameter :: phy_CMB_T = 2.72548D0
  !
  double precision, parameter :: phy_ratioDust2GasMass_ISM = 0.01D0
  double precision, parameter :: phy_Habing_photon_energy_CGS = 1.99D-11
  double precision, parameter :: phy_LyAlpha_energy_CGS = 1.64D-11
  double precision, parameter :: phy_UV_cont_energy_CGS = phy_Habing_photon_energy_CGS
  double precision, parameter :: phy_Habing_energy_density_CGS = 5.29D-14 ! Draine 2011 book, equation 12.6
  double precision, parameter :: phy_Habing_photon_flux_CGS = 6D7 ! cm-2 s-1
  double precision, parameter :: phy_Habing_energy_flux_CGS = 1.194D-3 ! erg cm-2 s-1
  double precision, parameter :: phy_UVext2Av = 2.6D0 ! Tielens 2005, eq 3.19
  !
  double precision, parameter :: phy_LyAlpha_nu0 = 2.4660718D15
  double precision, parameter :: phy_LyAlpha_l0 = 1215.668D0
  double precision, parameter :: phy_LyAlpha_dnul = 9.938D7
  double precision, parameter :: phy_LyAlpha_f12 = 0.4162D0
  !
  double precision, parameter :: const_LyAlpha_cross_H2O = 1.2D-17 ! Van Dishoeck 2006, Table 1
  double precision, parameter :: const_LyAlpha_cross_OH  = 1.8D-18 ! Van Dishoeck 2006, Table 1
  !
  double precision, parameter :: const_cosmicray_attenuate_N = 5.75D25 ! 96 g cm-2, Nomura 2007
end module phy_const
