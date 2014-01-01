module disk

use data_struct
use grid
use chemistry
use heating_cooling
use montecarlo
use load_Draine_dusts


implicit none

type :: type_a_dust_component
  integer itype
  double precision pmass_CGS
  type(type_dust_MRN) :: mrn
  type(type_Andrews_disk) :: andrews
end type type_a_dust_component


type :: type_disk_basic_info
  double precision star_luminosity_in_Lsun
  double precision star_mass_in_Msun, star_radius_in_Rsun, star_temperature
  double precision disk_mass_in_Msun
  double precision ratio_uv2total
  double precision ratio_lyman2uv
  double precision ratio_xray2total
  double precision Lyman_phlumi_star_surface, &
                   UV_cont_phlumi_star_surface, &
                   Xray_phlumi_star_surface
  character(len=32) filename_exe
  logical            :: backup_src = .true.
  character(len=128) :: backup_src_cmd = &
    'find *.f90 *.f *.py makefile | cpio -pdm --insecure '
  double precision :: geometric_factor_UV   = 0.01D0
  double precision :: geometric_factor_Xray = 0.001D0
  !
  double precision :: dust2gas_mass_bg = 1D-5
  type(type_Andrews_disk) andrews_gas, andrews_dust, andrews_dust_bg
  !
  integer ndustcompo
  type(type_a_dust_component), dimension(MaxNumOfDustComponents) :: dustcompo
  !
  !double precision :: colDen2Av_coeff = 1D-21 ! Sun Kwok, eq 10.21
  !double precision :: colDen2Av_coeff = 5.3D-22 ! Draine 2011, eq 21.7
end type type_disk_basic_info


type :: type_disk_iter_params
  integer :: n_iter=128, n_iter_used = 0
  integer :: nlocal_iter = 2
  !
  double precision :: rtol_T = 0.1D0,    atol_T = 2D0
  double precision :: rtol_abun = 0.2D0, atol_abun = 1D-12
  !
  logical flag_converged
  integer n_cell_converged
  real converged_cell_percentage_stop
  !
  logical :: redo_montecarlo = .true.
  logical :: flag_save_rates = .false.
  logical :: flag_shortcut_ini = .false.
  !
  integer :: nSpecies_check_refine = 0
  integer :: ncell_refine = 0, count_refine = 0
  integer :: nMax_refine = 2
  double precision :: threshold_ratio_refine = 10D0
  character(len=128) filename_list_check_refine
  !
  character(len=128) iter_files_dir
  !
  logical do_line_transfer
end type type_disk_iter_params


type :: type_simple_integer_list
  integer :: nlen = 0
  integer, dimension(:), allocatable :: vals
end type type_simple_integer_list


type :: type_ana_params
  logical :: do_analyse = .false.
  integer ana_i_incr
  character(len=128) analyse_points_inp_dir, analyse_out_dir
  character(len=128) file_list_analyse_points, file_list_analyse_species, &
                     file_analyse_res_ele, file_analyse_res_contri
  type(type_cell_rz_phy_basic) chempar
end type type_ana_params


type :: type_disk_iter_storage
  double precision, dimension(:), allocatable :: T_s
  double precision, dimension(:,:), allocatable :: abundances
end type type_disk_iter_storage


type :: type_mole_exc_conf
  character(len=128) :: dirname_mol_data=''
  character(len=128) :: fname_mol_data=''
  integer nfreq_window
  double precision, dimension(10) :: freq_mins, freq_maxs
  double precision abundance_factor
  double precision :: E_min = 50D0, E_max = 5D3
  logical :: useLTE = .true.
  !
  integer nf, nth, nx, ny
  double precision dist
  !
end type type_mole_exc_conf


type :: type_molecule_exc
  type(type_mole_exc_conf) :: conf
  type(type_molecule_energy_set), pointer :: p => null()
  integer nlevel_keep, ntran_keep
  integer, dimension(:), allocatable :: ilv_keep, ilv_reverse
  integer, dimension(:), allocatable :: itr_keep, itr_reverse
end type type_molecule_exc


type :: type_book_keeping
  integer fU
  character(len=128) dir, filename_log
end type type_book_keeping


type :: type_image
  integer nx, ny
  double precision xmin, xmax, dx, ymin, ymax, dy
  double precision view_theta
  double precision freq_min, freq_max
  double precision total_flux
  integer iTran
  type(type_rad_transition) rapar
  double precision, dimension(:,:), allocatable :: val
end type type_image


type :: type_cube_header
  integer iTran
  integer nx, ny, nz
  double precision dx, dy, dz
  double precision f0
  double precision theta
  double precision Eup, Elow
  double precision Aul, Bul, Blu
  double precision total_flux_max
end type type_cube_header


type :: type_cube
  type(type_cube_header) :: h
  double precision, dimension(:,:,:), allocatable :: val
end type type_cube


type :: type_fits_par
  character(len=256) :: filename
  integer stat, fU, blocksize, bitpix, naxis
  integer, dimension(3) :: naxes
  integer i, j, group, fpixel, nelements, decimals
  integer pcount, gcount
  logical simple, extend
  character(len=32) :: extname
  character(len=32) :: author, user
end type type_fits_par


! For logging
type(type_book_keeping) a_book_keeping

! Basic config params of the disk model
type(type_disk_basic_info) a_disk

! Initialization params for the disk model
!type(type_disk_basic_info) disk_params_ini

! Iteration params for the disk model
type(type_disk_iter_params) a_disk_iter_params

! Storage for the disk model
type(type_disk_iter_storage) a_iter_stor

! Initialization params for each cell, which are the same for all the cells.
type(type_cell_rz_phy_basic) cell_params_ini

! Params for doing analysis
type(type_ana_params) a_disk_ana_params

! Point (location) list and species list for analysis
type(type_simple_integer_list) :: ana_ptlist, ana_splist

! Columns of cells
type(type_leaves), dimension(:), allocatable :: columns

! Index list of the cells that are being calculated
integer, dimension(:), allocatable :: calculating_cells
integer n_calculating_cells, n_calculating_cells_max

! Filename for saving the results for each iteration; the name will be changed
! for different iterations.
character(len=128) :: filename_save_results
integer fU_save_results

! Energy of a typical X-ray particle, in kev.
! For calculating the X-ray number flux.
! Should find a better way to deal with this.
double precision, parameter, private :: xray_energy_kev = 1D0

! Index of species that are to be used for check whether cell refinement is
! needed.
integer, dimension(:), allocatable, private :: idx_Species_check_refine
double precision, dimension(:), allocatable, private :: &
    thr_Species_check_refine

! For displaying some text to the screen
character(len=256) str_disp

! Field length for certain output
integer, parameter :: len_item=14

type(type_mole_exc_conf) :: mole_line_conf

type(type_molecule_exc) :: mole_exc

! Namelists for reading config params
namelist /disk_configure/ &
  a_disk

namelist /cell_configure/ &
  cell_params_ini

namelist /iteration_configure/ &
  a_disk_iter_params

namelist /analyse_configure/ &
  a_disk_ana_params


namelist /mole_line_configure/ &
  mole_line_conf


contains


!subroutine make_images
!  integer ntr, itr
!  integer i, j, k
!  integer, parameter :: nf = 100, nth = 4
!  integer, parameter :: nx = 101, ny = 101
!  double precision, parameter :: xmin = -5D0, xmax = 5D0, ymin = -5D0, ymax = 5D0
!  double precision, parameter :: dist = 50D0 ! pc
!  double precision VeloHalfWidth
!  double precision delf, df, f0, fmin, dtheta
!  character(len=128) im_dir, fname
!  type(type_image) :: image
!  !
!  im_dir = trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'images/'))
!  if (.not. dir_exist(im_dir)) then
!    call my_mkdir(im_dir)
!  end if
!  !
!  dtheta = 90D0 / dble(nth-1)
!  !
!  ntr = mole_exc%ntran_keep
!  !
!  do i=1, ntr
!    itr = mole_exc%itr_keep(i)
!    f0 = a_mol_using%rad_data%list(itr)%freq
!    do j=1, nf
!      do k=1, nth
!        if (k .eq. 1) then
!          VeloHalfWidth = 20D3
!        else
!          VeloHalfWidth = 100D3
!        end if
!        delf = f0 * VeloHalfWidth / phy_SpeedOfLight_SI ! Todo
!        fmin = f0 - delf
!        df = delf * 2D0 / dble(nf)
!        !
!        image%iTran = itr
!        image%nx = nx
!        image%ny = ny
!        image%xmin = xmin
!        image%xmax = xmax
!        image%ymin = ymin
!        image%ymax = ymax
!        image%dx = (xmax - xmin) / dble(nx-1)
!        image%dy = (ymax - ymin) / dble(ny-1)
!        image%freq_min = fmin + dble(j-1) * df 
!        image%freq_max = fmin + dble(j)   * df 
!        image%view_theta = dtheta * dble(k-1)
!        !
!        image%rapar = a_mol_using%rad_data%list(itr)
!        !
!        write(*, '(3I4, " / ", 3I4)') i, j, k, ntr, nf, nth
!        write(fname, '(3(I0.4,"_"), ES14.4, "_", F9.2, ".dat")') i, j, k, &
!          image%freq_min, image%view_theta
!        call dropout_char(fname, ' ')
!        !
!        call make_a_channel_image(image)
!        !
!        image%total_flux = sum(image%val) * &
!                           (image%dx * image%dy * phy_AU2cm**2 / &
!                            (dist * phy_pc2cm)**2) / &
!                           phy_jansky2CGS
!        !
!        write(*, '(2A)') 'Saving image to ', fname
!        call save_a_image(combine_dir_filename(im_dir, fname), image)
!        !
!      end do
!    end do
!  end do
!end subroutine make_images


!subroutine save_a_image(fname, im)
!  character(len=*) fname
!  type(type_image), intent(in) :: im
!  integer fU
!  integer i, j
!  double precision x, y
!  !
!  call openFileSequentialWrite(fU, fname, 99)
!  write(fU, '(A, 2X, 2I6,     2X, A)') '!', im%nx, im%ny, '= nx, ny'
!  write(fU, '(A, 2X, F10.4,   2X, A)') '!', im%view_theta, '= theta'
!  write(fU, '(A, 2X, ES18.8,  2X, A)') '!', im%rapar%freq, '= frequency'
!  write(fU, '(A, 2X, 2ES18.8, 2X, A)') '!', im%freq_min, im%freq_max, '= freq_min, freq_max'
!  write(fU, '(A, 2X, ES18.8,  2X, A)') '!', im%rapar%lambda, '= wavelength (micron)'
!  write(fU, '(A, 2X, 2F10.2,  2X, A)') '!', im%rapar%Eup, im%rapar%Elow, '= Eup, Elow'
!  write(fU, '(A, 2X, 3ES12.2, 2X, A)') '!', im%rapar%Aul, im%rapar%Bul, im%rapar%Blu, '= Aul, Bul, Blu'
!  write(fU, '(A, 2X, ES12.4,  2X, A)') '!', im%total_flux, '= total flux in jy'
!  do j=1, im%ny
!    do i=1, im%nx
!      x = im%xmin + dble(i-1) * im%dx
!      y = im%ymin + dble(j-1) * im%dy
!      write(fU, '(2ES12.3, ES15.4E3)') x, y, im%val(i, j)
!    end do
!  end do
!  close(fU)
!end subroutine save_a_image




subroutine make_cubes
  use my_timer
  integer ntr, itr
  integer i, j, k, i1, j1
  !integer, parameter :: nf = 100, nth = 4
  !integer, parameter :: nx = 201, ny = 201
  !double precision, parameter :: dist = 50D0 ! pc
  integer nf, nth
  integer nx, ny
  double precision dist
  !
  double precision :: xmin, xmax, ymin, ymax
  double precision VeloHalfWidth
  double precision delf, df, f0, fmin, dtheta
  double precision dv, vmax
  character(len=128) im_dir, fname
  type(type_image) :: image
  type(type_cube) :: cube
  double precision, dimension(:,:), allocatable :: &
    arr_tau, arr_tau1, Ncol_up, Ncol_low
  double precision, dimension(:), allocatable :: vec_flux
  !
  type(type_fits_par) :: fp
  type(date_time) a_date_time
  !
  nf  = mole_exc%conf%nf
  nth = mole_exc%conf%nth
  nx  = mole_exc%conf%nx
  ny  = mole_exc%conf%ny
  dist = mole_exc%conf%dist
  !
  xmax = max(root%xmax, root%ymax)
  xmin = -xmax
  ymin = -xmax
  ymax = xmax
  !
  VeloHalfWidth = 1.3D0 * sqrt( &
    phy_GravitationConst_SI * a_disk%star_mass_in_Msun * phy_Msun_SI / &
      (root%xmin * phy_AU2m) + &
    phy_kBoltzmann_SI * 2D3 / (phy_mProton_SI * 2.8D0) * 3D0)
  !
  fp%stat = 0
  fp%blocksize = 1
  fp%pcount = 0
  fp%gcount = 1
  fp%group=1
  fp%fpixel=1
  fp%decimals = 16
  fp%author = 'Fujun Du (fdu@umich.edu)'
  fp%user = ''
  fp%simple=.true.
  fp%extend=.true.
  fp%bitpix=-64 ! double
  !
  im_dir = trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'images/'))
  if (.not. dir_exist(im_dir)) then
    call my_mkdir(im_dir)
  end if
  !
  dtheta = 90D0 / dble(nth-1)
  !
  ntr = mole_exc%ntran_keep
  !
  image%nx = nx
  image%ny = ny
  image%xmin = xmin
  image%xmax = xmax
  image%ymin = ymin
  image%ymax = ymax
  image%dx = (xmax - xmin) / dble(nx-1)
  image%dy = (ymax - ymin) / dble(ny-1)
  !
  allocate(cube%val(nx, ny, nf), &
           image%val(nx, ny), &
           arr_tau(nx, ny), &
           arr_tau1(nx, ny), &
           Ncol_up(nx, ny), &
           Ncol_low(nx, ny), &
           vec_flux(nf))
  !
  do i=1, ntr ! Transitions
    itr = mole_exc%itr_keep(i)
    image%iTran = itr
    f0 = a_mol_using%rad_data%list(itr)%freq
    delf = f0 * VeloHalfWidth / phy_SpeedOfLight_SI
    fmin = f0 - delf
    df = delf * 2D0 / dble(nf)
    image%rapar = a_mol_using%rad_data%list(itr)
    !
    do k=1, nth ! Viewing angles
      image%view_theta = dtheta * dble(k-1)
      !
      arr_tau = 0D0
      do j=1, nf ! Frequency channels
        !
        image%freq_min = fmin + dble(j-1) * df 
        image%freq_max = fmin + dble(j)   * df 
        !
        write(*, '(3I4, " / ", 3I4)') i, j, k, ntr, nf, nth
        !
        call make_a_channel_image(image, arr_tau1, Ncol_up, Ncol_low, nx, ny)
        !
        do j1=1, ny
        do i1=1, nx
          arr_tau(i1, j1) = max(arr_tau1(i1, j1), arr_tau(i1, j1))
        end do
        end do
        !
        cube%val(:, :, j) = image%val
        !
        image%total_flux = sum(image%val) * &
                           (image%dx * image%dy * phy_AU2cm**2 / &
                            (dist * phy_pc2cm)**2) / &
                           phy_jansky2CGS
        vec_flux(j) = image%total_flux
        !
      end do
      !
      write(fname, '(3(I0.5,"_"), ES14.5, "_", F09.2, ".fits")') i, itr, k, &
        f0, image%view_theta
      call dropout_char(fname, ' ')
      !
      fp%filename = trim(combine_dir_filename(im_dir, fname))
      !
      call ftgiou(fp%fU, fp%stat)
      call ftinit(fp%fU, fp%filename, fp%blocksize, fp%stat)
      !
      fp%naxis=3
      fp%naxes(1)=nx
      fp%naxes(2)=ny
      fp%naxes(3)=nf
      fp%nelements = nx * ny * nf
      call ftphpr(fp%fU, fp%simple, fp%bitpix, fp%naxis, fp%naxes, &
                  fp%pcount, fp%gcount, fp%extend, fp%stat)
      !
      call ftpprd(fp%fU, fp%group, fp%fpixel, fp%nelements, cube%val, fp%stat)
      !
      call ftpkyd(fp%fU, 'CDELT1', image%dx,  fp%decimals, 'dx (AU)', fp%stat)
      call ftpkyd(fp%fU, 'CDELT2', image%dy,  fp%decimals, 'dy (AU)', fp%stat)
      call ftpkyd(fp%fU, 'CDELT3', df      ,  fp%decimals, 'df (Hz)', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX3', 1.0D0,    fp%decimals, 'k0', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL1', xmin,   fp%decimals, 'xmin', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL2', ymin,   fp%decimals, 'ymin', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL3', fmin + 0.5D0 * df,   fp%decimals, 'fmin', fp%stat)
      call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
      call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
      call ftpkys(fp%fU, 'CTYPE3', 'F', 'Hz', fp%stat)
      !
      call ftpkyd(fp%fU, 'Dist',  mole_exc%conf%dist, fp%decimals, 'pc', fp%stat)
      call ftpkyd(fp%fU, 'Theta', image%view_theta, fp%decimals, 'deg', fp%stat)
      call ftpkyj(fp%fU, 'Itr',   itr   ,  'trans num', fp%stat)
      call ftpkyd(fp%fU, 'F0',    f0/1D9,  fp%decimals, 'GHz', fp%stat)
      call ftpkyd(fp%fU, 'lam0',  image%rapar%lambda, fp%decimals, 'micron', fp%stat)
      call ftpkyd(fp%fU, 'Eup',   image%rapar%Eup,  fp%decimals, 'K', fp%stat)
      call ftpkyd(fp%fU, 'Elow',  image%rapar%Elow,  fp%decimals, 'K', fp%stat)
      call ftpkyj(fp%fU, 'iup',   image%rapar%iup,  '', fp%stat)
      call ftpkyj(fp%fU, 'ilow',  image%rapar%ilow, '', fp%stat)
      call ftpkyd(fp%fU, 'Aul',   image%rapar%Aul,  fp%decimals, 's-1', fp%stat)
      call ftpkyd(fp%fU, 'Bul',   image%rapar%Bul,  fp%decimals, '', fp%stat)
      call ftpkyd(fp%fU, 'Blu',   image%rapar%Blu,  fp%decimals, '', fp%stat)
      call ftpkyd(fp%fU, 'MaxFlux', maxval(vec_flux),  fp%decimals, 'jy', fp%stat)
      call ftpkyd(fp%fU, 'MaxTau',  maxval(arr_tau),   fp%decimals, '', fp%stat)
      !
      call ftpkys(fp%fU, 'Author', fp%author, '', fp%stat)
      call ftpkys(fp%fU, 'User',   fp%user,   '', fp%stat)
      call ftpkys(fp%fU, 'SavedAt', trim(a_date_time%date_time_str()), '', fp%stat)
      !
      ! First extension: tau map
      call ftcrhd(fp%fU, fp%stat)
      fp%naxis = 2
      fp%naxes(1) = nx
      fp%naxes(2) = ny
      call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
      !
      call ftpprd(fp%fU, fp%group, fp%fpixel, nx*ny, arr_tau, fp%stat)
      !
      call ftpkyd(fp%fU, 'CDELT1', image%dx,  fp%decimals, 'dx', fp%stat)
      call ftpkyd(fp%fU, 'CDELT2', image%dy,  fp%decimals, 'dy', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL1', xmin,   fp%decimals, 'xmin', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL2', ymin,   fp%decimals, 'ymin', fp%stat)
      call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
      call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
      call ftpkys(fp%fU, 'ExtName', 'TauMap', 'peak values', fp%stat)
      !
      ! Second extension: integrated map
      call ftcrhd(fp%fU, fp%stat)
      fp%naxis = 2
      fp%naxes(1) = nx
      fp%naxes(2) = ny
      call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
      !
      call ftpprd(fp%fU, fp%group, fp%fpixel, nx*ny, sum(cube%val, 3) * df, fp%stat)
      !
      call ftpkyd(fp%fU, 'CDELT1', image%dx,  fp%decimals, 'dx', fp%stat)
      call ftpkyd(fp%fU, 'CDELT2', image%dy,  fp%decimals, 'dy', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL1', xmin,   fp%decimals, 'xmin', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL2', ymin,   fp%decimals, 'ymin', fp%stat)
      call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
      call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
      call ftpkys(fp%fU, 'ExtName', 'IntMap', 'Int(I, nu)', fp%stat)
      !
      ! Third extension: upper column density
      call ftcrhd(fp%fU, fp%stat)
      fp%naxis = 2
      fp%naxes(1) = nx
      fp%naxes(2) = ny
      call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
      !
      call ftpprd(fp%fU, fp%group, fp%fpixel, nx*ny, Ncol_up, fp%stat)
      !
      call ftpkyd(fp%fU, 'CDELT1', image%dx,  fp%decimals, 'dx', fp%stat)
      call ftpkyd(fp%fU, 'CDELT2', image%dy,  fp%decimals, 'dy', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL1', xmin,   fp%decimals, 'xmin', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL2', ymin,   fp%decimals, 'ymin', fp%stat)
      call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
      call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
      call ftpkys(fp%fU, 'ExtName', 'ColumnDensityUp', 'cm-2', fp%stat)
      !
      ! Fourth extension: lower column density
      call ftcrhd(fp%fU, fp%stat)
      fp%naxis = 2
      fp%naxes(1) = nx
      fp%naxes(2) = ny
      call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
      !
      call ftpprd(fp%fU, fp%group, fp%fpixel, nx*ny, Ncol_low, fp%stat)
      !
      call ftpkyd(fp%fU, 'CDELT1', image%dx,  fp%decimals, 'dx', fp%stat)
      call ftpkyd(fp%fU, 'CDELT2', image%dy,  fp%decimals, 'dy', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,    fp%decimals, 'i0', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX2', 1.0D0,    fp%decimals, 'j0', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL1', xmin,   fp%decimals, 'xmin', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL2', ymin,   fp%decimals, 'ymin', fp%stat)
      call ftpkys(fp%fU, 'CTYPE1', 'X', 'AU', fp%stat)
      call ftpkys(fp%fU, 'CTYPE2', 'Y', 'AU', fp%stat)
      call ftpkys(fp%fU, 'ExtName', 'ColumnDensityLow', 'cm-2', fp%stat)
      !
      ! Fifth extension: spectrum integrated over the whole region
      call ftcrhd(fp%fU, fp%stat)
      fp%naxis = 2
      fp%naxes(1) = nf
      fp%naxes(2) = 1
      call ftiimg(fp%fU, fp%bitpix, fp%naxis, fp%naxes(1:2), fp%stat)
      !
      call ftpprd(fp%fU, fp%group, fp%fpixel, nf, vec_flux, fp%stat)
      !
      dv = -df / f0 * phy_SpeedOfLight_SI / 1D3
      vmax = (1D0 - (fmin + 0.5D0 * df)/f0) * phy_SpeedOfLight_SI / 1D3
      call ftpkyd(fp%fU, 'CDELT1', dv   ,  fp%decimals, 'dv', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX1', 1.0D0,  fp%decimals, 'k0', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL1', vmax, fp%decimals, 'vmax', fp%stat)
      call ftpkys(fp%fU, 'CTYPE1', 'V', 'km s-1', fp%stat)
      call ftpkyd(fp%fU, 'CDELT2', 0D0   ,  fp%decimals, 'dumb', fp%stat)
      call ftpkyd(fp%fU, 'CRPIX2', 0.0D0,  fp%decimals, 'dumb', fp%stat)
      call ftpkyd(fp%fU, 'CRVAL2', 0D0, fp%decimals, 'dumb', fp%stat)
      call ftpkys(fp%fU, 'CTYPE2', 'dumb', 'To make ds9 work.', fp%stat)
      call ftpkys(fp%fU, 'ExtName', 'FluxSpec', 'jy', fp%stat)
      !
      call ftclos(fp%fU, fp%stat)
      call ftfiou(fp%fU, fp%stat)
    end do
  end do
end subroutine make_cubes



subroutine make_a_channel_image(im, arr_tau, Ncol_up, Ncol_low, nx, ny)
  integer, intent(in) :: nx, ny
  type(type_image), intent(inout) :: im
  double precision, dimension(nx, ny), intent(out) :: arr_tau, Ncol_up, Ncol_low
  integer i, j, k, i1, j1
  integer xy_sub_div
  integer f_sub_div
  double precision dx, dy, df
  double precision x, y, z, f
  double precision x_ll, y_ll, x_rr, y_rr
  double precision costheta, sintheta
  double precision nave
  double precision tau
  double precision I_0
  double precision Nup, Nlow, Ncol, tau_tot
  type(type_photon_packet) ph
  double precision :: min_tau = 1D-8
  !
  z = -max(root%xmax, root%ymax, &
           abs(im%xmax), abs(im%xmin), &
           abs(im%ymax), abs(im%ymin)) * 5D0
  !
  costheta = cos(im%view_theta / 180D0 * phy_Pi)
  sintheta = sin(im%view_theta / 180D0 * phy_Pi)
  !
  if (.not. allocated(im%val)) then
    allocate(im%val(im%nx, im%ny))
  end if
  !
  a_mol_using => mole_exc%p
  !
  write(*,*)
  do j=1, im%ny
    y_ll = im%ymin + (dble(j-1) - 0.5D0) * im%dy
    y_rr = im%ymin + (dble(j-1) + 0.5D0) * im%dy
    do i=1, im%nx
      x_ll = im%xmin + (dble(i-1) - 0.5D0) * im%dx
      x_rr = im%xmin + (dble(i-1) + 0.5D0) * im%dx
      !
      im%val(i, j)   = 0D0
      arr_tau(i, j)  = 0D0
      Ncol_up(i, j)  = 0D0
      Ncol_low(i, j) = 0D0
      !
      Ncol = max( &
        colden_along_a_direction(x_ll, y_ll, z, &
            costheta, sintheta, a_mol_using%iSpe), &
        colden_along_a_direction(x_ll, y_rr, z, &
            costheta, sintheta, a_mol_using%iSpe), &
        colden_along_a_direction(x_rr, y_ll, z, &
            costheta, sintheta, a_mol_using%iSpe), &
        colden_along_a_direction(x_rr, y_rr, z, &
            costheta, sintheta, a_mol_using%iSpe)) &
        * a_mol_using%abundance_factor
      !
      tau_tot = phy_hPlanck_CGS * im%freq_min / (4D0*phy_Pi) * &
           Ncol / (phy_sqrt2Pi * im%freq_min * a_mol_using%dv &
                   / phy_SpeedOfLight_CGS) &
           * a_mol_using%rad_data%list(im%iTran)%Blu
      if (tau_tot .gt. min_tau) then
        xy_sub_div = 3 + int(10e0/(x_ll*x_ll + y_ll*y_ll + 1.0))
        f_sub_div  = 3 + int(10e0/(x_ll*x_ll + y_ll*y_ll + 1.0))
        nave = dble(f_sub_div * xy_sub_div * xy_sub_div)
      else
        xy_sub_div = 3
        f_sub_div  = 2
        nave = dble(f_sub_div * xy_sub_div * xy_sub_div)
      end if
      !
      dx = im%dx / dble(xy_sub_div-1)
      dy = im%dy / dble(xy_sub_div-1)
      df = (im%freq_max - im%freq_min) / dble(f_sub_div - 1)
      !
      write(*, '(A, 10X, 2I6, " / (", 2I6, ")")') &
        CHAR(27)//'[A', i, j, im%nx, im%ny
      !
      do k=1, f_sub_div
        f = im%freq_min + dble(k - 1) * df
        ph%f = f
        !
        I_0 = planck_B_nu(phy_CMB_T, f)
        !
        ph%ray%vx = 0D0
        ph%ray%vy = -sintheta
        ph%ray%vz =  costheta
        !
        ph%lam = phy_SpeedOfLight_CGS / (f * phy_micron2cm)
        ph%iKap = get_idx_for_kappa(ph%lam, dust_0)
        ph%iTran = im%iTran
        !
        y = y_ll
        do j1=1, xy_sub_div
          x = x_ll
          do i1=1, xy_sub_div
            !
            ph%ray%x =  x
            ph%ray%y =  y * costheta - z * sintheta
            ph%ray%z =  y * sintheta + z * costheta
            !
            ph%Inu = I_0
            !
            call integerate_a_ray(ph, tau, Nup, Nlow)
            !
            if (tau .gt. arr_tau(i, j)) then
              arr_tau(i, j) = tau
            end if
            !
            im%val(i, j) = im%val(i, j) + ph%Inu
            Ncol_up(i, j)  = Ncol_up(i, j) + Nup
            Ncol_low(i, j) = Ncol_low(i, j) + Nlow
            !
            x = x + dx
          end do
          y = y + dy
          !
        end do
      end do
      im%val(i, j)   = im%val(i, j) / nave
      Ncol_up(i, j)  = Ncol_up(i, j) / nave
      Ncol_low(i, j) = Ncol_low(i, j) / nave
    end do
  end do
end subroutine make_a_channel_image



subroutine integerate_a_ray(ph, tau, Nup, Nlow)
  ! ph must be guaranteed to be inside c.
  ! An intersection between ph and c must exist, unless there is a
  ! numerical error.
  type(type_photon_packet), intent(inout) :: ph
  double precision, intent(out) :: tau
  double precision, intent(out) :: Nup, Nlow
  type(type_cell), pointer :: c
  type(type_cell), pointer :: cnext
  logical found
  double precision length, r, z, eps
  integer dirtype
  integer itr, ilow, iup, iL, iU
  double precision ylow, yup
  double precision f0, del_nu, line_alpha, line_J
  double precision tau_this
  double precision t1
  !
  double precision cont_alpha, cont_J
  !
  integer i
  !
  tau = 0D0
  Nup = 0D0
  Nlow = 0D0
  !
  call enter_the_domain_mirror(ph, root, c, found)
  if (.not. found) then
    return
  end if
  !
  do i=1, root%nOffspring*2
    ! Get the intersection between the photon ray and the boundary of the cell
    ! that this photon resides in
    call calc_intersection_ray_cell_mirror(ph%ray, c, &
      length, r, z, eps, found, dirtype)
    if (.not. found) then
      write(*,'(A, I6, 10ES10.2/)') 'In integerate_a_ray, ph not cross c: ', &
        i, &
        ph%ray%x, ph%ray%y, ph%ray%z, &
        ph%ray%vx, ph%ray%vy, ph%ray%vz, &
        c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    !
    if (c%using) then
      !
      call set_using_mole_params(a_mol_using, c)
      !
      itr = ph%iTran
      ilow = a_mol_using%rad_data%list(itr)%ilow
      iup  = a_mol_using%rad_data%list(itr)%iup
      iL = mole_exc%ilv_reverse(ilow)
      iU = mole_exc%ilv_reverse(iup)
      !
      ylow = c%focc%vals(iL)
      yup  = c%focc%vals(iU)
      del_nu = ph%f * a_mol_using%dv / phy_SpeedOfLight_CGS
      f0 = a_mol_using%rad_data%list(itr)%freq
      !
      Nup  = Nup  + a_mol_using%density_mol * length * phy_AU2cm * yup
      Nlow = Nlow + a_mol_using%density_mol * length * phy_AU2cm * ylow
      !
      ! Rybicki & Lightman, p31
      t1 = phy_hPlanck_CGS * f0 / (4D0*phy_Pi) * &
           a_mol_using%density_mol &
           / (phy_sqrt2Pi * del_nu)
      line_alpha = t1 * &
                   (ylow * a_mol_using%rad_data%list(itr)%Blu - &
                    yup  * a_mol_using%rad_data%list(itr)%Bul)
      line_J     = t1 * yup * &
                   a_mol_using%rad_data%list(itr)%Aul
      !
      if ((ph%iKap .gt. 0) .and. c%using) then
        call make_local_cont_lut(c)
        cont_alpha = cont_lut%alpha(ph%iKap)
        cont_J = cont_lut%J(ph%iKap) * cont_alpha
      else
        cont_alpha = 0D0
        cont_J = 0D0
      end if
      call integrate_within_one_cell(ph, length, f0, del_nu, &
        line_alpha, line_J, cont_alpha, cont_J, tau_this)
      tau = tau + tau_this
    end if
    !
    ph%ray%x = ph%ray%x + ph%ray%vx * (length + eps)
    ph%ray%y = ph%ray%y + ph%ray%vy * (length + eps)
    ph%ray%z = ph%ray%z + ph%ray%vz * (length + eps)
    !
    call locate_photon_cell_mirror(r, z, c, cnext, found)
    if (.not. found) then! Not entering a neighboring cell
      ! May be entering a non-neighboring cell?
      call enter_the_domain_mirror(ph, root, cnext, found)
      if (.not. found) then ! Escape
        return
      end if
    end if
    !
    c => cnext
    !
  end do
  !
  write(*, '(A)') 'In integerate_a_ray:'
  write(*, '(A)') 'Should not reach here!'
  write(*,'(I6, 10ES10.2/)') &
    i, &
    ph%ray%x, ph%ray%y, ph%ray%z, &
    ph%ray%vx, ph%ray%vy, ph%ray%vz, &
    c%xmin, c%xmax, c%ymin, c%ymax
  write(*,'(4ES10.2, I4, L4/)') r, z, length, eps, dirtype, found
  stop
end subroutine integerate_a_ray


function colden_along_a_direction(x, y, z, costheta, sintheta, iSpe) result(N)
  double precision N
  double precision, intent(in) :: x, y, z
  double precision, intent(in) :: costheta, sintheta
  integer, intent(in) :: iSpe
  type(type_photon_packet) :: ph
  !
  ph%ray%vx = 0D0
  ph%ray%vy = -sintheta
  ph%ray%vz =  costheta
  !
  ph%ray%x =  x
  ph%ray%y =  y * costheta - z * sintheta
  ph%ray%z =  y * sintheta + z * costheta
  !
  N = colden_along_a_ray(ph, iSpe)
end function colden_along_a_direction



function colden_along_a_ray(ph0, iSpe) result(N)
  ! ph must be guaranteed to be inside c.
  ! An intersection between ph and c must exist, unless there is a
  ! numerical error.
  double precision N
  type(type_photon_packet), intent(in) :: ph0
  integer, intent(in) :: iSpe
  type(type_cell), pointer :: c
  type(type_cell), pointer :: cnext
  type(type_photon_packet) :: ph
  double precision length, r, z, eps
  logical found
  integer dirtype
  !
  integer i
  !
  N = 0D0
  !
  ph = ph0
  !
  call enter_the_domain_mirror(ph, root, c, found)
  if (.not. found) then
    return
  end if
  !
  do i=1, root%nOffspring*2
    ! Get the intersection between the photon ray and the boundary of the cell
    ! that this photon resides in
    call calc_intersection_ray_cell_mirror(ph%ray, c, &
      length, r, z, eps, found, dirtype)
    if (.not. found) then
      write(*,'(A, I6, 10ES10.2/)') 'In colden_along_a_ray, ph not cross c: ', &
        i, &
        ph%ray%x, ph%ray%y, ph%ray%z, &
        ph%ray%vx, ph%ray%vy, ph%ray%vz, &
        c%xmin, c%xmax, c%ymin, c%ymax
      return
    end if
    !
    if (c%using) then
      !
      N = N + c%par%n_gas * c%abundances(iSpe) * length * phy_AU2cm
      !
    end if
    !
    ph%ray%x = ph%ray%x + ph%ray%vx * (length + eps)
    ph%ray%y = ph%ray%y + ph%ray%vy * (length + eps)
    ph%ray%z = ph%ray%z + ph%ray%vz * (length + eps)
    !
    call locate_photon_cell_mirror(r, z, c, cnext, found)
    if (.not. found) then! Not entering a neighboring cell
      ! May be entering a non-neighboring cell?
      call enter_the_domain_mirror(ph, root, cnext, found)
      if (.not. found) then ! Escape
        return
      end if
    end if
    !
    c => cnext
    !
  end do
  !
  write(*, '(A)') 'In colden_along_a_ray:'
  write(*, '(A)') 'Should not reach here!'
  write(*,'(I6, 10ES10.2/)') &
    i, &
    ph%ray%x, ph%ray%y, ph%ray%z, &
    ph%ray%vx, ph%ray%vy, ph%ray%vz, &
    c%xmin, c%xmax, c%ymin, c%ymax
  write(*,'(4ES10.2, I4, L4/)') r, z, length, eps, dirtype, found
  stop
end function colden_along_a_ray



subroutine integrate_within_one_cell(ph, length, f0, del_nu, &
        line_alpha, line_J, cont_alpha, cont_J, tau)
  type(type_photon_packet), intent(inout) :: ph
  double precision, intent(in) :: length, f0, del_nu, line_alpha, line_J, cont_alpha, cont_J
  double precision nu, nu1, nu2, dnu, x, dl, dtau, t1
  double precision, intent(out) :: tau
  double precision jnu, knu
  type(type_ray) ray
  integer i, ndiv
  !
  ray%x = ph%ray%x + ph%ray%vx * length
  ray%y = ph%ray%y + ph%ray%vy * length
  ray%z = ph%ray%z + ph%ray%vz * length
  ray%vx = ph%ray%vx
  ray%vy = ph%ray%vy
  ray%vz = ph%ray%vz
  !
  nu1 = get_doppler_nu(star_0%mass, ph%f, ph%ray)
  nu2 = get_doppler_nu(star_0%mass, ph%f, ray)
  !
  ndiv = 1 + &
    min(int(10D0 * abs(nu1 - nu2) / del_nu), &
        int(1D2 * (line_alpha + cont_alpha) * length * phy_AU2cm))
  dnu = (nu2 - nu1) / dble(ndiv)
  dl = length * phy_AU2cm / dble(ndiv)
  nu = nu1
  !
  tau = 0D0
  !
  do i=1, ndiv
    x = (nu - f0) / del_nu
    !
    if ((x .gt. 20D0) .or. (x .lt. -20D0)) then
      t1 = 0D0
    else
      t1 = exp(-x*x*0.5D0)
    end if
    !
    jnu = t1 * line_J + cont_J
    knu = t1 * line_alpha + cont_alpha
    !
    dtau = knu * dl
    tau = tau + dtau
    !
    if (dtau .ge. 1D-4) then
      if (dtau .gt. 100D0) then
        t1 = 0D0
        ph%Inu = jnu/knu
      else
        t1 = exp(-dtau)
        ph%Inu = ph%Inu * t1 + jnu/knu * (1D0 - t1)
      end if
    else
      ph%Inu = ph%Inu * (1D0 - dtau) + dl * jnu
    end if
    nu = nu + dnu
  end do
end subroutine integrate_within_one_cell




subroutine line_excitation_do
  type(type_cell), pointer :: c
  integer i
  if (.not. a_disk_iter_params%do_line_transfer) then
    return
  end if
  write(*, '(/A/)') 'Doing energy level excitation calculation.'
  do i=1, leaves%nlen
    c => leaves%list(i)%p
    call do_exc_calc(c)
    write(*, '(I4, 4ES10.2)') i, c%xmin, c%xmax, c%ymin, c%ymax
  end do
end subroutine line_excitation_do


subroutine line_tran_prep
  if (a_disk_iter_params%do_line_transfer) then
    call load_exc_molecule
  end if
  statistic_equil_params%NEQ = mole_exc%p%n_level
  statistic_equil_params%LIW = 20 + statistic_equil_params%NEQ
  statistic_equil_params%LRW = 22 + 9*statistic_equil_params%NEQ + &
                               statistic_equil_params%NEQ*statistic_equil_params%NEQ
  if (statistic_equil_params%NEQ .gt. mole_exc%p%n_level) then
    if (allocated(statistic_equil_params%IWORK)) then
      deallocate(statistic_equil_params%IWORK, statistic_equil_params%RWORK)
    end if
  end if
  if (.not. allocated(statistic_equil_params%IWORK)) then
    allocate(statistic_equil_params%IWORK(statistic_equil_params%LIW), &
             statistic_equil_params%RWORK(statistic_equil_params%LRW))
  end if
end subroutine line_tran_prep



subroutine load_exc_molecule
  integer i, i0, i1, j
  character(len=const_len_species_name) str, str1
  integer, dimension(:), allocatable :: itmp, itmp1
  double precision freq, en
  integer iup, ilow
  logical in_freq_window
  !
  mole_exc%conf = mole_line_conf
  allocate(mole_exc%p)
  !
  a_mol_using => mole_exc%p
  !
  mole_exc%p%abundance_factor = mole_exc%conf%abundance_factor
  !
  call load_moldata_LAMBDA(&
    combine_dir_filename(mole_exc%conf%dirname_mol_data, &
    mole_exc%conf%fname_mol_data))
  !
  a_mol_using%iType = -1
  !
  i = index(a_mol_using%name_molecule, '(')
  if (i .eq. 0) then
    str = a_mol_using%name_molecule
    a_mol_using%iType = 0
    str1 = ''
  else
    str = a_mol_using%name_molecule(1:(i-1))
    i = index(a_mol_using%name_molecule, 'ortho')
    if (i .ne. 0) then
      a_mol_using%iType = 1
      str1 = 'ortho'
    else
      i = index(a_mol_using%name_molecule, 'para')
      if (i .ne. 0) then
        a_mol_using%iType = 2
        str1 = 'para'
      end if
    end if
  end if
  !
  a_mol_using%iSpe = -1
  !
  do i=1, chem_species%nSpecies
    if (str .eq. chem_species%names(i)) then
      a_mol_using%iSpe = i
      exit
    end if
  end do
  if ((a_mol_using%iSpe .eq. -1) .or. (a_mol_using%iType .eq. -1)) then
    write(*, '(A)') 'In load_exc_molecule:'
    write(*, '(A)') 'Unidentified molecule name and/or type:'
    write(*, '(A)') a_mol_using%name_molecule
    write(*, '(A)') 'In file:'
    write(*, '(A)') combine_dir_filename( &
      mole_exc%conf%dirname_mol_data, &
      mole_exc%conf%fname_mol_data)
    stop
  end if
  write(*, '(A, 2A16)') 'Molecule: ', trim(str), str1
  write(*, '(A, I6)') 'Total number of levels: ', mole_exc%p%n_level
  write(*, '(A, I6)') 'Total number of radiative transitions: ', mole_exc%p%rad_data%n_transition
  write(*, '(A, I6)') 'Total number of collisional partners: ', mole_exc%p%colli_data%n_partner
  do i=1, mole_exc%p%colli_data%n_partner
    write(*, '(I2, 2X, 2A)') i, 'Partner name: ', mole_exc%p%colli_data%list(i)%name_partner
    write(*, '(I2, 2X, A, I6)') i, 'Total number of collisional transitions: ', &
      mole_exc%p%colli_data%list(i)%n_transition
    write(*, '(I2, 2X, A, I6)') i, 'Total number of collisional temperatures: ', &
      mole_exc%p%colli_data%list(i)%n_T
  end do
  !write(*,*)
  !write(*, '(A, 2ES12.4)') 'Frequency range to consider: ', &
  !     mole_exc%conf%freq_min, mole_exc%conf%freq_max
  !
  allocate(itmp(a_mol_using%n_level), &
           itmp1(a_mol_using%rad_data%n_transition), &
           mole_exc%ilv_reverse(a_mol_using%n_level))
  mole_exc%ilv_reverse = 0
  i0 = 0
  i1 = 0
  do i=1, a_mol_using%rad_data%n_transition
    freq = a_mol_using%rad_data%list(i)%freq
    en   = a_mol_using%rad_data%list(i)%Eup
    !
    in_freq_window = .false.
    do j=1, mole_exc%conf%nfreq_window
      if ((mole_exc%conf%freq_mins(j) .le. freq) .and. &
          (freq .le. mole_exc%conf%freq_maxs(j))) then
        in_freq_window = .true.
        exit
      end if
    end do
    !
    if (in_freq_window .and. &
        (en .ge. mole_exc%conf%E_min) .and. &
        (en .le. mole_exc%conf%E_max)) then
      i1 = i1 + 1
      itmp1(i1) = i
      iup = a_mol_using%rad_data%list(i)%iup
      ilow = a_mol_using%rad_data%list(i)%ilow
      if (.not. is_in_list_int(ilow, i0, itmp(1:i0))) then
        i0 = i0 + 1
        itmp(i0) = ilow
        mole_exc%ilv_reverse(ilow) = i0
      end if
      if (.not. is_in_list_int(iup, i0, itmp(1:i0))) then
        i0 = i0 + 1
        itmp(i0) = iup
        mole_exc%ilv_reverse(iup) = i0
      end if
    end if
  end do
  !
  mole_exc%nlevel_keep = i0
  mole_exc%ntran_keep  = i1
  allocate(mole_exc%ilv_keep(i0), &
           mole_exc%itr_keep(i1))
  mole_exc%ilv_keep = itmp(1:i0)
  mole_exc%itr_keep = itmp1(1:i1)
  deallocate(itmp, itmp1)
  write(*, '(A, I6)') 'Number of levels to keep:', mole_exc%nlevel_keep
  write(*, '(A, I6)') 'Number of transitions to keep:', mole_exc%ntran_keep
  !do i=1, mole_exc%nlevel_keep
  !  i0 = mole_exc%ilv_keep(i)
  !  write(*, '(I4, ES12.4)') i, mole_exc%p%level_list(i0)%energy
  !end do
  write(*,*)
end subroutine load_exc_molecule


subroutine set_using_mole_params(mole, c)
  type(type_cell), intent(in), pointer :: c
  type(type_molecule_energy_set), intent(inout), pointer :: mole
  select case (mole%iType)
  case (0)
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe)
  case (1)
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe) * 0.75D0
  case (2)
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe) * 0.25D0
  case default
    write(*, '(A)') 'In set_using_mole_params:'
    write(*, '(A, I4)') 'Unknown molecule type: ', mole%iType
    write(*, '(A)') 'Will use the full abundance.'
    mole%density_mol = c%par%n_gas * c%abundances(mole%iSpe)
  end select
  !
  mole%density_mol = mole%density_mol * mole%abundance_factor
  !
  mole%Tkin = c%par%Tgas
  mole%dv = c%par%velo_width_turb
  mole%length_scale = c%par%coherent_length
end subroutine set_using_mole_params


subroutine do_exc_calc(c)
  type(type_cell), intent(inout), pointer :: c
  integer i
  !
  a_mol_using => mole_exc%p
  !
  call set_using_mole_params(a_mol_using, c)
  !
  a_mol_using%f_occupation = a_mol_using%level_list%weight * &
      exp(-a_mol_using%level_list%energy / a_mol_using%Tkin)
  a_mol_using%f_occupation = a_mol_using%f_occupation / sum(a_mol_using%f_occupation)
  !
  if (.not. mole_exc%conf%useLTE) then
    do i=1, a_mol_using%colli_data%n_partner
      select case (a_mol_using%colli_data%list(i)%name_partner)
      case ('H2')
        a_mol_using%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_H2
      case ('o-H2')
        a_mol_using%colli_data%list(i)%dens_partner = &
          0.75D0 * c%par%n_gas * c%par%X_H2
      case ('p-H2')
        a_mol_using%colli_data%list(i)%dens_partner = &
          0.25D0 * c%par%n_gas * c%par%X_H2
      case ('H')
        a_mol_using%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_HI
      case ('H+')
        a_mol_using%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_Hplus
      case ('e')
        a_mol_using%colli_data%list(i)%dens_partner = &
          c%par%n_gas * c%par%X_E
      case default
        write(*, '(A)') 'In do_exc_calc:'
        write(*, '(A)') 'Unknown collision partner:'
        write(*, '(A)') a_mol_using%colli_data%list(i)%name_partner
        write(*, '(A)') 'Will use zero abundance for this partner.'
        a_mol_using%colli_data%list(i)%dens_partner = 0D0
      end select
    end do
    !
    call make_local_cont_lut(c)
    !
    call statistic_equil_solve
  end if
  !
  if (.not. allocated(c%focc)) then
    allocate(c%focc)
    c%focc%nlevels = mole_exc%nlevel_keep
    allocate(c%focc%vals(c%focc%nlevels))
  end if
  c%focc%vals = a_mol_using%f_occupation(mole_exc%ilv_keep)
  !
  !nullify(a_mol_using)
end subroutine do_exc_calc



subroutine make_local_cont_lut(c)
  type(type_cell), intent(in), pointer :: c
  integer i
  double precision dlam, lam
  !
  if (.not. allocated(cont_lut%lam)) then
    cont_lut%n = dust_0%n
    allocate(cont_lut%lam(dust_0%n), &
             cont_lut%alpha(dust_0%n), &
             cont_lut%J(dust_0%n))
  end if
  !
  do i=1, cont_lut%n
    cont_lut%lam(i) = dust_0%lam(i)
    cont_lut%alpha(i) = c%optical%summed(i)
    !
    if (i .lt. cont_lut%n) then
      dlam = cont_lut%lam(i+1) - cont_lut%lam(i)
      lam = (cont_lut%lam(i+1) + cont_lut%lam(i)) * 0.5D0
      ! Energy per unit area per unit frequency per second per sqradian
      cont_lut%J(i) = c%optical%flux(i) &
        / dlam * lam * lam * phy_micron2cm / phy_SpeedOfLight_CGS &
        / (4D0 * phy_Pi)
    else
      cont_lut%J(i) = cont_lut%J(i-1)
    end if
  end do
end subroutine make_local_cont_lut



subroutine montecarlo_prep
  integer n, n0
  double precision lam_max, lam_start
  double precision, dimension(:), allocatable :: ltmp, vtmp
  type(type_stellar_spectrum) stmp
  double precision, parameter :: T_Lya = 1000D0
  !
  mc_conf%minw = sin(mc_conf%min_ang*phy_Deg2Rad)
  mc_conf%maxw = sin(mc_conf%max_ang*phy_Deg2Rad)
  mc_conf%maxw = get_surf_max_angle()
  write(*,'(/A, 2ES12.4)') 'minw,maxw = ', mc_conf%minw, mc_conf%maxw
  !
  dust_0 = dusts%list(1)
  !
  call load_H2O_ab_crosssection( &
    combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_water), &
    water_0)
  !
  call make_H_Lya(T_Lya, HI_0)
  !
  ! Let the dust and H2O, H data share the same wavelength axis
  call align_optical_data
  !
  !call make_LUT_Tdust(dust_0, lut_0)
  !
  ! Create the lookup table for finding temperature and new wavelength during
  ! the Monte Carlo
  call make_luts
  !
  ! Prepare for the stellar spectrum
  lam_max = min(1D6, dust_0%lam(dust_0%n)) ! in angstrom
  if (mc_conf%use_blackbody_star) then
    call make_stellar_spectrum(dust_0%lam(1), &
      lam_max, 10000, star_0)
  else
    !
    call load_stellar_spectrum( &
      trim(combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_star)), &
      star_0)
    if (star_0%lam(star_0%n) .lt. lam_max) then
      ! Fill the remainder with blackbody radiation
      stmp%mass = star_0%mass
      stmp%radius = star_0%radius
      stmp%T = star_0%T
      lam_start = 2D0 * star_0%lam(star_0%n) - star_0%lam(star_0%n - 1)
      n = max(10, int(lam_max / lam_start * 10))
      call make_stellar_spectrum(lam_start, &
        lam_max, n, stmp)
      !
      n0 = star_0%n
      allocate(ltmp(n0), vtmp(n0))
      ltmp = star_0%lam
      vtmp = star_0%vals
      !
      star_0%n = star_0%n + n
      deallocate(star_0%lam, star_0%vals)
      allocate(star_0%lam(star_0%n), star_0%vals(star_0%n))
      !
      star_0%lam(1:n0) = ltmp
      star_0%vals(1:n0) = vtmp
      star_0%lam((n0+1) : star_0%n) = stmp%lam
      star_0%vals((n0+1): star_0%n) = stmp%vals
      !
      deallocate(ltmp, vtmp, stmp%lam, stmp%vals)
    end if
    !
  end if
  !
  star_0%lumi = get_stellar_luminosity(star_0)
  star_0%lumi_UV = get_stellar_luminosity(star_0, lam_range_UV(1), &
    lam_range_UV(2))
  write(*,'(A, ES16.6, A)') 'Stellar total luminosity: ', star_0%lumi, &
    ' erg s-1.'
  write(*,'(A, ES16.6, A)') 'Stellar UV luminosity: ', star_0%lumi_UV, &
    ' erg s-1.'
  !
  allocate(star_0%vals0(star_0%n))
  star_0%vals0 = star_0%vals
  star_0%lumi0 = star_0%lumi
  star_0%lumi_UV0 =  star_0%lumi_UV
  !
  call get_mc_stellar_par(mc_conf)
  !
  ! Global optical property collection
  call make_global_coll
  !
  p4lam%n = luts%list(1)%m
  allocate(p4lam%pvals(0:p4lam%n))
  !
end subroutine montecarlo_prep



subroutine make_dusts_data
  integer i, j, itype, nradius, nradius_prev, nlam
  double precision rmin, rmax, ind, swei, m
  double precision, dimension(:), allocatable :: t1, t2
  !
  dusts%n = a_disk%ndustcompo
  allocate(dusts%list(dusts%n))
  nradius = 0
  !
  do i=1, dusts%n
    !
    call calc_dust_MRN_par(a_disk%dustcompo(i)%mrn)
    !
    rmin  = a_disk%dustcompo(i)%mrn%rmin
    rmax  = a_disk%dustcompo(i)%mrn%rmax
    ind   = a_disk%dustcompo(i)%mrn%n
    itype = a_disk%dustcompo(i)%itype
    !
    rmax = max(rmax, rmin*1.0001D0)
    !
    nradius_prev = nradius
    nradius = dustmix_data%list(itype)%nradius
    nlam = dustmix_data%list(itype)%nlam
    !
    if (i .gt. 1) then
      if (dusts%list(i-1)%n .ne. nlam) then
        write(*,'(A)') 'In make_dusts_data:'
        write(*,'(A)') 'Arrays for different dust types not'// &
                       'having the same dimension!'
        stop
      end if
      if (nradius .ne. nradius_prev) then
        write(*, '(A)') 'In make_dusts_data:'
        write(*, '(A)') 'Inconsistent radius array size!'
        stop
      end if
    end if
    !
    dusts%list(i)%n = nlam
    !
    allocate(dusts%list(i)%lam(nlam), &
             dusts%list(i)%ab(nlam), &
             dusts%list(i)%sc(nlam), &
             dusts%list(i)%g(nlam))
    !
    if (.not. allocated(t1)) then
      allocate(t1(nradius), t2(nradius))
    end if
    !
    do j=1, nlam
      !
      dusts%list(i)%lam(j) = &
        dustmix_data%list(itype)%w(j) / phy_Angstrom2micron
      !
      t1 = exp(-ind*log(dustmix_data%list(itype)%r)) ! = r**(-ind)
      !
      swei = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t1, rmin, rmax)
      !
      t2 = t1 * dustmix_data%list(itype)%ab(j, :)
      dusts%list(i)%ab(j) = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t2, rmin, rmax)
      !
      t2 = t1 * dustmix_data%list(itype)%sc(j, :)
      dusts%list(i)%sc(j) = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t2, rmin, rmax)
      !
      t2 = t1 * dustmix_data%list(itype)%g(j, :)
      dusts%list(i)%g(j) = discrete_integral( &
        nradius, dustmix_data%list(itype)%r, t2, rmin, rmax)
      !
      m = 4D0*phy_Pi/3D0 * a_disk%dustcompo(i)%mrn%r3av * &
          phy_micron2cm**3 * dustmix_info%mix(itype)%rho
      a_disk%dustcompo(i)%pmass_CGS = m ! dust particle mass in gram
      !
      ! Now the unit of ab and sc become cm2 g-1.
      dusts%list(i)%ab(j) = dusts%list(i)%ab(j) / swei * phy_micron2cm**2 / m
      dusts%list(i)%sc(j) = dusts%list(i)%sc(j) / swei * phy_micron2cm**2 / m
      dusts%list(i)%g(j) = dusts%list(i)%g(j) / swei
      !write(*, '(2I4, 10ES12.4)') i, j, dusts%list(i)%lam(j), &
      !  dusts%list(i)%ab(j), dusts%list(i)%sc(j), dusts%list(i)%g(j), &
      !  m, a_disk%dustcompo(i)%mrn%r3av, swei, rmin, rmax, maxval(dustmix_data%list(itype)%g(j, :))
    end do
  end do
end subroutine make_dusts_data



subroutine disk_iteration
  use my_timer
  type(date_time) a_date_time
  integer i, i0, i_count, l_count, ii
  !
  call disk_iteration_prepare
  !
  call montecarlo_prep
  !
  call line_tran_prep
  !
  ! call montecarlo_reset_cells
  ! !
  ! call montecarlo_do(mc_conf, root)
  ! !
  ! call post_montecarlo
  ! !
  ! call openFileSequentialWrite(ii, &
  !   combine_dir_filename(mc_conf%mc_dir_out, 'cmp_radmc.dat'), 999)
  ! do i=1, leaves%nlen
  !   write(ii, '(2F9.2, 2I11, 18ES14.6)') &
  !     leaves%list(i)%p%par%Tdust1, leaves%list(i)%p%par%Tdust, &
  !     leaves%list(i)%p%optical%ab_count_dust, &
  !     leaves%list(i)%p%optical%cr_count, &
  !     leaves%list(i)%p%optical%en_gain_dust, &
  !     leaves%list(i)%p%optical%en_gain_abso, &
  !     leaves%list(i)%p%par%n_gas, &
  !     leaves%list(i)%p%par%mdust_tot, &
  !     leaves%list(i)%p%par%Ncol_toISM, &
  !     leaves%list(i)%p%par%Ncol_toStar, &
  !     leaves%list(i)%p%par%flux_UV, &
  !     leaves%list(i)%p%par%flux_Lya, &
  !     leaves%list(i)%p%par%dir_UV_r, &
  !     leaves%list(i)%p%par%dir_UV_z, &
  !     leaves%list(i)%p%par%aniso_UV, &
  !     leaves%list(i)%p%par%dir_Lya_r, &
  !     leaves%list(i)%p%par%dir_Lya_z, &
  !     leaves%list(i)%p%par%aniso_Lya, &
  !     leaves%list(i)%p%par%rmin, &
  !     leaves%list(i)%p%par%rmax, &
  !     leaves%list(i)%p%par%zmin, &
  !     leaves%list(i)%p%par%zmax
  ! end do
  ! close(ii)
  ! return
  !
  call save_post_config_params
  !
  ! Now start the major big loop.
  !
  a_disk_iter_params%count_refine = 0
  !
  write(str_disp, '("! ", A)') "Iteration begins."
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  do ii = 1, a_disk_iter_params%n_iter
    !
    a_disk_iter_params%n_iter_used = ii
    !
    write(str_disp, '("! ", A)') "Monte Carlo begins."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    call montecarlo_reset_cells
    !
    call montecarlo_do(mc_conf, root)
    !
    ! Retrieve physical parameters from the monte carlo results
    call post_montecarlo
    !
    write(str_disp, '("! ", A)') "Monte Carlo finished."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    ! Write header to the file
    call disk_save_results_pre
    !
    ! Calculate layer by layer.
    ! Start from the surface layer.
    n_calculating_cells = surf_cells%nlen
    calculating_cells(1:surf_cells%nlen) = surf_cells%idx
    i_count = 0 ! Counter for cells
    l_count = 0 ! Counter for layers
    do
      l_count = l_count + 1
      !
      do i=1, n_calculating_cells
        i_count = i_count + 1
        i0 = calculating_cells(i)
        !
        write(*, '(3(A, I5, A, I5, ",", 2X), (A, I4, ","), 2X, A, 4F8.3)') &
          "Iter:", a_disk_iter_params%n_iter_used, "/", a_disk_iter_params%n_iter, &
          "Cell:", i_count, '/', leaves%nlen, &
          "cell:", i, '/', n_calculating_cells, &
          "Layer:", l_count, &
          'rz:', &
          leaves%list(i0)%p%par%rmin, &
          leaves%list(i0)%p%par%rmax, &
          leaves%list(i0)%p%par%zmin, &
          leaves%list(i0)%p%par%zmax
        write(*, '(2(A, ES10.3, 2X), 2X, 2A, 2X, 2A)') &
          'n_gas: ', leaves%list(i0)%p%par%n_gas, &
          'Tdust: ', leaves%list(i0)%p%par%Tdust, &
          'exe: ', trim(a_disk%filename_exe), &
          'dir: ', trim(a_disk_iter_params%iter_files_dir)
        !
        call calc_this_cell(i0)
        !
        call check_convergency_cell(i0)
        !
        write(*, '(12X, 10A10)') chem_idx_some_spe%names(1:10)
        write(*, '(A, 2X, 10ES10.3, L3/)') 'Abundances:',  &
          leaves%list(i0)%p%abundances(chem_idx_some_spe%idx(1:10)), &
          leaves%list(i0)%p%converged
        !
        a_iter_stor%T_s(i0) = leaves%list(i0)%p%par%Tgas
        !
        a_iter_stor%abundances(:, i0) = &
          leaves%list(i0)%p%abundances(chem_idx_some_spe%idx)
        !
        call disk_save_results_write(fU_save_results, leaves%list(i0)%p)
        flush(fU_save_results)
      end do
      !
      call update_calculating_cells
      !
      if (n_calculating_cells .eq. 0) then
        exit
      end if
    end do
    !
    write(str_disp, '("! ", A, I4, A)') "Iteration ", ii, " finished."
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
    call display_string_both(str_disp, a_book_keeping%fU)
    !
    ! At this point all the layers have been walked through.
    call check_convergency_whole_disk
    !
    write(fU_save_results, '(A, L6)') &
      '! flag_converged = ', a_disk_iter_params%flag_converged
    write(fU_save_results, '(A)') &
      '! Finish saving ' // trim(filename_save_results)
    write(fU_save_results, '(A)') '! at ' // trim(a_date_time%date_time_str())
    flush(fU_save_results)
    close(fU_save_results)
    !
    if (a_disk_iter_params%flag_converged) then
      ! Converged.  Finish iteration.
      exit
    !
    else if (a_disk_iter_params%redo_montecarlo) then
      !if (ii .ge. min(4, a_disk_iter_params%n_iter/3)) then
      if (mod(ii, 7) .eq. 6) then
        ! Adjust the vertical structure
        call vertical_pressure_gravity_balance
        write(str_disp, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
        call display_string_both(str_disp, a_book_keeping%fU)
        !
      end if
      cycle
    !
    else if (a_disk_iter_params%count_refine .gt. &
             a_disk_iter_params%nMax_refine) then
      write(str_disp, '(A, I4, " > ", I4)') &
        '! Will not refine any more. count_refine: ', &
        a_disk_iter_params%count_refine, a_disk_iter_params%nMax_refine
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    !
    else
      write(*, '(/A)') 'Doing refinements where necessary.'
      !
      call do_refine
      !
      if (a_disk_iter_params%ncell_refine .ge. 1) then
        !
        a_disk_iter_params%count_refine = a_disk_iter_params%count_refine + 1
        a_disk_iter_params%flag_converged = .false.
        !
        write(*, '(I5, " out of ", I5, " cells are refined.", /)') &
          a_disk_iter_params%ncell_refine, leaves%nlen
        !
        call remake_index
        !
        call load_ana_points_list ! Reload, actually
        !
        if (allocated(a_iter_stor%T_s)) then
          deallocate(a_iter_stor%T_s, a_iter_stor%abundances)
        end if
        allocate(a_iter_stor%T_s(leaves%nlen), &
                 a_iter_stor%abundances(chem_idx_some_spe%nItem, &
                                                leaves%nlen))
        do i=1, leaves%nlen
          a_iter_stor%T_s(i) = leaves%list(i)%p%par%Tgas
          a_iter_stor%abundances(:,i) = &
            leaves%list(i)%p%abundances(chem_idx_some_spe%idx)
        end do
        !
        if (allocated(calculating_cells)) then
          deallocate(calculating_cells)
        end if
        n_calculating_cells_max = leaves%nlen
        allocate(calculating_cells(n_calculating_cells_max))
        !
        write(str_disp, '("!", A, 2X, I5)') 'New number of cells (leaf):', leaves%nlen
        call display_string_both(str_disp, a_book_keeping%fU)
        write(str_disp, '("!", A, 2X, I5)') 'New number of cells (total):', root%nOffspring
        call display_string_both(str_disp, a_book_keeping%fU)
      else
        write(str_disp, '("! ", A)') "No further refinement needed."
        call display_string_both(str_disp, a_book_keeping%fU)
      end if
    end if
  end do
  !
  if (a_disk_iter_params%flag_converged) then
    write(*, '(A/)') "Iteration has converged!"
  else
    write(*, '(A/)') "Iteration hasn't converged. :("
  end if
  !
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=iteration_configure)
  end if
  write(str_disp, '("!Final number of cells =", I4)') leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  call line_excitation_do
  !
  ! call disk_iteration_postproc
  !
  call make_cubes
  !
end subroutine disk_iteration




subroutine post_montecarlo
  integer i, j
  integer i1, i2
  integer, parameter :: cr_TH = 10
  double precision vx, vy, vz
  double precision RR
  !
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      c%par%Tdusts = 0D0
      if (c%optical%cr_count .ge. cr_TH) then
        do j=1, dusts%n
          c%par%Tdusts(j) = get_Tdust_from_LUT( &
              c%par%en_gains(j) / (4*phy_Pi*c%par%mdusts_cell(j)), &
              luts%list(j), i1)
          !
        end do
      else
        call calc_T_diff_approx(c)
      end if
      c%par%Tdust = dot_product(c%par%Tdusts, c%par%mdusts_cell) / &
                    sum(c%par%mdusts_cell)
      !
      c%par%en_gain_tot = sum(c%par%en_gains)
      c%par%en_gain_abso_tot = sum(c%par%en_gains_abso)
      !
      ! Flux of each cell as a function of wavelength
      c%optical%flux = c%optical%flux * (phy_AU2cm / c%par%volume)
      call fill_blank(dust_0%lam, c%optical%flux, c%optical%phc, &
                      c%optical%nlam, 1, 3+c%optical%nlam/100)
      !
      ! Get some properties of the radiation field
      ! Only use the wavelength vector of dust_0
      !
      ! UV
      i1 = max(1, get_idx_for_kappa(lam_range_UV(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_UV(2), dust_0))
      c%par%flux_UV = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_UV)
      c%par%dir_UV_r = vx
      c%par%dir_UV_z = vz
      c%par%aniso_UV = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! Lya
      i1 = max(1, get_idx_for_kappa(lam_range_LyA(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_LyA(2), dust_0))
      c%par%flux_Lya = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_Lya)
      c%par%dir_Lya_r = vx
      c%par%dir_Lya_z = vz
      c%par%aniso_Lya = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! NIR
      i1 = max(1, get_idx_for_kappa(lam_range_NIR(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_NIR(2), dust_0))
      c%par%flux_NIR = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_NIR)
      c%par%dir_NIR_r = vx
      c%par%dir_NIR_z = vz
      c%par%aniso_NIR = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! MIR
      i1 = max(1, get_idx_for_kappa(lam_range_MIR(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_MIR(2), dust_0))
      c%par%flux_MIR = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_MIR)
      c%par%dir_MIR_r = vx
      c%par%dir_MIR_z = vz
      c%par%aniso_MIR = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! FIR
      i1 = max(1, get_idx_for_kappa(lam_range_FIR(1), dust_0))
      i2 = min(dust_0%n, get_idx_for_kappa(lam_range_FIR(2), dust_0))
      c%par%flux_FIR = sum(c%optical%flux(i1:i2))
      vx = sum(c%optical%dir_wei(i1:i2)%u) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      vy = sum(c%optical%dir_wei(i1:i2)%v) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      vz = sum(c%optical%dir_wei(i1:i2)%w) / c%par%volume * phy_AU2cm / (1D-100 + c%par%flux_FIR)
      c%par%dir_FIR_r = vx
      c%par%dir_FIR_z = vz
      c%par%aniso_FIR = sqrt(vx**2 + vy**2 + vz**2)
      !
      ! Local number flux of Lyman alpha
      c%par%phflux_Lya = c%par%flux_Lya / phy_LyAlpha_energy_CGS
      c%par%G0_Lya_atten = c%par%flux_Lya / phy_Habing_energy_flux_CGS
      !
      ! Calculate the total column density to the star and to the ISM
      call calc_Ncol_to_ISM(leaves%list(i)%p)
      call calc_Ncol_to_Star(leaves%list(i)%p)
      !
      RR = (c%par%rcen**2 + c%par%zcen**2) * phy_AU2cm**2
      c%par%flux_UV_star_unatten = star_0%lumi_UV0 / (4D0*phy_Pi*RR)
      !
      ! Calculate the G0 factors
      ! The G0 is the unattenuated one, so a further
      ! exp(-k*Av) should be applied.
      c%par%G0_UV_toStar = c%par%flux_UV_star_unatten / phy_Habing_energy_flux_CGS
      c%par%G0_UV_toISM  = c%par%UV_G0_factor_background
      !
      c%par%Av_toStar = max(0D0, &
        -log(c%par%flux_UV / c%par%flux_UV_star_unatten) / phy_UVext2Av)
      ! The Av to ISM is a simple scaling of the dust column density
      c%par%Av_toISM = 1.086D0 * (phy_Pi * c%par%GrainRadius_CGS**2 * 2D0) * &
                       calc_Ncol_from_cell_to_point(c, c%par%rcen, root%ymax*2D0, -5)
    end associate
  end do
end subroutine post_montecarlo



subroutine calc_T_diff_approx(c)
  type(type_cell) :: c
end subroutine calc_T_diff_approx


subroutine fill_blank(x, v, mask, n, nth, nrange)
  integer, intent(in) :: n, nth, nrange
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(inout) :: v
  integer, dimension(n), intent(in) :: mask
  integer i, j, jmin, jmax
  double precision s, smean
  do i=1, n
    if (mask(i) .lt. nth) then
      jmin = n
      jmax = 1
      do j=i-1, 1, -1
        if (mask(j) .ge. nth) then
          jmin = j
          exit
        end if
      end do
      do j=i+1, n
        if (mask(j) .ge. nth) then
          jmax = j
          exit
        end if
      end do
      jmin = min(jmin, max(1, i-nrange))
      jmax = max(jmax, min(n, i+nrange))
      s = 0D0
      do j=jmin, jmax-1
        s = s + v(j)
      end do
      smean = s / abs(x(jmax) - x(jmin))
      do j=jmin, jmax-1
        v(j) = smean * abs(x(j+1) - x(j))
      end do
    end if
  end do
end subroutine fill_blank



subroutine montecarlo_reset_cells
  integer i
  do i=1, leaves%nlen
    associate(c => leaves%list(i)%p)
      !
      c%par%Tdusts = 0D0
      !
      c%par%X_HI  = c%abundances(chem_idx_some_spe%i_HI)
      c%par%X_H2O = c%abundances(chem_idx_some_spe%i_H2O)
      !
      call calc_Ncol_to_ISM(leaves%list(i)%p)
      call calc_Ncol_to_Star(leaves%list(i)%p)
      !
      call allocate_local_optics(leaves%list(i)%p, &
                                 opmaterials%ntype, dust_0%n)
      call reset_local_optics(leaves%list(i)%p)
    end associate
  end do
  !
end subroutine montecarlo_reset_cells



subroutine disk_iteration_prepare
  integer i
  !
  ! The density structure is needed for making the grid.
  a_andrews_4ini = a_disk%andrews_gas
  a_andrews_4ini%particlemass = 1.4D0 * phy_mProton_CGS
  !
  call make_grid
  !
  n_calculating_cells_max = leaves%nlen
  allocate(calculating_cells(n_calculating_cells_max))
  !
  call prep_dust_data ! Load dust data, and create the mixtures
  !
  !call load_dust_data( &
  !  combine_dir_filename(mc_conf%mc_dir_in, mc_conf%fname_dust), dust_0)
  call make_dusts_data ! Prepare the dust optical data for use
  !
  ! Prepare the chemical stuff
  chemsol_params%fU_log = a_book_keeping%fU
  !
  call chem_read_reactions()
  call chem_load_reactions()
  call chem_parse_reactions()
  call chem_get_dupli_reactions()
  call chem_get_idx_for_special_species()
  call load_species_enthalpies
  call get_reaction_heat
  !
  call load_refine_check_species
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  call chem_load_initial_abundances
  !
  write(str_disp, '("!", A, 2X, I5)') 'Number of cells (leaf):', leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I5)') 'Number of cells (total):', root%nOffspring
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I5)') 'Number of reactions:', chem_net%nReactions
  call display_string_both(str_disp, a_book_keeping%fU)
  write(str_disp, '("!", A, 2X, I5)') 'Number of species ', chem_species%nSpecies
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  ! Set the disk and cell parameters
  call disk_set_disk_params
  call disk_set_gridcell_params
  call make_columns
  !
  if (.NOT. allocated(a_iter_stor%T_s)) then
    allocate(a_iter_stor%T_s(leaves%nlen), &
             a_iter_stor%abundances(chem_idx_some_spe%nItem, &
                                    leaves%nlen))
  end if
  !
  do i=1, leaves%nlen
    leaves%list(i)%p%abundances = chemsol_stor%y(1:chem_species%nSpecies)
    a_iter_stor%T_s(i) = leaves%list(i)%p%par%Tgas
    a_iter_stor%abundances(:, i) = chemsol_stor%y(chem_idx_some_spe%idx)
  end do
  !
  call disk_calc_disk_mass
  !
  call heating_cooling_prepare
  !
  if (a_disk_ana_params%do_analyse) then
    call load_ana_species_list
    call load_ana_points_list
    call get_species_produ_destr
    a_disk_ana_params%analyse_out_dir = &
      trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'ana/'))
    if (.not. dir_exist(a_disk_ana_params%analyse_out_dir)) then
      call my_mkdir(a_disk_ana_params%analyse_out_dir)
    end if
  end if
  !
  mc_conf%mc_dir_out = trim( &
    combine_dir_filename( &
      a_disk_iter_params%iter_files_dir, &
      mc_conf%mc_dir_out))
  if (.not. dir_exist(mc_conf%mc_dir_out)) then
    call my_mkdir(mc_conf%mc_dir_out)
  end if
  !
  star_0%mass   = a_disk%star_mass_in_Msun
  star_0%radius = a_disk%star_radius_in_Rsun
  star_0%T      = a_disk%star_temperature
  !
end subroutine disk_iteration_prepare



subroutine calc_this_cell(id)
  integer, intent(in) :: id
  integer j
  double precision tmp
  double precision R3, Z
  !
  leaves%list(id)%p%iIter = a_disk_iter_params%n_iter_used
  !
  do j=1, a_disk_iter_params%nlocal_iter
    !
    write(*, '("Local iter: ", I4, " of ", I4)') j, a_disk_iter_params%nlocal_iter
    !
    ! Set the initial condition for chemical evolution
    call set_initial_condition_4solver(id, j)
    !
    ! leaves%list(id)%p%par%Tgas = 96D0
    ! chemsol_stor%y(chem_species%nSpecies+1) = 96D0
    ! leaves%list(id)%p%par%n_gas = max(1.2D10, leaves%list(id)%p%par%n_gas)
    ! write(*, '("Tgas,ngas=", 2ES12.2, //)') leaves%list(id)%p%par%Tgas, &
    !     leaves%list(id)%p%par%n_gas
    !
    write(*, '(4X, A, F12.3/)') 'Tgas_old: ', leaves%list(id)%p%par%Tgas
    !
    call update_params_above_alt(id)
    !
    call set_chemistry_params_from_cell(id)
    call chem_cal_rates
    !
    call set_heatingcooling_params_from_cell(id)
    !
    if (chemsol_params%flag_chem_evol_save) then
      chemsol_params%chem_evol_save_filename = &
        trim(combine_dir_filename( &
          a_disk_iter_params%iter_files_dir, 'chem_evol_tmp.dat'))
    end if
    !
    call chem_set_solver_flags_alt(j)
    if (j .eq. 1) then
      chemsol_params%evolT = .true.
      chemsol_params%maySwitchT = .false.
    else if (chem_params%n_gas .gt. 1D11) then
      chemsol_params%evolT = .false.
      chemsol_params%maySwitchT = .false.
    else
      chemsol_stor%y(chem_species%nSpecies+1) = &
        (0.5D0 + dble(j)*0.1D0) * &
        (leaves%list(id)%p%par%Tgas + leaves%list(id)%p%par%Tdust)
      chemsol_params%evolT = .true.
      chemsol_params%maySwitchT = .false.
    end if
    !!!  chemsol_params%evolT = .true.
    !!!  chemsol_params%maySwitchT = .false.
    !
    call chem_evol_solve
    !
    !!!! For testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !chemsol_params%n_record_real = 1
    !if ((chem_params%rcen .le. 2D0) .and. (chem_params%zcen .le. 1D0)) then
    !  chemsol_stor%y(chem_idx_some_spe%i_H2O) = 1D-4
    !end if
    !
    if ((j .gt. 1) .and. (chemsol_params%ISTATE .eq. -3) .and. &
        (chemsol_stor%touts(chemsol_params%n_record_real) .le. &
         0.3D0 * chemsol_params%t_max)) then
      ! An unsuccessful run; will not update params
      exit
    end if
    !
    leaves%list(id)%p%abundances = chemsol_stor%y(1:chem_species%nSpecies)
    leaves%list(id)%p%par%Tgas = chemsol_stor%y(chem_species%nSpecies+1)
    leaves%list(id)%p%quality = chemsol_params%quality
    leaves%list(id)%p%par%t_final = chemsol_stor%touts(chemsol_params%n_record_real)
    !
    write(*, '(4X, A, F12.3)') 'Tgas_new: ', leaves%list(id)%p%par%Tgas
    !
    call update_params_above_alt(id)
    !
    leaves%list(id)%p%par%pressure_thermal = &
      leaves%list(id)%p%par%n_gas * &
      (leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI) + &
       leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)) * &
      leaves%list(id)%p%par%Tgas * phy_kBoltzmann_CGS
    !
    ! R3 = R^3 = (sqrt(r^2 + z^2))^3
    R3 = (sqrt((leaves%list(id)%p%par%rcen)**2 + &
               (leaves%list(id)%p%par%zcen)**2))**3
    Z = leaves%list(id)%p%par%zcen
    leaves%list(id)%p%par%gravity_z = &
      phy_GravitationConst_CGS * (star_0%mass * phy_Msun_CGS) * &
      (leaves%list(id)%p%par%mgas_cell + &
       leaves%list(id)%p%par%mdust_tot) * &
      (-Z / R3 / (phy_AU2cm**2))
    !
    if ((leaves%list(id)%p%quality .eq. 0) .or. &
        ((j .ge. 2) .and. &
         (leaves%list(id)%p%par%t_final .ge. &
          0.3D0 * chemsol_params%t_max))) then
      exit
    end if
    !
    ! Update local velocity width
    call calc_local_dynamics(leaves%list(id)%p)
    !
  end do
  !
  if (.not. chemsol_params%evolT) then
    call chem_cal_rates
    call realtime_heating_cooling_rate(tmp, chemsol_params%NEQ, chemsol_stor%y)
  end if
  leaves%list(id)%p%h_c_rates = heating_cooling_rates
  !
  if (a_disk_iter_params%flag_save_rates) then
    call save_chem_rates(id)
  end if
  !
  if (a_disk_ana_params%do_analyse) then
    if ((leaves%list(id)%p%quality .gt. 0) .or. &
        ((a_disk_iter_params%n_iter_used .gt. 1) .and. &
         (is_in_list_int(id, ana_ptlist%nlen, ana_ptlist%vals) .or. &
          need_to_refine(leaves%list(id)%p)))) then
      call chem_analyse(id)
    end if
  end if
end subroutine calc_this_cell



subroutine update_params_above_alt(i0)
  use load_Visser_CO_selfshielding
  integer, intent(in) :: i0
  integer i
  !
  ! Calculate the column density of a few species to the star and to the ISM
  do i=1, chem_idx_some_spe%nItem
    call calc_Ncol_to_ISM(leaves%list(i0)%p, i)
    call calc_Ncol_to_Star(leaves%list(i0)%p, i)
  end do
  associate(c => leaves%list(i0)%p)
    ! Kwok eq 10.20
    ! c%par%Av_toISM = 1.086D0 * c%par%ratioDust2HnucNum * &
    !   (phy_Pi * c%par%GrainRadius_CGS**2) * 2D0 * c%par%Ncol_toISM
    ! c%par%Av_toStar = 1.086D0 * c%par%ratioDust2HnucNum * &
    !   (phy_Pi * c%par%GrainRadius_CGS**2) * 2D0 * c%par%Ncol_toStar
    !
    c%par%f_selfshielding_toISM_H2  = min(1D0, get_H2_self_shielding( &
      c%col_den_toISM(chem_idx_some_spe%iiH2), c%par%velo_width_turb))
    c%par%f_selfshielding_toStar_H2  = min(1D0, get_H2_self_shielding( &
      c%col_den_toStar(chem_idx_some_spe%iiH2), c%par%velo_width_turb))
    !
    ! H2O and OH self shielding are already taken into account in the radiative transfer.
    ! Only for output; not used.
    c%par%f_selfshielding_toISM_H2O = &
      min(1D0, exp(-(c%col_den_toISM(chem_idx_some_spe%iiH2O) * const_LyAlpha_cross_H2O)))
    c%par%f_selfshielding_toStar_H2O = &
      min(1D0, exp(-(c%col_den_toStar(chem_idx_some_spe%iiH2O) * const_LyAlpha_cross_H2O)))
    !
    c%par%f_selfshielding_toISM_OH = &
      min(1D0, exp(-(c%col_den_toISM(chem_idx_some_spe%iiOH) * const_LyAlpha_cross_OH)))
    c%par%f_selfshielding_toStar_OH = &
      min(1D0, exp(-(c%col_den_toStar(chem_idx_some_spe%iiOH) * const_LyAlpha_cross_OH)))
    !
    c%par%f_selfshielding_toISM_CO = min(1D0, max(0D0, get_12CO_shielding( &
      c%col_den_toISM(chem_idx_some_spe%iiH2), &
      c%col_den_toISM(chem_idx_some_spe%iiCO))))
    !
    c%par%f_selfshielding_toStar_CO = min(1D0, max(0D0, get_12CO_shielding( &
      c%col_den_toStar(chem_idx_some_spe%iiH2), &
      c%col_den_toStar(chem_idx_some_spe%iiCO))))
    !
    ! Calculate the gravitational force from above
    if (c%above%n .eq. 0) then
      c%par%gravity_acc_z = &
        phy_GravitationConst_CGS * (star_0%mass * phy_Msun_CGS) * &
          (calc_Ncol_from_cell_to_point(leaves%list(i0)%p, &
                                        c%par%rcen, root%ymax*2D0, -4) * &
           c%par%area_T * phy_mProton_CGS * c%par%MeanMolWeight) * &
          (-c%ymax / (sqrt(c%xmax**2 + c%ymax**2))**3 / (phy_AU2cm**2))
    else
      c%par%gravity_acc_z = &
        calc_Ncol_from_cell_to_point(leaves%list(i0)%p, &
                                     c%par%rcen, root%ymax*2D0, -1)
    end if
  end associate
end subroutine update_params_above_alt



function get_H2_self_shielding(N_H2, dv_turb)
  ! Draine 1996, equation 37
  double precision get_H2_self_shielding
  double precision, intent(in) :: N_H2, dv_turb
  double precision x, b5
  x = N_H2 / 5D14
  b5 = dv_turb / 1D5
  get_H2_self_shielding = 0.965D0 / (1D0 + x/b5)**2 + &
    0.035 / sqrt(1D0 + x) * exp(-8.5D-4 * sqrt(1D0 + x))
end function get_H2_self_shielding



subroutine check_convergency_cell(i0)
  integer, intent(in) :: i0
  ! Temperature is not considered.
  if (maxval(abs(leaves%list(i0)%p%abundances(chem_idx_some_spe%idx) &
                 - a_iter_stor%abundances(:, i0)) &
             - (a_disk_iter_params%atol_abun + &
                a_disk_iter_params%rtol_abun * &
                abs(leaves%list(i0)%p%abundances(chem_idx_some_spe%idx) &
                  + a_iter_stor%abundances(:, i0))) &
            ) .le. 0D0) then
    leaves%list(i0)%p%converged = .true.
  else
    leaves%list(i0)%p%converged = .false.
  end if
end subroutine check_convergency_cell



subroutine check_convergency_whole_disk
  integer i
  a_disk_iter_params%n_cell_converged = 0
  do i=1, leaves%nlen
    if (leaves%list(i)%p%converged) then
      a_disk_iter_params%n_cell_converged = a_disk_iter_params%n_cell_converged + 1
    end if
  end do
  a_disk_iter_params%flag_converged = &
    a_disk_iter_params%n_cell_converged .ge. &
    int(a_disk_iter_params%converged_cell_percentage_stop * real(leaves%nlen))
  write(str_disp, '("! Iter", I4, 4X, "Number of cells converged: ", I6, "/", I6)') &
    a_disk_iter_params%n_iter_used, a_disk_iter_params%n_cell_converged, leaves%nlen
  call display_string_both(str_disp, a_book_keeping%fU)
end subroutine check_convergency_whole_disk



subroutine update_calculating_cells
  integer, dimension(:), allocatable :: list_tmp
  integer i, i0, itmp, j, k, n
  logical flag_notyet
  allocate(list_tmp(n_calculating_cells_max))
  n = 0
  do i=1, n_calculating_cells
    i0 = calculating_cells(i)
    do j=1, leaves%list(i0)%p%below%n
      itmp = leaves%list(i0)%p%below%idx(j)
      flag_notyet = .true.
      do k=1, n
        if (list_tmp(k) .eq. itmp) then
          flag_notyet = .false.
          exit
        end if
      end do
      if (flag_notyet) then
        n = n + 1
        list_tmp(n) = itmp
      end if
    end do
  end do
  n_calculating_cells = n
  if (n .ge. 1) then
    calculating_cells(1:n) = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine update_calculating_cells


subroutine set_initial_condition_4solver(id, iloc_iter)
  integer, intent(in) :: id, iloc_iter
  !double precision tmp
  !logical found_neighbor
  !integer i, i0, ntmp
  !
  !if (a_disk_iter_params%flag_shortcut_ini) then
  !  if (a_disk_iter_params%n_iter_used .eq. 1) then
  !    found_neighbor = .false.
  !    do i=1, leaves%list(id)%p%around%n
  !      i0 = leaves%list(id)%p%around%idx(i)
  !      if (leaves%list(i0)%p%iIter .gt. leaves%list(id)%p%iIter) then
  !        chemsol_stor%y(1:chem_species%nSpecies) = leaves%list(i0)%p%abundances
  !        found_neighbor = .true.
  !        exit
  !      end if
  !    end do
  !    if (.not. found_neighbor) then
  !      chemsol_stor%y(1:chem_species%nSpecies) = &
  !        chemsol_stor%y0(1:chem_species%nSpecies)
  !    end if
  !  else
  !    chemsol_stor%y(1:chem_species%nSpecies) = leaves%list(id)%p%abundances
  !  end if
  !  if (chemsol_params%neutralize) then
  !    tmp = sum(chemsol_stor%y(1:chem_species%nSpecies) * &
  !                      dble(chem_species%elements(1,:)))
  !    if (abs(tmp) .ge. 1D-2*chemsol_stor%y(chem_idx_some_spe%i_E)) then
  !      chemsol_stor%y(1:chem_species%nSpecies) = &
  !        chemsol_stor%y0(1:chem_species%nSpecies)
  !    else
  !      chemsol_stor%y(chem_idx_some_spe%i_E) = &
  !        chemsol_stor%y(chem_idx_some_spe%i_E) + tmp
  !      if (chemsol_stor%y(chem_idx_some_spe%i_E) .lt. 0D0) then
  !        ! When it is not possible to neutralize the composition by artificially
  !        ! changing the electron abundance, then use the general initial abundances,
  !        ! which should be absolutely neutral.
  !        chemsol_stor%y(1:chem_species%nSpecies) = &
  !          chemsol_stor%y0(1:chem_species%nSpecies)
  !        write(str_disp, '("! Cannot neutralize: X(E-) = ", ES12.4)') &
  !          chemsol_stor%y(chem_idx_some_spe%i_E)
  !        call display_string_both(str_disp, a_book_keeping%fU)
  !        write(str_disp, '("! Use y0 as initial abundance.")')
  !        call display_string_both(str_disp, a_book_keeping%fU)
  !        write(str_disp, '("! x, y = ", 2ES10.2, " iIter = ", I4)') &
  !          leaves%list(id)%p%xmin, leaves%list(id)%p%ymin, leaves%list(id)%p%iIter
  !        call display_string_both(str_disp, a_book_keeping%fU)
  !      end if
  !    end if
  !  end if
  !else
  !  chemsol_stor%y(1:chem_species%nSpecies) = &
  !    chemsol_stor%y0(1:chem_species%nSpecies)
  !end if
  chemsol_stor%y(1:chem_species%nSpecies) = &
    chemsol_stor%y0(1:chem_species%nSpecies)
  !
  ! Initial abundance of *neutral* dust
  ! There should be no dust in the input initial abundance file
  if (chem_idx_some_spe%i_Grain0 .ne. 0) then
    chemsol_stor%y(chem_idx_some_spe%i_Grain0) = &
        leaves%list(id)%p%par%ratioDust2HnucNum
  end if
  !
  ! Always use the temperature of the above cell as init
  !if ((a_disk_iter_params%n_iter_used .eq. 1) .and. (iloc_iter .eq. 1)) then
  !  if (leaves%list(id)%p%above%n .gt. 0) then
  !    leaves%list(id)%p%par%Tgas = 0D0
  !    ntmp = 0
  !    do i=1, leaves%list(id)%p%above%n
  !      i0 = leaves%list(id)%p%above%idx(i)
  !      if (leaves%list(i0)%p%iIter .lt. leaves%list(id)%p%iIter) then
  !        cycle
  !      end if 
  !      ntmp = ntmp + 1
  !      leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tgas + &
  !                                    leaves%list(i0)%p%par%Tgas
  !    end do
  !    if (ntmp .gt. 0) then
  !      leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tgas / &
  !                                    dble(ntmp)
  !    else
  !      leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tdust
  !    end if
  !  else
  !    leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tdust
  !  end if
  !end if
  leaves%list(id)%p%par%Tgas = leaves%list(id)%p%par%Tdust
  !
  ! Set the initial temeprature
  chemsol_stor%y(chem_species%nSpecies+1) = leaves%list(id)%p%par%Tgas
  !
end subroutine set_initial_condition_4solver



subroutine vertical_pressure_gravity_balance
  integer i
  double precision dznew, pnew, nnew, vnew
  if (.not. grid_config%columnwise) then
    return
  end if
  !
  write(str_disp, '(A)') '! Adjusting the vertical structure.'
  call display_string_both(str_disp, a_book_keeping%fU)
  !
  ! Calculate the new density and size of each cell
  do i=1, leaves%nlen
    pnew = -leaves%list(i)%p%par%gravity_acc_z / leaves%list(i)%p%par%area_T
    ! Not change too much in one go
    pnew = sqrt(pnew * leaves%list(i)%p%par%pressure_thermal)
    !
    nnew = pnew / &
      (leaves%list(i)%p%par%Tgas * phy_kBoltzmann_CGS * &
       (leaves%list(i)%p%abundances(chem_idx_some_spe%i_HI) + &
        leaves%list(i)%p%abundances(chem_idx_some_spe%i_H2)))
    vnew = leaves%list(i)%p%par%mgas_cell / nnew / &
      (phy_mProton_CGS * leaves%list(i)%p%par%MeanMolWeight)
    dznew = vnew / leaves%list(i)%p%par%area_T / phy_AU2cm
    !
    leaves%list(i)%p%par%n_dusts = leaves%list(i)%p%par%n_dusts * nnew / &
                                   leaves%list(i)%p%par%n_gas
    leaves%list(i)%p%par%n_gas = nnew
    leaves%list(i)%p%val(1) = nnew
    leaves%list(i)%p%par%volume = vnew
    leaves%list(i)%p%par%dz = dznew
    leaves%list(i)%p%par%pressure_thermal = pnew
  end do
  !
  ! Move the cells column by column from bottom to top.
  call shift_and_scale_above
  !
  ! Remake neigbor information
  call grid_make_neighbors
  !
  call grid_make_surf_bott
  !
  mc_conf%maxw = get_surf_max_angle()
  call get_mc_stellar_par(mc_conf)
  !
  write(str_disp, '(A, 2F9.4)') 'New minw,maxw: ', mc_conf%minw, mc_conf%maxw
  call display_string_both(str_disp, a_book_keeping%fU)
  write(*, *)
  !
end subroutine vertical_pressure_gravity_balance



subroutine shift_and_scale_above
  type(type_cell), pointer :: cthis
  double precision ybelow
  integer i, ic
  !
  do ic=1, bott_cells%nlen ! Loop over the columns
    ybelow = columns(ic)%list(1)%p%ymin
    do i=1, columns(ic)%nlen
      cthis => columns(ic)%list(i)%p
      if (cthis%using) then
        cthis%ymin = ybelow
        cthis%ymax = ybelow + cthis%par%dz
        !
        cthis%par%zmin = cthis%ymin
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        ! dz is already set
        !
        cthis%par%area_I = phy_2Pi * cthis%par%rmin * cthis%par%dz * phy_AU2cm**2
        cthis%par%area_O = phy_2Pi * cthis%par%rmax * cthis%par%dz * phy_AU2cm**2
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
        !
      else
        cthis%ymax = ybelow + (cthis%ymax - cthis%ymin)
        cthis%ymin = ybelow
      end if
      !
      root%ymax = max(root%ymax, cthis%ymax)
      !
      ybelow = cthis%ymax
      !
    end do
  end do
  !
  ! Align the top edge to the domain upper boundary
  do ic=1, bott_cells%nlen
    i = columns(ic)%nlen
    cthis => columns(ic)%list(i)%p
    if (cthis%ymax .le. root%ymax) then
      ! Rescale the density
      cthis%val(1) = cthis%val(1) * &
        (cthis%ymax - cthis%ymin) / (root%ymax - cthis%ymin)
      cthis%ymax = root%ymax
      !
      if (cthis%using) then
        cthis%par%zmax = cthis%ymax
        cthis%par%zcen = (cthis%ymax + cthis%ymin) * 0.5D0
        cthis%par%dz = cthis%par%zmax - cthis%par%zmin
        !
        cthis%par%n_dusts = cthis%par%n_dusts * cthis%val(1) / cthis%par%n_gas
        cthis%par%n_gas  = cthis%val(1)
        !
        cthis%par%area_I = phy_2Pi * cthis%par%rmin * cthis%par%dz * phy_AU2cm**2
        cthis%par%area_O = phy_2Pi * cthis%par%rmax * cthis%par%dz * phy_AU2cm**2
        cthis%par%surf_area = cthis%par%area_T + cthis%par%area_B + &
          cthis%par%area_I + cthis%par%area_O
      end if
    else
      write(*,'(A)') 'Should not have this case:'
      write(*,'(A/)') 'in shift_and_scale_above.'
    end if
  end do
end subroutine shift_and_scale_above


subroutine make_columns
  integer i, j
  type(type_cell), pointer :: cthis, cnext
  double precision length, r, z, eps
  integer dirtype
  logical found
  type(type_ray) ray
  !
  allocate(columns(bott_cells%nlen))
  do i=1, bott_cells%nlen
    cthis => leaves%list(bott_cells%idx(i))%p
    r = cthis%par%rcen
    z = root%ymax * 2D0
    columns(i)%nlen = int(calc_Ncol_from_cell_to_point(cthis, r, z, -2))
    if (.not. allocated(columns(i)%list)) then
      allocate(columns(i)%list(columns(i)%nlen))
    end if
    !
    ray%x = cthis%par%rcen
    ray%y = 0D0
    ray%z = cthis%par%zcen
    ray%vx = 0D0
    ray%vy = 0D0
    ray%vz = 1D0
    !
    j = 0
    ! First make a list of all the cells (including null ones) above the bottom
    ! (mid-plane) one.  They are ordered from bottom to top.
    do
      call calc_intersection_ray_cell(ray, cthis, &
        length, r, z, eps, found, dirtype)
      if (.not. found) then
        write(str_disp,'(A)') 'In make_columns:'
        call display_string_both(str_disp, a_book_keeping%fU)
        write(str_disp,'(A, 9ES16.6)') 'ray does not intersect cthis: ', &
          sqrt(ray%x**2+ray%y**2), ray%z, ray%vx, ray%vy, ray%vz, &
            cthis%xmin, cthis%xmax, cthis%ymin, cthis%ymax
        call display_string_both(str_disp, a_book_keeping%fU)
        write(*, *)
        return
      end if
      !
      j = j + 1
      columns(i)%list(j)%p => cthis
      !
      ray%x = ray%x + ray%vx * (length + eps)
      ray%y = ray%y + ray%vy * (length + eps)
      ray%z = ray%z + ray%vz * (length + eps)
      !
      call locate_photon_cell_alt(r, z, cthis, dirtype, cnext, found)
      if (found) then
        cthis => cnext
      else
        exit
      end if
    end do
  end do
end subroutine make_columns



subroutine calc_Ncol_to_ISM(c, iSp)
  ! iSp is the index in chem_idx_some_spe, not in the range 1 to
  ! chem_species%nSpecies
  type(type_cell), intent(inout), pointer :: c
  integer, intent(in), optional :: iSp
  if (present(iSp)) then
    c%col_den_toISM(iSp) = calc_Ncol_from_cell_to_point( &
      c, c%par%rcen, root%ymax * 2D0, &
      chem_idx_some_spe%idx(iSp))
  else
    c%par%Ncol_toISM = calc_Ncol_from_cell_to_point( &
      c, c%par%rcen, root%ymax * 2D0)
  end if
end subroutine calc_Ncol_to_ISM



subroutine calc_Ncol_to_Star(c, iSp)
  ! iSp is the index in chem_idx_some_spe, not in the range 1 to
  ! chem_species%nSpecies
  type(type_cell), intent(inout), pointer :: c
  integer, intent(in), optional :: iSp
  if (present(iSp)) then
    c%col_den_toStar(iSp) = calc_Ncol_from_cell_to_point( &
      c, 0D0, 0D0, chem_idx_some_spe%idx(iSp))
  else
    c%par%Ncol_toStar = calc_Ncol_from_cell_to_point( &
      c, 0D0, 0D0)
  end if
end subroutine calc_Ncol_to_Star



function calc_Ncol_from_cell_to_point(c, r, z, iSpe) result(N)
  double precision N
  type(type_cell), intent(in), target :: c ! Todo: pointer or not?
  double precision, intent(in) :: r, z
  integer, intent(in), optional :: iSpe
  type(type_ray) ray
  type(type_cell), pointer :: cthis, cnext
  double precision t, length, r1, z1, eps
  logical found
  integer dirtype
  double precision, parameter :: small_dist = 1D-50
  !
  N = 0D0
  !
  ray%x = c%par%rcen
  ray%y = 0D0
  ray%z = c%par%zcen
  !
  ray%vx = r - ray%x
  ray%vy = 0D0
  ray%vz = z - ray%z
  t = sqrt(ray%vx**2 + ray%vy**2 + ray%vz**2)
  if (t .le. small_dist) then
    return
  end if
  ray%vx = ray%vx / t
  ray%vy = ray%vy / t
  ray%vz = ray%vz / t
  !
  cthis => c
  !
  do
    call calc_intersection_ray_cell(ray, cthis, &
        length, r1, z1, eps, found, dirtype)
    if (.not. found) then
      write(str_disp,'(A, 9ES13.4)') &
        'In calc_Ncol_from_cell_to_point: ray does not intersect cthis: ', &
        sqrt(ray%x**2+ray%y**2), ray%z, ray%vx, ray%vy, ray%vz, &
            cthis%xmin, cthis%xmax, cthis%ymin, cthis%ymax
      call display_string_both(str_disp, a_book_keeping%fU)
      exit
    end if
    !
    if (present(iSpe)) then
      if ((iSpe .ge. 1) .and. (iSpe .le. chem_species%nSpecies)) then
        ! Calculate column density of a species
        if (cthis%using) then
          N = N + cthis%par%n_gas * cthis%abundances(iSpe) * &
              (length * phy_AU2cm)
        end if
      !
      else if (iSpe .eq. 0) then
        ! Calculate column density
        if (cthis%using) then
          N = N + cthis%par%n_gas * length * phy_AU2cm
        end if
      !
      else if (iSpe .eq. -1) then
        ! Calculate gravity force (accumulated), not including the gravity of
        ! the starting cell.
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          if (cthis%using) then
            N = N + cthis%par%gravity_z
          end if
        end if
      !
      else if (iSpe .eq. -2) then
        ! Calculate the number of cells crossed by the ray.
        N = N + 1D0
      !
      else if (iSpe .eq. -3) then
        ! Calculate column density, including null cells.
        ! Density is taken from the initial value
        N = N + cthis%val(1) * length * phy_AU2cm
      !
      else if (iSpe .eq. -4) then
        ! Calculate column density, including null cells, but not itself.
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          N = N + cthis%val(1) * length * phy_AU2cm
        end if
      !
      else if (iSpe .eq. -5) then
        ! Calculate the dust column density
        if ((cthis%xmin .ne. c%xmin) .or. (cthis%xmax .ne. c%xmax) .or. &
            (cthis%ymin .ne. c%ymin) .or. (cthis%ymax .ne. c%ymax)) then
          if (cthis%using) then
            N = N + cthis%par%ndust_tot * length * phy_AU2cm
          end if
        end if
      else
        write(str_disp,'(A)') 'I do not know what to do!'
        call display_string_both(str_disp, a_book_keeping%fU)
        stop
      end if
    else
      ! Calculate column density by default
      if (cthis%using) then
        N = N + cthis%par%n_gas * length * phy_AU2cm
      end if
    end if
    !
    ray%x = ray%x + ray%vx * (length + eps)
    ray%y = ray%y + ray%vy * (length + eps)
    ray%z = ray%z + ray%vz * (length + eps)
    !
    call locate_photon_cell_alt(r1, z1, cthis, dirtype, cnext, found)
    if (found) then
      cthis => cnext
    else
      exit
    end if
  end do
  !
  if (present(iSpe)) then
    if (iSpe .eq. -1) then
      if (cthis%using) then
        N = N + cthis%par%gravity_acc_z
      end if
    end if
  end if
end function calc_Ncol_from_cell_to_point



subroutine disk_save_results_pre
  if (.NOT. getFileUnit(fU_save_results)) then
    write(*,*) 'Cannot get a file unit for output!'
    stop
  end if
  write(filename_save_results, '("iter_", I4.4, ".dat")') &
    a_disk_iter_params%n_iter_used
  filename_save_results = trim(combine_dir_filename( &
    a_disk_iter_params%iter_files_dir, filename_save_results))
  call openFileSequentialWrite(fU_save_results, filename_save_results, 99999)
  !
  call write_header(fU_save_results)
end subroutine disk_save_results_pre


subroutine write_header(fU)
  integer, intent(in) :: fU
  character(len=64) fmt_str
  character(len=8192) tmp_str
  write(fmt_str, '("(", I4, "A14)")') chem_species%nSpecies
  write(tmp_str, fmt_str) chem_species%names
  write(fU, '(A)') &
    '!' // &
    str_pad_to_len('cvg', 4) // &
    str_pad_to_len('qual', 5) // &
    str_pad_to_len('arnd', 5) // &
    str_pad_to_len('abov', 5) // &
    str_pad_to_len('belo', 5) // &
    str_pad_to_len('innr', 5) // &
    str_pad_to_len('outr', 5) // &
    str_pad_to_len('cr_count',len_item) // &
    str_pad_to_len('abc_dus', len_item) // &
    str_pad_to_len('scc_HI',  len_item) // &
    str_pad_to_len('abc_wat', len_item) // &
    str_pad_to_len('t_final', len_item) // &
    str_pad_to_len('rmin',    len_item) // &
    str_pad_to_len('rmax',    len_item) // &
    str_pad_to_len('zmin',    len_item) // &
    str_pad_to_len('zmax',    len_item) // &
    str_pad_to_len('Tgas',    len_item) // &
    str_pad_to_len('Tdust',   len_item) // &
    str_pad_to_len('n_gas',   len_item) // &
    str_pad_to_len('ndust_t', len_item) // &
    str_pad_to_len('d2gmas',  len_item) // &
    str_pad_to_len('d2gnum',  len_item) // &
    str_pad_to_len('deplet',  len_item) // &
    str_pad_to_len('mg_cell', len_item) // &
    str_pad_to_len('md_cell', len_item) // &
    str_pad_to_len('presr_t', len_item) // &
    str_pad_to_len('grav_z',  len_item) // &
    str_pad_to_len('egain_d', len_item) // &
    str_pad_to_len('egain_ab',len_item) // &
    str_pad_to_len('flx_UV',  len_item) // &
    str_pad_to_len('flx_Lya', len_item) // &
    str_pad_to_len('flx_NIR', len_item) // &
    str_pad_to_len('flx_MIR', len_item) // &
    str_pad_to_len('flx_FIR', len_item) // &
    str_pad_to_len('vr_UV',   len_item) // &
    str_pad_to_len('vz_UV',   len_item) // &
    str_pad_to_len('ani_UV',  len_item) // &
    str_pad_to_len('vr_Lya',  len_item) // &
    str_pad_to_len('vz_Lya',  len_item) // &
    str_pad_to_len('ani_Lya', len_item) // &
    str_pad_to_len('vr_NIR',  len_item) // &
    str_pad_to_len('vz_NIR',  len_item) // &
    str_pad_to_len('ani_NIR', len_item) // &
    str_pad_to_len('vr_MIR',  len_item) // &
    str_pad_to_len('vz_MIR',  len_item) // &
    str_pad_to_len('ani_MIR', len_item) // &
    str_pad_to_len('vr_FIR',  len_item) // &
    str_pad_to_len('vz_FIR',  len_item) // &
    str_pad_to_len('ani_FIR', len_item) // &
    str_pad_to_len('Av_ISM',  len_item) // &
    str_pad_to_len('Av_Star', len_item) // &
    str_pad_to_len('UV_G0_I', len_item) // &
    str_pad_to_len('UV_G0_S', len_item) // &
    str_pad_to_len('LyAG0_a', len_item) // &
    str_pad_to_len('LyANF0',  len_item) // &
    str_pad_to_len('XRay0',   len_item) // &
    str_pad_to_len('Ncol_I',  len_item) // &
    str_pad_to_len('Ncol_S',  len_item) // &
    str_pad_to_len('N_H2_I',  len_item) // &
    str_pad_to_len('N_H2O_I', len_item) // &
    str_pad_to_len('N_OH_I',  len_item) // &
    str_pad_to_len('N_CO_I',  len_item) // &
    str_pad_to_len('N_H2_S',  len_item) // &
    str_pad_to_len('N_H2O_S', len_item) // &
    str_pad_to_len('N_OH_S',  len_item) // &
    str_pad_to_len('N_CO_S',  len_item) // &
    str_pad_to_len('f_H2_I',  len_item) // &
    str_pad_to_len('f_H2O_I', len_item) // &
    str_pad_to_len('f_OH_I',  len_item) // &
    str_pad_to_len('f_CO_I',  len_item) // &
    str_pad_to_len('f_H2_S',  len_item) // &
    str_pad_to_len('f_H2O_S', len_item) // &
    str_pad_to_len('f_OH_S',  len_item) // &
    str_pad_to_len('f_CO_S',  len_item) // &
    str_pad_to_len('R_H2_fo', len_item) // &
    str_pad_to_len('hc_net',  len_item) // &
    str_pad_to_len('h_ph_gr', len_item) // &
    str_pad_to_len('h_fo_H2', len_item) // &
    str_pad_to_len('h_cosmi', len_item) // &
    str_pad_to_len('h_vi_H2', len_item) // &
    str_pad_to_len('h_io_CI', len_item) // &
    str_pad_to_len('h_ph_H2', len_item) // &
    str_pad_to_len('h_ph_wa', len_item) // &
    str_pad_to_len('h_ph_OH', len_item) // &
    str_pad_to_len('h_Xray ', len_item) // &
    str_pad_to_len('h_visco', len_item) // &
    str_pad_to_len('h_chem',  len_item) // &
    str_pad_to_len('c_el_gr', len_item) // &
    str_pad_to_len('c_vi_H2', len_item) // &
    str_pad_to_len('c_gg_co', len_item) // &
    str_pad_to_len('c_OI   ', len_item) // &
    str_pad_to_len('c_CII  ', len_item) // &
    str_pad_to_len('c_wa_ro', len_item) // &
    str_pad_to_len('c_wa_vi', len_item) // &
    str_pad_to_len('c_CO_ro', len_item) // &
    str_pad_to_len('c_CO_vi', len_item) // &
    str_pad_to_len('c_H2_ro', len_item) // &
    str_pad_to_len('c_LyAlp', len_item) // &
    str_pad_to_len('c_fb   ', len_item) // &
    str_pad_to_len('c_ff   ', len_item) // &
    trim(tmp_str)
end subroutine write_header


subroutine disk_save_results_write(fU, c)
  character(len=64) fmt_str
  integer, intent(in) :: fU
  type(type_cell), pointer, intent(in) :: c
  integer converged
  !
  write(fmt_str, '(", ", I4, "ES14.4E4)")') chem_species%nSpecies
  if (c%converged) then
    converged = 1
  else
    converged = 0
  end if
  write(fU, '(7I5, 4I14, 89ES14.5E3' // trim(fmt_str)) &
  converged                                              , &
  c%quality                                              , &
  c%around%n                                             , &
  c%above%n                                              , &
  c%below%n                                              , &
  c%inner%n                                              , &
  c%outer%n                                              , &
  c%optical%cr_count                                     , &
  c%par%ab_count_dust                                    , &
  c%par%sc_count_HI                                      , &
  c%par%ab_count_water                                   , &
  c%par%t_final                                          , &
  c%par%rmin                                             , &
  c%par%rmax                                             , &
  c%par%zmin                                             , &
  c%par%zmax                                             , &
  c%par%Tgas                                             , &
  c%par%Tdust                                            , &
  c%par%n_gas                                            , &
  c%par%ndust_tot                                        , &
  c%par%ratioDust2GasMass                                , &
  c%par%ratioDust2HnucNum                                , &
  c%par%dust_depletion                                   , &
  c%par%mgas_cell                                        , &
  c%par%mdust_tot                                        , &
  c%par%pressure_thermal                                 , &
  c%par%gravity_acc_z                                    , &
  c%par%en_gain_tot                                      , &
  c%par%en_gain_abso_tot                                 , &
  c%par%flux_UV                                          , &
  c%par%flux_Lya                                         , &
  c%par%flux_NIR                                         , &
  c%par%flux_MIR                                         , &
  c%par%flux_FIR                                         , &
  c%par%dir_UV_r                                         , &
  c%par%dir_UV_z                                         , &
  c%par%aniso_UV                                         , &
  c%par%dir_Lya_r                                        , &
  c%par%dir_Lya_z                                        , &
  c%par%aniso_Lya                                        , &
  c%par%dir_NIR_r                                        , &
  c%par%dir_NIR_z                                        , &
  c%par%aniso_NIR                                        , &
  c%par%dir_MIR_r                                        , &
  c%par%dir_MIR_z                                        , &
  c%par%aniso_MIR                                        , &
  c%par%dir_FIR_r                                        , &
  c%par%dir_FIR_z                                        , &
  c%par%aniso_FIR                                        , &
  c%par%Av_toISM                                         , &
  c%par%Av_toStar                                        , &
  c%par%G0_UV_toISM                                      , &
  c%par%G0_UV_toStar                                     , &
  c%par%G0_Lya_atten                                     , &
  c%par%phflux_Lya                                       , &
  c%par%Xray_flux_0                                      , &
  c%par%Ncol_toISM                                       , &
  c%par%Ncol_toStar                                      , &
  c%col_den_toISM(chem_idx_some_spe%iiH2)                , &
  c%col_den_toISM(chem_idx_some_spe%iiH2O)               , &
  c%col_den_toISM(chem_idx_some_spe%iiOH)                , &
  c%col_den_toISM(chem_idx_some_spe%iiCO)                , &
  c%col_den_toStar(chem_idx_some_spe%iiH2)               , &
  c%col_den_toStar(chem_idx_some_spe%iiH2O)              , &
  c%col_den_toStar(chem_idx_some_spe%iiOH)               , &
  c%col_den_toStar(chem_idx_some_spe%iiCO)               , &
  c%par%f_selfshielding_toISM_H2                         , &
  c%par%f_selfshielding_toISM_H2O                        , &
  c%par%f_selfshielding_toISM_OH                         , &
  c%par%f_selfshielding_toISM_CO                         , &
  c%par%f_selfshielding_toStar_H2                        , &
  c%par%f_selfshielding_toStar_H2O                       , &
  c%par%f_selfshielding_toStar_OH                        , &
  c%par%f_selfshielding_toStar_CO                        , &
  c%par%R_H2_form_rate                                   , &
  c%h_c_rates%hc_net_rate                                , &
  c%h_c_rates%heating_photoelectric_small_grain_rate     , &
  c%h_c_rates%heating_formation_H2_rate                  , &
  c%h_c_rates%heating_cosmic_ray_rate                    , &
  c%h_c_rates%heating_vibrational_H2_rate                , &
  c%h_c_rates%heating_ionization_CI_rate                 , &
  c%h_c_rates%heating_photodissociation_H2_rate          , &
  c%h_c_rates%heating_photodissociation_H2O_rate         , &
  c%h_c_rates%heating_photodissociation_OH_rate          , &
  c%h_c_rates%heating_Xray_Bethell_rate                  , &
  c%h_c_rates%heating_viscosity_rate                     , &
  c%h_c_rates%heating_chem                               , &
  c%h_c_rates%cooling_photoelectric_small_grain_rate     , &
  c%h_c_rates%cooling_vibrational_H2_rate                , &
  c%h_c_rates%cooling_gas_grain_collision_rate           , &
  c%h_c_rates%cooling_OI_rate                            , &
  c%h_c_rates%cooling_CII_rate                           , &
  c%h_c_rates%cooling_Neufeld_H2O_rate_rot               , &
  c%h_c_rates%cooling_Neufeld_H2O_rate_vib               , &
  c%h_c_rates%cooling_Neufeld_CO_rate_rot                , &
  c%h_c_rates%cooling_Neufeld_CO_rate_vib                , &
  c%h_c_rates%cooling_Neufeld_H2_rot_rate                , &
  c%h_c_rates%cooling_LymanAlpha_rate                    , &
  c%h_c_rates%cooling_free_bound_rate                    , &
  c%h_c_rates%cooling_free_free_rate                     , &
  c%abundances
end subroutine disk_save_results_write


subroutine disk_calc_disk_mass
  integer i
  a_disk%disk_mass_in_Msun = 0D0
  do i=1, leaves%nlen
    associate(p => leaves%list(i)%p%par)
      a_disk%disk_mass_in_Msun = &
        a_disk%disk_mass_in_Msun + &
        p%n_gas * p%MeanMolWeight * (phy_2Pi * p%rcen * p%dr * p%dz)
    end associate
  end do
  a_disk%disk_mass_in_Msun = a_disk%disk_mass_in_Msun * &
    phy_AU2cm**3 * phy_mProton_CGS / phy_Msun_CGS
end subroutine disk_calc_disk_mass


subroutine set_heatingcooling_params_from_cell(id)
  integer id
  hc_params%type_cell_rz_phy_basic = leaves%list(id)%p%par
  hc_params%Neufeld_dv_dz = leaves%list(id)%p%par%velo_gradient * 1D-5 ! cm s-1 to km s-1
  hc_params%Neufeld_G     = 1D0
  hc_params%X_H2    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2)
  hc_params%X_HI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_HI)
  hc_params%X_CI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CI)
  hc_params%X_Cplus = leaves%list(id)%p%abundances(chem_idx_some_spe%i_Cplus)
  hc_params%X_OI    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_OI)
  hc_params%X_CO    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_CO)
  hc_params%X_H2O   = leaves%list(id)%p%abundances(chem_idx_some_spe%i_H2O)
  hc_params%X_OH    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_OH)
  hc_params%X_E     = leaves%list(id)%p%abundances(chem_idx_some_spe%i_E)
  hc_params%X_Hplus = leaves%list(id)%p%abundances(chem_idx_some_spe%i_Hplus)
  hc_params%X_gH    = leaves%list(id)%p%abundances(chem_idx_some_spe%i_gH)
  !
  hc_params%R_H2_form_rate = &
    get_H2_form_rate( &
      hc_params%R_H2_form_rate_coeff, &
      hc_params%X_gH, &
      hc_params%X_HI, &
      hc_params%n_gas)
  !
  leaves%list(id)%p%par%R_H2_form_rate = hc_params%R_H2_form_rate
end subroutine set_heatingcooling_params_from_cell


subroutine set_chemistry_params_from_cell(id)
  integer id
  chem_params => leaves%list(id)%p%par
end subroutine set_chemistry_params_from_cell


subroutine disk_set_a_cell_params(c, cell_params_copy)
  integer i
  type(type_cell), target :: c
  type(type_cell_rz_phy_basic), intent(in) :: cell_params_copy
  if (.not. associated(c%par)) then
    allocate(c%par)
  end if
  if (.not. allocated(c%h_c_rates)) then
    allocate(c%h_c_rates)
  end if
  !
  if (.not. allocated(c%abundances)) then
    allocate(c%abundances(chem_species%nSpecies), &
             c%col_den_toISM(chem_idx_some_spe%nItem), &
             c%col_den_toStar(chem_idx_some_spe%nItem))
  end if
  !
  c%iIter = 0
  c%quality = 0
  !
  c%par = cell_params_copy
  !
  c%par%rmin = c%xmin
  c%par%rmax = c%xmax
  c%par%rcen = (c%xmax + c%xmin) * 0.5D0
  c%par%dr   = c%xmax - c%xmin
  !
  c%par%zmin = c%ymin
  c%par%zmax = c%ymax
  c%par%zcen = (c%ymax + c%ymin) * 0.5D0
  c%par%dz   = c%ymax - c%ymin
  !
  c%par%daz  = 0D0
  !
  c%par%volume = phy_Pi * (c%par%rmax + c%par%rmin) * c%par%dr * c%par%dz * phy_AU2cm**3
  c%par%area_T = phy_Pi * (c%par%rmax + c%par%rmin) * c%par%dr * phy_AU2cm**2
  c%par%area_B = phy_Pi * (c%par%rmax + c%par%rmin) * c%par%dr * phy_AU2cm**2
  c%par%area_I = phy_2Pi * c%par%rmin * c%par%dz * phy_AU2cm**2
  c%par%area_O = phy_2Pi * c%par%rmax * c%par%dz * phy_AU2cm**2
  c%par%surf_area = c%par%area_T + c%par%area_B + c%par%area_I + c%par%area_O
  !
  ! Get gas number density and gas mass in each cell
  a_disk%andrews_gas%particlemass = c%par%MeanMolWeight * phy_mProton_CGS
  c%par%n_gas = Andrews_dens(c%par%rcen, c%par%zcen, a_disk%andrews_gas)
  c%par%mgas_cell = c%par%n_gas * c%par%volume * &
                    (phy_mProton_CGS * c%par%MeanMolWeight)
  !
  c%par%mdust_tot = 0D0
  c%par%ndust_tot = 0D0
  c%par%sigdust_ave = 0D0
  !
  c%par%ndustcompo = a_disk%ndustcompo
  !
  do i=1, a_disk%ndustcompo
    ! Dust mass density
    c%par%rho_dusts(i) = Andrews_dens(c%par%rcen, c%par%zcen, &
                                      a_disk%dustcompo(i)%andrews)
    c%par%mp_dusts(i)  = a_disk%dustcompo(i)%pmass_CGS ! Dust particle mass in gram
    c%par%n_dusts(i)  = c%par%rho_dusts(i) / c%par%mp_dusts(i)
    c%par%mdusts_cell(i)  = c%par%rho_dusts(i) * c%par%volume
    !
    c%par%mdust_tot = c%par%mdust_tot + c%par%mdusts_cell(i)
    c%par%ndust_tot = c%par%ndust_tot + c%par%n_dusts(i)
    c%par%sigdust_ave = c%par%sigdust_ave + &
      c%par%n_dusts(i) * phy_Pi * a_disk%dustcompo(i)%mrn%r2av * phy_micron2cm**2
  end do
  !
  c%par%sigdust_ave = c%par%sigdust_ave / c%par%ndust_tot
  !
  c%par%GrainRadius_CGS = sqrt(c%par%sigdust_ave / phy_Pi)
  !
  c%par%SitesPerGrain = 4D0 * c%par%sigdust_ave * SitesDensity_CGS
  !
  c%par%ratioDust2GasMass = c%par%mdust_tot / c%par%mgas_cell
  !
  c%par%ratioDust2HnucNum = c%par%ndust_tot / c%par%n_gas
  !
  c%par%dust_depletion = c%par%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
  !
  c%par%abso_wei = c%par%mdusts_cell / c%par%mdust_tot
  !
  ! Local dust size distribution; preliminary
  !c%mrn%rmin = c%par%aGrainMin_micron ! micron
  !c%mrn%rmax = c%par%aGrainMax_micron
  !c%mrn%n    = c%par%mrn_ind
  !call calc_dust_MRN_par(c%mrn)
  !c%par%GrainRadius_CGS = sqrt(c%mrn%r2av) * 1D-4 ! = sqrt(<r**2>)
  !
  ! Dust particle mass
  !c%par%mdust = 4.0D0*phy_Pi/3.0D0 * (1D-4)**3 * c%mrn%r3av * &
  !              c%par%GrainMaterialDensity_CGS
  ! Here n_dust is actually the mass density
  ! Two dust component: a background one, and a settled one.
  !c%par%n_dust = Andrews_dens(c%par%rcen, c%par%zcen, a_disk%andrews_dust) + &
  !               Andrews_dens(c%par%rcen, c%par%zcen, a_disk%andrews_dust_bg)
  !c%par%mdust_tot = c%par%n_dust * c%par%volume
  ! Now convert to particle density for dust
  !c%par%n_dust = c%par%n_dust / c%par%mdust
  !
  c%val(1) = c%par%n_gas
  !
  if (grid_config%use_data_file_input) then
    c%par%Tgas    = c%val(2)
    c%par%Tdust   = c%val(2)
  else
    c%par%Tgas    = 400D0 / (1D0 + c%par%rcen) * (1D0 + c%par%zcen)
    c%par%Tdust   = 0D0 ! instead of c%par%Tgas
  end if
  !
  c%par%pressure_thermal = 0D0
  c%par%gravity_z = 0D0
  c%par%gravity_acc_z = 0D0
  !
  ! c%par%UV_G0_factor = c%par%UV_G0_factor_background + &
  !   a_disk%UV_cont_phlumi_star_surface &
  !      / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
  !      / phy_Habing_photon_flux_CGS &
  !      * a_disk%geometric_factor_UV
  ! c%par%LymanAlpha_number_flux_0 = &
  !   a_disk%Lyman_phlumi_star_surface &
  !      / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2)
  ! c%par%LymanAlpha_energy_flux_0 = c%par%LymanAlpha_number_flux_0 * phy_LyAlpha_energy_CGS
  ! c%par%LymanAlpha_G0_factor = c%par%LymanAlpha_energy_flux_0 / phy_Habing_energy_flux_CGS
  c%par%Xray_flux_0 = &
    a_disk%Xray_phlumi_star_surface &
      / (4D0*phy_Pi * (c%par%rcen * phy_AU2cm)**2) &
      * a_disk%geometric_factor_Xray
  !
  ! Calculate the local velocity gradient, thermal velocity width, turbulent
  ! width, coherent length
  !
  call calc_local_dynamics(c)
  !
end subroutine disk_set_a_cell_params


subroutine calc_local_dynamics(c)
  type(type_cell), intent(inout) :: c
  associate( &
          G     => phy_GravitationConst_CGS, &
          M     => a_disk%star_mass_in_Msun * phy_Msun_CGS, &
          r     => c%par%rcen * phy_AU2cm, &
          !
          v     => c%par%velo_Kepler, &
          w     => c%par%omega_Kepler, &
          dv_dr => c%par%velo_gradient, &
          cs    => c%par%sound_speed, &
          delv  => c%par%velo_width_turb, &
          l     => c%par%coherent_length)
    v = sqrt(G * M / r)
    w = v / r
    dv_dr = 0.5D0 * v / r
    cs = sqrt(phy_kBoltzmann_CGS*c%par%Tgas / (phy_mProton_CGS * c%par%MeanMolWeight*2D0))
    delv = cs ! Todo
    l = delv / dv_dr
  end associate
end subroutine calc_local_dynamics


subroutine disk_set_gridcell_params
  integer i
  do i=1, leaves%nlen
    call disk_set_a_cell_params(leaves%list(i)%p, cell_params_ini)
  end do
end subroutine disk_set_gridcell_params


subroutine disk_set_disk_params
  !a_disk = disk_params_ini
  !
  ! Background dust
  a_disk%andrews_dust_bg = a_disk%andrews_gas
  a_disk%andrews_dust_bg%useNumDens = .false.
  a_disk%andrews_dust_bg%Md = &
    a_disk%andrews_gas%Md * a_disk%dust2gas_mass_bg
  !
  associate( &
    Lstar => a_disk%star_luminosity_in_Lsun * phy_Lsun_CGS, &
    uv2total => a_disk%ratio_uv2total, &
    lyman2uv => a_disk%ratio_lyman2uv, &
    xray2total => a_disk%ratio_xray2total)
    a_disk%UV_cont_phlumi_star_surface = &
      Lstar * uv2total * (1D0 - lyman2uv) / phy_UV_cont_energy_CGS
    a_disk%Lyman_phlumi_star_surface = &
      Lstar * uv2total * lyman2uv         / phy_LyAlpha_energy_CGS
    a_disk%Xray_phlumi_star_surface  = &
      Lstar * xray2total / (xray_energy_kev*1D3*phy_eV2erg)
    write(str_disp, '("!Stellar total luminosity = ", ES12.4, " erg s-1")') Lstar
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar UV cont luminosity = ", ES12.4, " erg s-1")') &
      a_disk%UV_cont_phlumi_star_surface * phy_UV_cont_energy_CGS
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar UV cont photon count rate = ", ES12.4, " s-1")') &
      a_disk%UV_cont_phlumi_star_surface
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar LyA luminosity = ", ES12.4, " erg s-1")') &
      a_disk%Lyman_phlumi_star_surface * phy_LyAlpha_energy_CGS
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar LyA photon count rate = ", ES12.4, " s-1")') &
      a_disk%Lyman_phlumi_star_surface
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar X-ray luminosity = ", ES12.4, " erg s-1")') &
      a_disk%Xray_phlumi_star_surface * (xray_energy_kev*1D3*phy_eV2erg)
    call display_string_both(str_disp, a_book_keeping%fU)
    write(str_disp, '("!Stellar X-ray photon count rate = ", ES12.4, " s-1")') &
      a_disk%Xray_phlumi_star_surface
    call display_string_both(str_disp, a_book_keeping%fU)
  end associate
end subroutine disk_set_disk_params


function get_local_doppler_kepler_scale(M, r, dv, factor)
  double precision :: get_local_doppler_kepler_scale
  double precision, intent(in) :: M, r, dv
  double precision, optional :: factor
  if (.NOT. present(factor)) then
    factor = 1D0
  end if
  get_local_doppler_kepler_scale = factor * 2D0 * r * dv / &
    sqrt(phy_GravitationConst_CGS * M * phy_Msun_CGS / (r * phy_AU2cm))
end function get_local_doppler_kepler_scale


subroutine calc_dust_MRN_par(mrn)
  type(type_dust_MRN), intent(inout) :: mrn
  double precision tmp1, tmp2, norm
  double precision, parameter :: smallnum = 1D-6
  !
  tmp1 = exp((1D0 - mrn%n) * log(mrn%rmin)) ! = rmin**(1-n)
  tmp2 = exp((1D0 - mrn%n) * log(mrn%rmax)) ! = rmax**(1-n)
  if (abs(mrn%n - 1D0) .le. smallnum) then
    norm = log(mrn%rmax/mrn%rmin)
  else
    norm = (tmp2 - tmp1) / (1D0 - mrn%n)
  end if
  if (abs(mrn%n - 2D0) .le. smallnum) then
    mrn%rav = log(mrn%rmax/mrn%rmin) / norm
  else
    mrn%rav  = (tmp2 * mrn%rmax    - tmp1 * mrn%rmin) &
        / ((2D0 - mrn%n) * norm)
  end if
  if (abs(mrn%n - 3D0) .le. smallnum) then
    mrn%r2av = log(mrn%rmax/mrn%rmin) / norm
  else
    mrn%r2av = (tmp2 * mrn%rmax**2 - tmp1 * mrn%rmin**2) &
        / ((3D0 - mrn%n) * norm)
  end if
  if (abs(mrn%n - 4D0) .le. smallnum) then
    mrn%r3av = log(mrn%rmax/mrn%rmin) / norm
  else
    mrn%r3av = (tmp2 * mrn%rmax**3 - tmp1 * mrn%rmin**3) &
        / ((4D0 - mrn%n) * norm)
  end if
end subroutine calc_dust_MRN_par


subroutine save_chem_rates(i0)
  integer, intent(in) :: i0
  integer fU, k
  character(len=128) filename, dir
  type(type_heating_cooling_parameters) heat_cool_log
  ! Use namelist for output some logging infomation.
  ! Not very readable, but easy to implement.
  namelist /cell_par_log/ heat_cool_log
  !
  write(filename, '("reac_rates_cell_", I4.4, ".dat")') i0
  if(.NOT. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit for output!'
    stop
  end if
  dir = trim(combine_dir_filename(a_book_keeping%dir, 'rates_log/'))
  if (.NOT. dir_exist(dir)) then
    call my_mkdir(dir)
  end if
  call openFileSequentialWrite(fU, combine_dir_filename(dir, filename), 99999)
  !
  heat_cool_log = hc_params
  write(fU, nml=cell_par_log)
  !
  do k=1, chem_net%nReactions
    write(fU, '(A135, ES16.4E4)') chem_reac_str%list(k), chem_net%rates(k)
  end do
  close(fU)
end subroutine save_chem_rates


subroutine save_post_config_params
  type(type_disk_basic_info) disk_params_tmp
  namelist /disk_params_log/ disk_params_tmp
  disk_params_tmp = a_disk
  if (FileUnitOpened(a_book_keeping%fU)) then
    write(a_book_keeping%fU, nml=disk_params_log)
    flush(a_book_keeping%fU)
  end if
end subroutine save_post_config_params


subroutine load_refine_check_species
  integer fU, i, i1, ios, n
  character(len=const_len_init_abun_file_row) str
  n = GetFileLen_comment_blank( &
      combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
        a_disk_iter_params%filename_list_check_refine), '!')
  allocate(idx_Species_check_refine(n), &
           thr_Species_check_refine(n))
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialRead(fU, &
    combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
      a_disk_iter_params%filename_list_check_refine), 99)
  i1 = 0
  do
    read(fU, FMT='(A)', IOSTAT=ios) str
    if (ios .NE. 0) then
      exit
    end if
    do i=1, chem_species%nSpecies
      if (trim(str(1:const_len_species_name)) .EQ. chem_species%names(i)) then
        i1 = i1 + 1
        idx_Species_check_refine(i1) = i
        read(str(const_len_species_name+1:const_len_init_abun_file_row), &
          '(F7.1)') thr_Species_check_refine(i1)
        exit
      end if
    end do
  end do
  close(fU)
  a_disk_iter_params%nSpecies_check_refine = i1
  write(str_disp, '("! Species used for checking refinement:")')
  call display_string_both(str_disp, a_book_keeping%fU)
  do i=1, a_disk_iter_params%nSpecies_check_refine
    write(str_disp, '("! ", A12, ES12.2)') &
      chem_species%names(idx_Species_check_refine(i)), &
      thr_Species_check_refine(i)
    call display_string_both(str_disp, a_book_keeping%fU)
  end do
end subroutine load_refine_check_species


subroutine do_refine
  integer i, n_refine
  a_disk_iter_params%ncell_refine = 0
  do i=1, leaves%nlen
    if (need_to_refine(leaves%list(i)%p, n_refine)) then
      a_disk_iter_params%ncell_refine = a_disk_iter_params%ncell_refine + 1
      write(str_disp, '("!", I4, A, 4ES12.2, " into ", I4, " parts.")') &
        a_disk_iter_params%ncell_refine, ' Refining ', &
        leaves%list(i)%p%xmin, leaves%list(i)%p%xmax, &
        leaves%list(i)%p%ymin, leaves%list(i)%p%ymax, n_refine
      call display_string_both(str_disp, a_book_keeping%fU)
      call refine_this_cell_vertical(leaves%list(i)%p, n_refine)
    end if
  end do
end subroutine do_refine


subroutine remake_index
  call get_number_of_leaves(root)
  leaves%nlen = root%nleaves
  call grid_make_leaves(root)
  call grid_make_neighbors
  call grid_make_surf_bott
end subroutine remake_index


function need_to_refine(c, n_refine)
  logical need_to_refine
  type(type_cell), target :: c
  integer, intent(out), optional :: n_refine
  integer i, i0, i1, j
  double precision val_max, val_min
  logical flag1, flag2
  flag1 = .false.
  flag2 = .false.
  if (present(n_refine)) then
    n_refine = 0
  end if
  if (c%par%dz .le. grid_config%smallest_cell_size) then
    need_to_refine = .false.
    return
  end if
  do i=1, c%above%n
    i0 = c%above%idx(i)
    !if (leaves%list(i0)%p%iIter .lt. c%iIter) then
    if (leaves%list(i0)%p%iIter .lt. 1) then
      cycle
    end if
    do j=1, a_disk_iter_params%nSpecies_check_refine
      i1 = idx_Species_check_refine(j)
      val_max = max(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      val_min = min(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      if (val_max .gt. thr_Species_check_refine(j)) then
        if (val_max / val_min .gt. a_disk_iter_params%threshold_ratio_refine) then
          flag1 = .true.
          if (present(n_refine)) then
            n_refine = max(n_refine, int(log10(val_max / val_min)) * 2)
          end if
        end if
      end if
    end do
  end do
  do i=1, c%below%n
    i0 = c%below%idx(i)
    !if (leaves%list(i0)%p%iIter .lt. c%iIter) then
    if (leaves%list(i0)%p%iIter .lt. 1) then
      cycle
    end if
    do j=1, a_disk_iter_params%nSpecies_check_refine
      i1 = idx_Species_check_refine(j)
      val_max = max(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      val_min = min(leaves%list(i0)%p%abundances(i1), c%abundances(i1))
      if (val_max .gt. thr_Species_check_refine(j)) then
        if (val_max / val_min .gt. a_disk_iter_params%threshold_ratio_refine) then
          flag2 = .true.
          if (present(n_refine)) then
            n_refine = max(n_refine, int(log10(val_max / val_min)) * 2)
          end if
        end if
      end if
    end do
  end do
  need_to_refine = flag1 .or. flag2
  return
end function need_to_refine


subroutine refine_this_cell_vertical(c, n)
  ! c is a working cell that needs to be refined.
  type(type_cell), target :: c
  double precision dy
  integer, intent(in), optional :: n
  integer i, ndivide
  !
  if (present(n)) then
    ndivide = n
  else
    ndivide = 3
  end if
  !
  if (ndivide .lt. 2) then
    return
  end if
  !
  c%nleaves = ndivide
  c%nChildren = ndivide
  call init_children(c, ndivide)
  !
  dy = (c%ymax - c%ymin) / dble(ndivide)
  !
  do i=1, c%nChildren
    associate(cc => c%children(i)%p)
      cc%xmin = c%xmin
      cc%xmax = c%xmax
      cc%ymin = c%ymin + dble(i-1) * dy
      cc%ymax = c%ymin + dble(i)   * dy
      !
      ! Re-interpolate density from the input data.
      call set_cell_par_preliminary(cc)
      cc%using = .true.
      cc%converged = .false.
      cc%nOffspring = 0
      cc%nChildren = 0
      cc%nleaves = 1
      !
      cc%iIter = c%iIter
      !
      call disk_set_a_cell_params(cc, c%par)
      cc%par%Tgas = c%par%Tgas
      !
      cc%h_c_rates = c%h_c_rates
      cc%abundances = c%abundances
      !cc%col_den = c%col_den
      !cc%col_den_acc = c%col_den_acc
    end associate
  end do
  ! Avoid numerical roundings
  c%children(1)%p%ymin       = c%ymin
  c%children(ndivide)%p%ymax = c%ymax
  do i=1, c%nChildren-1
    c%children(i)%p%ymax = c%children(i+1)%p%ymin
  end do
  !
  ! Deactivate c
  c%using = .false.
  c%converged = .false.
  !deallocate(c%par, c%h_c_rates, c%abundances, c%col_den, c%col_den_acc)
  !deallocate(c%inner%idx, c%inner%fra)
  !deallocate(c%outer%idx, c%outer%fra)
  !deallocate(c%above%idx, c%above%fra)
  !deallocate(c%below%idx, c%below%fra)
  !deallocate(c%around%idx, c%around%fra)
  !deallocate(c%inner, c%outer, c%above, c%below, c%around)
end subroutine refine_this_cell_vertical


subroutine disk_iteration_postproc
  integer fU, fU1, fU2, ios, i, i0, j, idx, idx_diff
  double precision r, z
  double precision sum_prod, sum_dest, accum
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  call get_species_produ_destr
  !
  write(*,*) 'Trying to find out where are the elements.'
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_points), 99)
  if (.not. getFileUnit(fU1)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialWrite(fU1, &
       combine_dir_filename(a_book_keeping%dir, &
         a_disk_ana_params%file_analyse_res_ele), 999)
  if (.not. getFileUnit(fU2)) then
    write(*,*) 'Cannot get a file unit in disk_iteration_postproc.'
    return
  end if
  call openFileSequentialWrite(fU2, &
       combine_dir_filename(a_book_keeping%dir, &
         a_disk_ana_params%file_analyse_res_contri), 999)
  do
    read(fU, '(2F6.2)', iostat=ios) r, z
    if (ios .ne. 0) then
      exit
    end if
    idx = 0
    do i=1, leaves%nlen
      if ((leaves%list(i)%p%par%rmin .le. r) .and. (leaves%list(i)%p%par%rmax .ge. r) .and. &
          (leaves%list(i)%p%par%zmin .le. z) .and. (leaves%list(i)%p%par%zmax .ge. z)) then
        idx = i
        exit
      end if
    end do
    if (idx .eq. 0) then
      write(*, '("Point (", 2F6.2, ")", A)') r, z, ' not in any cells!'
      cycle
    end if
    chemsol_stor%y(1:chem_species%nSpecies) = &
        leaves%list(idx)%p%abundances(1:chem_species%nSpecies)
    call chem_elemental_residence
    write(fU1, '("(", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') r, z, &
      leaves%list(idx)%p%par%Tgas, leaves%list(idx)%p%par%Tdust, &
      leaves%list(idx)%p%par%n_gas, leaves%list(idx)%p%par%Ncol, &
      leaves%list(idx)%p%par%Av
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    if (leaves%list(idx)%p%above%n .gt. 0) then
      idx_diff = leaves%list(idx)%p%above%idx(1)
    else if (leaves%list(idx)%p%below%n .gt. 0) then
      idx_diff = leaves%list(idx)%p%below%idx(1)
    else
      idx_diff = idx
    end if
    !
    call set_chemistry_params_from_cell(idx_diff)
    call chem_cal_rates
    call get_contribution_each
    !
    write(fU2, '("This (", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') r, z, &
      leaves%list(idx)%p%par%Tgas, leaves%list(idx)%p%par%Tdust, &
      leaves%list(idx)%p%par%n_gas, leaves%list(idx)%p%par%Ncol, &
      leaves%list(idx)%p%par%Av
    write(fU2, '("Diff (", 2F6.2, ")", 2F7.1, 2ES12.2, F9.1)') &
      leaves%list(idx_diff)%p%par%rcen, &
      leaves%list(idx_diff)%p%par%zcen, &
      leaves%list(idx_diff)%p%par%Tgas,  leaves%list(idx_diff)%p%par%Tdust, &
      leaves%list(idx_diff)%p%par%n_gas, leaves%list(idx_diff)%p%par%Ncol, &
      leaves%list(idx_diff)%p%par%Av
    do i=1, chem_species%nSpecies
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(A12, ": ", ES12.2, " Diff: ", ES12.2, " Rate: ", ES12.2)') chem_species%names(i), &
        chemsol_stor%y(i), leaves%list(idx_diff)%p%abundances(i), &
        sum_prod - sum_dest
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%produ(i)%contri(j) .le. &
            chem_species%produ(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
      write(fU2, '(2X, A, 2X, ES12.2)') 'Destruction', sum_dest
      accum = 0D0
      do j=1, min(chem_species%destr(i)%nItem, 20)
        i0 = chem_species%destr(i)%list(j)
        accum = accum + chem_species%destr(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
  end do
  close(fU)
  close(fU1)
  close(fU2)
end subroutine disk_iteration_postproc




subroutine load_ana_species_list
  integer fU, ios, i, n
  integer, dimension(:), allocatable :: list_tmp
  character(len=12) str
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in load_ana_species_list.'
    return
  end if
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_species), 99)
  allocate(list_tmp(chem_species%nSpecies))
  n = 0
  do
    read(fU, '(A12)', iostat=ios) str
    if (ios .ne. 0) then
      exit
    end if
    do i=1, chem_species%nSpecies
      if (chem_species%names(i) .eq. str) then
        !
        if (.not. is_in_list_int(i, n, list_tmp(1:n))) then
          n = n + 1
          list_tmp(n) = i
        end if
        !
        exit
        !
      end if
    end do
  end do
  close(fU)
  !
  ana_splist%nlen = n
  if (n .gt. 0) then
    if (allocated(ana_splist%vals)) then
      deallocate(ana_splist%vals)
    end if
    allocate(ana_splist%vals(n))
    ana_splist%vals = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine load_ana_species_list



subroutine load_ana_points_list
  integer fU, ios, i, n
  double precision r, z
  integer, dimension(:), allocatable :: list_tmp
  if (.not. a_disk_ana_params%do_analyse) then
    return
  end if
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit in load_ana_points_list.'
    return
  end if
  call openFileSequentialRead(fU, &
       combine_dir_filename(a_disk_ana_params%analyse_points_inp_dir, &
         a_disk_ana_params%file_list_analyse_points), 99)
  allocate(list_tmp(leaves%nlen))
  n = 0
  do
    read(fU, '(2F6.2)', iostat=ios) r, z
    if (ios .ne. 0) then
      exit
    end if
    do i=1, leaves%nlen
      if ((leaves%list(i)%p%par%rmin .le. r) .and. (leaves%list(i)%p%par%rmax .ge. r) .and. &
          (leaves%list(i)%p%par%zmin .le. z) .and. (leaves%list(i)%p%par%zmax .ge. z)) then
        !
        if (.not. is_in_list_int(i, n, list_tmp(1:n))) then
          n = n + 1
          list_tmp(n) = i
        end if
        !
        exit
        !
      end if
    end do
  end do
  close(fU)
  !
  ana_ptlist%nlen = n
  if (n .gt. 0) then
    if (allocated(ana_ptlist%vals)) then
      deallocate(ana_ptlist%vals)
    end if
    allocate(ana_ptlist%vals(n))
    ana_ptlist%vals = list_tmp(1:n)
  end if
  deallocate(list_tmp)
end subroutine load_ana_points_list



subroutine chem_analyse(id)
  integer, intent(in) :: id
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=128) fname_pre
  character(len=32) FMTstryHistory
  double precision dy_y, dt_t
  double precision frac
  frac = 0.1D0
  !
  write(*, '(/A/)') 'Doing some analysis... Might be very slow.'
  !
  if (.not. getFileUnit(fU3)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  write(fname_pre, &
        '(I4.4, "_rz_", F0.6, "_", F0.6, "_iter_", I3.3)') &
        id, leaves%list(id)%p%xmin, leaves%list(id)%p%ymin, &
        leaves%list(id)%p%iIter
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'evol_'//trim(fname_pre)//'.dat'), 999999)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 2
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, '  Tgas        '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 2
  do i=1, chemsol_params%n_record
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i)
  end do
  close(fU3)
  !
  if (.not. getFileUnit(fU1)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  !
  call openFileSequentialWrite(fU1, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'ele_'//trim(fname_pre)//'.dat'), 999)
  !
  if (.not. getFileUnit(fU2)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  call openFileSequentialWrite(fU2, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'contri_'//trim(fname_pre)//'.dat'), 999)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = 1+chemsol_params%n_record/20
  end if
  !
  write(fU1, '(2F10.1, 2ES12.2, F9.1, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, chem_params%Ncol, &
    chem_params%Av, &
    id, leaves%list(id)%p%iIter, &
    leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
    leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
  write(fU2, '(2F10.1, 2ES12.2, F9.1, 2I5, 4ES16.6)') &
    chem_params%Tgas,  chem_params%Tdust, &
    chem_params%n_gas, chem_params%Ncol, &
    chem_params%Av, &
    id, leaves%list(id)%p%iIter, &
    leaves%list(id)%p%xmin, leaves%list(id)%p%xmax, &
    leaves%list(id)%p%ymin, leaves%list(id)%p%ymax
  do k=1, chemsol_params%n_record, a_disk_ana_params%ana_i_incr
    !
    if (k .ge. 2) then
      dy_y = maxval( &
        abs((chemsol_stor%record(:, k) - chemsol_stor%record(:, k-1))) / &
        (chemsol_stor%record(:, k) + chemsol_stor%record(:, k-1) + 1D-15))
      dt_t = (chemsol_stor%touts(k) - chemsol_stor%touts(k-1)) / &
             (chemsol_stor%touts(k) + chemsol_stor%touts(k-1))
      if (dy_y .lt. frac * dt_t) then
        cycle
      end if
    end if
    !
    write(fU1, '("time = ", ES14.4)') chemsol_stor%touts(k)
    !
    chemsol_stor%y = chemsol_stor%record(:, k)
    !
    call chem_elemental_residence
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("time = ", ES14.4)') chemsol_stor%touts(k)
    if (ana_splist%nlen .le. 0) then
      cycle
    end if
    call get_contribution_each
    do i=1, chem_species%nSpecies
      if (.not. is_in_list_int(i, ana_splist%nlen, ana_splist%vals)) then
        cycle
      end if
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chemsol_stor%y(i)
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%produ(i)%contri(j) .le. &
            chem_species%produ(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
      write(fU2, '(2X, A, 2X, ES12.2)') 'Destruction', sum_dest
      accum = 0D0
      do j=1, min(chem_species%destr(i)%nItem, 20)
        i0 = chem_species%destr(i)%list(j)
        accum = accum + chem_species%destr(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
  end do
  close(fU1)
  close(fU2)
end subroutine chem_analyse




subroutine a_test_case
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=64) FMTstryHistory, fname_pre
  !
  allocate(chem_params)
  chem_params = a_disk_ana_params%chempar
  !
  call chem_read_reactions
  call chem_load_reactions
  call chem_parse_reactions
  call chem_get_dupli_reactions
  call chem_get_idx_for_special_species
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  call chem_load_initial_abundances
  !chemsol_stor%y(chem_idx_some_spe%i_E) = &
  !  chemsol_stor%y(chem_idx_some_spe%i_E) + &
  !    sum(chemsol_stor%y(1:chem_species%nSpecies) * &
  !        dble(chem_species%elements(1,:)))
  write(*, '(A)') 'Initial abundances:'
  do i=1, chem_species%nSpecies
    if (chemsol_stor%y(i) .ge. 1D-20) then
      write(*, '(I5, 2X, A, ES12.5)') i, chem_species%names(i), &
        chemsol_stor%y(i)
    end if
  end do
  write(*, '(//)')
  !
  associate(ch => chem_params)
    !ch%Tgas = 10D0
    !ch%Tdust = 10D0
    !ch%n_gas = 1D6
    !ch%UV_G0_factor = 0D0
    !ch%UV_G0_factor_background = 1D0
    !ch%Av = 10D0
    !ch%LymanAlpha_number_flux_0 = 0D0
    !ch%Xray_flux_0 = 0D0
    !ch%Ncol = 1D22
    !ch%dNcol = 1D21
    !ch%f_selfshielding_H2 = 0D0
    !ch%f_selfshielding_CO = 0D0
    !ch%f_selfshielding_H2O = 0D0
    !ch%f_selfshielding_OH = 0D0
    ch%GrainMaterialDensity_CGS = 2D0
    ch%ratioDust2GasMass = 0.01D0
    ch%MeanMolWeight = 1.4D0
    ch%ratioDust2HnucNum = &
          ch%ratioDust2GasMass * (phy_mProton_CGS * ch%MeanMolWeight) &
          / (4.0D0*phy_Pi/3.0D0 * (ch%GrainRadius_CGS)**3 * &
             ch%GrainMaterialDensity_CGS)
    ch%dust_depletion = ch%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
    ch%ndust_tot = ch%n_gas * ch%ratioDust2HnucNum
    chemsol_stor%y(chem_species%nSpecies+1) = ch%Tgas
  end associate
  !
  call chem_cal_rates
  call chem_set_solver_flags
  chemsol_params%evolT = .false.
  call chem_evol_solve
  !
  write(*,*) 'Doing some analysis... Might be very slow.'
  !
  call get_species_produ_destr
  !
  if (.not. getFileUnit(fU3)) then
    write(*,*) 'Cannot get a file unit.'
    return
  end if
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_iter_params%iter_files_dir, 'func_of_time.dat'), 999999)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 2
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, &
    '  Tgas        '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 2
  do i=1, chemsol_params%n_record
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i)
  end do
  close(fU3)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = chemsol_params%n_record / 4
  end if
  do k=1, chemsol_params%n_record, a_disk_ana_params%ana_i_incr
    write(fname_pre, '(I4.4, "_")') k
    !
    if (.not. getFileUnit(fU1)) then
      write(*,*) 'Cannot get a file unit.'
      return
    end if
    call openFileSequentialWrite(fU1, &
      combine_dir_filename(a_disk_iter_params%iter_files_dir, &
        trim(fname_pre)//'elemental_residence.dat'), 999)
    !
    if (.not. getFileUnit(fU2)) then
      write(*,*) 'Cannot get a file unit.'
      return
    end if
    call openFileSequentialWrite(fU2, &
      combine_dir_filename(a_disk_iter_params%iter_files_dir, &
        trim(fname_pre)//'contribution_reactions.dat'), 999)
    !
    chemsol_stor%y = chemsol_stor%record(:, k)
    !
    call chem_elemental_residence
    write(fU1, '(ES12.2, 2F7.1, 2ES12.2, F9.1)') &
      chemsol_stor%touts(k), &
      chem_params%Tgas,  chem_params%Tdust, &
      chem_params%n_gas, chem_params%Ncol, &
      chem_params%Av
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    call get_contribution_each
    !
    write(fU2, '(ES12.2, 2F7.1, 2ES12.2, F9.1)') &
      chemsol_stor%touts(k), &
      chem_params%Tgas,  chem_params%Tdust, &
      chem_params%n_gas, chem_params%Ncol, &
      chem_params%Av
    do i=1, chem_species%nSpecies
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chemsol_stor%y(i)
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, 3ES12.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%produ(i)%contri(j) .le. &
            chem_species%produ(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
      write(fU2, '(2X, A, 2X, ES12.2)') 'Destruction', sum_dest
      accum = 0D0
      do j=1, min(chem_species%destr(i)%nItem, 20)
        i0 = chem_species%destr(i)%list(j)
        accum = accum + chem_species%destr(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, 3ES11.2, 2F8.1)') &
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
    close(fU1)
    close(fU2)
  end do
end subroutine a_test_case



subroutine b_test_case
  integer i, j, fU
  type(type_cell_rz_phy_basic), pointer :: ch => null()
  double precision Tmin, Tmax, dT, ratio
  double precision n_gas_min, n_gas_max, dn
  double precision h_c_net_rate
  character(len=128) filename, fname_pre, header
  type(type_cell), pointer :: c => null()
  !
  filename = 'Tgas_hc_abundances.dat'
  !
  allocate(ch)
  ch = a_disk_ana_params%chempar
  chem_params => ch
  !
  call chem_read_reactions
  call chem_load_reactions
  call chem_parse_reactions
  call chem_get_dupli_reactions
  call chem_get_idx_for_special_species
  call load_species_enthalpies
  call get_reaction_heat
  !
  call chem_make_sparse_structure
  call chem_prepare_solver_storage
  call chem_evol_solve_prepare
  !
  call chem_load_initial_abundances
  !
  call heating_cooling_prepare
  !
  call load_ana_species_list
  call get_species_produ_destr
  !
  Tmin = 1D2
  Tmax = 1200D0
  dT = 10D0
  n_gas_min = 3.5D5
  n_gas_max = 3.6D5
  dn = 1D3
  ratio = 1D0
  !
  ch%Tgas = Tmin
  ch%n_gas = n_gas_min
  !
  allocate(c)
  allocate(&!c%col_den_acc(chem_idx_some_spe%nItem), &
           !c%col_den(chem_idx_some_spe%nItem), &
           c%abundances(chem_species%nSpecies))
  allocate(c%around, c%above, c%below, c%inner, c%outer)
  allocate(c%h_c_rates)
  !
  c%par => ch
  !
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit!'
    stop
  end if
  call openFileSequentialWrite(fU, &
    combine_dir_filename(a_disk_iter_params%iter_files_dir, filename), 99999)
  !
  call write_header(fU)
  !
  a_disk_ana_params%analyse_out_dir = &
    trim(combine_dir_filename(a_disk_iter_params%iter_files_dir, 'ana/'))
  if (.not. dir_exist(a_disk_ana_params%analyse_out_dir)) then
    call my_mkdir(a_disk_ana_params%analyse_out_dir)
  end if
  !
  do i=1, 1
    do j=1, 299
      !
      chemsol_params%evolT = .true.
      !
      ch%Tgas = ch%Tdust
      !
      ch%GrainMaterialDensity_CGS = 2D0
      ch%ratioDust2GasMass = 0.01D0
      ch%MeanMolWeight = 1.4D0
      ch%ratioDust2HnucNum = &
            ch%ratioDust2GasMass * (phy_mProton_CGS * ch%MeanMolWeight) &
            / (4.0D0*phy_Pi/3.0D0 * (ch%GrainRadius_CGS)**3 * &
               ch%GrainMaterialDensity_CGS)
      ch%dust_depletion = ch%ratioDust2GasMass / phy_ratioDust2GasMass_ISM
      ch%ndust_tot = ch%n_gas * ch%ratioDust2HnucNum
      write(*,*) 'Dust density ', ch%ndust_tot
      write(*,*) ch%ratioDust2HnucNum
      write(*,*) ch%n_gas
      !
      ch%velo_Kepler = 30D5
      ch%omega_Kepler = ch%velo_Kepler / phy_AU2cm
      ch%velo_gradient = 0.5D0 * ch%velo_Kepler / phy_AU2cm
      ch%velo_width_turb = ch%velo_Kepler
      ch%coherent_length = ch%velo_width_turb / ch%velo_gradient
      !
      write(*,'(I4, F9.1, ES12.4, F9.1/)') i, ch%Tgas, ch%n_gas, ch%Tdust
      !
      chemsol_stor%y(1:chem_species%nSpecies) = chemsol_stor%y0(1:chem_species%nSpecies)
      chemsol_stor%y(chem_species%nSpecies+1) = ch%Tgas
      !
      call chem_cal_rates
      write(*,'(2ES12.2/)') chem_params%f_selfshielding_H2, chem_params%Av
      call chem_set_solver_flags_alt(1)
      !
      hc_params%type_cell_rz_phy_basic = ch
      !
      hc_params%Neufeld_dv_dz = 10D0/phy_AU2cm
      hc_params%Neufeld_G     = 1D0
      !
      hc_params%X_H2    = chemsol_stor%y(chem_idx_some_spe%i_H2)
      hc_params%X_HI    = chemsol_stor%y(chem_idx_some_spe%i_HI)
      hc_params%X_CI    = chemsol_stor%y(chem_idx_some_spe%i_CI)
      hc_params%X_Cplus = chemsol_stor%y(chem_idx_some_spe%i_Cplus)
      hc_params%X_OI    = chemsol_stor%y(chem_idx_some_spe%i_OI)
      hc_params%X_CO    = chemsol_stor%y(chem_idx_some_spe%i_CO)
      hc_params%X_H2O   = chemsol_stor%y(chem_idx_some_spe%i_H2O)
      hc_params%X_OH    = chemsol_stor%y(chem_idx_some_spe%i_OH)
      hc_params%X_E     = chemsol_stor%y(chem_idx_some_spe%i_E)
      hc_params%X_Hplus = chemsol_stor%y(chem_idx_some_spe%i_Hplus)
      hc_params%X_gH    = chemsol_stor%y(chem_idx_some_spe%i_gH)
      !
      hc_params%R_H2_form_rate = &
        get_H2_form_rate( &
          hc_params%R_H2_form_rate_coeff, &
          hc_params%X_gH, &
          hc_params%X_HI, &
          hc_params%n_gas)
      ch%R_H2_form_rate = hc_params%R_H2_form_rate
      !
      !call realtime_heating_cooling_rate(tmp, chemsol_params%NEQ, chemsol_stor%y)
      !write(*,'(2ES16.6/)') tmp, heating_minus_cooling()
      !call disp_h_c_rates
      call chem_evol_solve
      !
      c%abundances  = chemsol_stor%y(1:chem_species%nSpecies)
      !c%col_den     = c%abundances(chem_idx_some_spe%idx) * c%par%dNcol
      !c%col_den_acc = c%abundances(chem_idx_some_spe%idx) * c%par%Ncol
      !
      hc_Tgas = ch%Tgas
      hc_Tdust = ch%Tdust
      h_c_net_rate = heating_minus_cooling()
      !
      c%h_c_rates = heating_cooling_rates
      c%par%t_final = chemsol_stor%touts(chemsol_params%n_record_real)
      !
      call disk_save_results_write(fU, c)
      !
      write(fname_pre, '(I4.4, "_", I4.4)') i, j
      write(header, '("n_gas = ", ES13.6)') ch%n_gas
      !
      a_disk_ana_params%ana_i_incr = 1
      call do_a_analysis(fname_pre, header)
      !
      ch%n_gas = ch%n_gas + dn
      dn = dn * ratio
      if (ch%n_gas .gt. n_gas_max) then
        exit
      end if
    end do
    ch%Tgas = ch%Tgas + dT
    dT = dT * ratio
    if (ch%Tgas .gt. Tmax) then
      exit
    end if
  end do
  close(fU)
  !
end subroutine b_test_case



function get_H2_form_rate(c, XgH, XH, ngas) result(r)
  ! dn(H2)/dt
  double precision r
  double precision, intent(in) :: c, XgH, XH, ngas
  if (chemsol_params%H2_form_use_moeq) then
    r = c * XgH * XH * ngas
  else
    r = c * XgH * XgH * ngas
    !r = c * XH * ngas
  end if
end function get_H2_form_rate



subroutine do_a_analysis(fname_pre, header)
  integer i, j, k, i0, fU1, fU2, fU3
  double precision sum_prod, sum_dest, accum
  character(len=128), intent(in) :: fname_pre, header
  character(len=32) FMTstryHistory
  double precision dy_y, dt_t
  double precision frac
  double precision r
  frac = 0.1D0
  !
  write(*, '(/A/)') 'Doing some analysis... Might be slow.'
  !
  call openFileSequentialWrite(fU3, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'evol_'//trim(fname_pre)//'.dat'), maxRowLen=999999, getu=1)
  !
  write(FMTstryHistory, '("(", I4, "A14)")') chem_species%nSpecies + 3
  write(fU3, FMTstryHistory) '!Time_(yr)    ', chem_species%names, &
    '  Tgas        ', &
    '  hc          '
  write(FMTstryHistory, '("(", I4, "ES14.4E4)")') chem_species%nSpecies + 3
  do i=1, chemsol_params%n_record
    call realtime_heating_cooling_rate(r, chemsol_params%NEQ, chemsol_stor%record(:, i))
    write(fU3, FMTstryHistory) chemsol_stor%touts(i), chemsol_stor%record(:, i), r
  end do
  close(fU3)
  !
  call openFileSequentialWrite(fU1, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'ele_'//trim(fname_pre)//'.dat'), maxRowLen=999, getu=1)
  call openFileSequentialWrite(fU2, &
    combine_dir_filename(a_disk_ana_params%analyse_out_dir, &
      'contri_'//trim(fname_pre)//'.dat'), maxRowLen=999, getu=1)
  !
  if (a_disk_ana_params%ana_i_incr .le. 0) then
    a_disk_ana_params%ana_i_incr = 1+chemsol_params%n_record/20
  end if
  !
  write(fU1, '(A)') trim(header)
  write(fU2, '(A)') trim(header)
  !
  do k=1, chemsol_params%n_record, a_disk_ana_params%ana_i_incr
    !+++
    if ((chemsol_stor%touts(k) .le. 1D2) .or. (chemsol_stor%touts(k) .ge. 1D5)) then
      cycle
    end if
    if (chemsol_stor%touts(k) .ge. 4D2) then
      if (mod(k, 50) .ne. 0) then
        cycle
      end if
    end if
    !---
    if (k .ge. 2) then
      dy_y = maxval( &
        abs((chemsol_stor%record(1:chem_species%nSpecies, k) - &
             chemsol_stor%record(1:chem_species%nSpecies, k-1))) / &
        (chemsol_stor%record(1:chem_species%nSpecies, k) + &
         chemsol_stor%record(1:chem_species%nSpecies, k-1) + 1D-15))
      dt_t = (chemsol_stor%touts(k) - chemsol_stor%touts(k-1)) / &
             (chemsol_stor%touts(k) + chemsol_stor%touts(k-1))
      if (dy_y .lt. frac * dt_t) then
        cycle
      end if
    end if
    !
    write(fU1, '("time = ", ES14.4)') chemsol_stor%touts(k)
    !
    chemsol_stor%y(1:chem_species%nSpecies) = chemsol_stor%record(1:chem_species%nSpecies, k)
    chem_params%Tgas = chemsol_stor%record(chem_species%nSpecies+1, k)
    call chem_cal_rates
    !
    call chem_elemental_residence
    write(fU1, '(4X, "Total net charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * dble(chem_species%elements(1,:)))
    write(fU1, '(4X, "Total free charge: ", ES10.2)') &
        sum(chemsol_stor%y(1:chem_species%nSpecies) * abs(dble(chem_species%elements(1,:)))) / 2D0
    do i=1, const_nElement
      write(fU1, '(4X, A8)') const_nameElements(i)
      do j=1, chem_ele_resi(i)%n_nonzero
        i0 = chem_ele_resi(i)%iSpecies(j)
        write(fU1, '(6X, A12, 3ES10.2)') chem_species%names(i0), chemsol_stor%y(i0), &
          chem_ele_resi(i)%ele_frac(j), chem_ele_resi(i)%ele_accu(j)
      end do
    end do
    !
    write(fU2, '("time = ", ES14.4)') chemsol_stor%touts(k)
    if (ana_splist%nlen .le. 0) then
      cycle
    end if
    !
    call get_contribution_each
    !
    do i=1, chem_species%nSpecies
      if (.not. is_in_list_int(i, ana_splist%nlen, ana_splist%vals)) then
        cycle
      end if
      write(fU2, '(A12, ES12.2)') chem_species%names(i), chemsol_stor%y(i)
      sum_prod = sum(chem_species%produ(i)%contri)
      sum_dest = sum(chem_species%destr(i)%contri)
      write(fU2, '(2X, A, 2X, ES12.2)') 'Production', sum_prod
      accum = 0D0
      do j=1, min(chem_species%produ(i)%nItem, 20)
        i0 = chem_species%produ(i)%list(j)
        accum = accum + chem_species%produ(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%produ(i)%contri(j), accum, accum/sum_prod, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%produ(i)%contri(j) .le. &
            chem_species%produ(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
      write(fU2, '(2X, A, 2X, ES12.2)') 'Destruction', sum_dest
      accum = 0D0
      do j=1, min(chem_species%destr(i)%nItem, 20)
        i0 = chem_species%destr(i)%list(j)
        accum = accum + chem_species%destr(i)%contri(j)
        write(fU2, '(4X, I4, 2ES12.2, F8.2, ES12.2, 2X, 6A12, ES12.2, 2F9.2, 2F8.1)') &
          j, chem_species%destr(i)%contri(j), accum, accum/sum_dest, chem_net%rates(i0), &
          chem_net%reac_names(1:2, i0), chem_net%prod_names(1:4, i0), &
          chem_net%ABC(1:3, i0), chem_net%T_range(1:2, i0)
        if (chem_species%destr(i)%contri(j) .le. &
            chem_species%destr(i)%contri(1) * 1D-6) then
          exit
        end if
      end do
    end do
  end do
  close(fU1)
  close(fU2)
end subroutine do_a_analysis



!subroutine save_fits_cube(filename, im)
!  character(len=*), intent(in) :: filename
!  type(type_image), intent(in) :: im
!  integer stat, fU, blocksize, bitpix, naxis
!  integer, dimension(3) :: naxes
!  integer i, j, group, fpixel, nelements, decimals
!  integer pcount, gcount
!  logical simple, extend
!  !
!  stat=0
!  blocksize = 1
!  pcount = 0
!  gcount = 1
!  group=1
!  fpixel=1
!  decimals = 1
!  author_info_fits = 'fdu@umich.edu'
!  !
!  call ftgiou(fU, stat)
!  !
!  call ftinit(fU, filename, blocksize, stat)
!  !
!  simple=.true.
!  bitpix=-64 ! double
!  naxis=3
!  naxes(1)=im%nx
!  naxes(2)=im%ny
!  naxes(3)=im%nz
!  extend=.true.
!  !
!  call ftphpr(fU, simple, bitpix, naxis, naxes, pcount, gcount, extend, stat)
!  !
!  nelements=naxes(1)*naxes(2)*naxes(3)
!  !
!  call ftpprd(fU, group, fpixel, nelements, im%val, stat)
!  !
!  call ftpkyd(fU, 'BZERO',  0.0D0,  decimals, 'Zero point', stat)
!  call ftpkyd(fU, 'BSCALE', 1.0D0,  decimals, 'Scaling factor', stat)
!  call ftpkyd(fU, 'CDELT1', 1.0D0,  decimals, 'dx', stat)
!  call ftpkyd(fU, 'CDELT2', 1.0D0,  decimals, 'dy', stat)
!  call ftpkyd(fU, 'CDELT3', 1.0D0,  decimals, 'dz', stat)
!  call ftpkyf(fU, 'CRPIX1', 51.0,   decimals, 'i0', stat)
!  call ftpkyf(fU, 'CRPIX2', 51.0,   decimals, 'j0', stat)
!  call ftpkyf(fU, 'CRPIX3', 51.0,   decimals, 'k0', stat)
!  call ftpkyf(fU, 'CRVAL1',  0.0,   decimals, 'x0', stat)
!  call ftpkyf(fU, 'CRVAL2',  0.0,   decimals, 'y0', stat)
!  call ftpkyf(fU, 'CRVAL3',  0.0,   decimals, 'z0', stat)
!  call ftpkys(fU, 'CTYPE1', 'X', '', stat)
!  call ftpkys(fU, 'CTYPE2', 'Y', '', stat)
!  call ftpkys(fU, 'CTYPE3', 'Z', '', stat)
!  call ftpkys(fU, 'AUTHOR', author_info_fits, '', stat)
!  !
!  call ftclos(fU, stat)
!  call ftfiou(fU, stat)
!end subroutine save_fits_cube



end module disk



subroutine chem_ode_f(NEQ, t, y, ydot)
  use chemistry
  use heating_cooling
  implicit none
  integer NEQ, i, j, i1
  double precision t, y(NEQ), ydot(NEQ), rtmp, tmp
  ydot = 0D0
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    if (y(chem_species%nSpecies+1) .ne. chem_params%Tgas) then
      chem_params%Tgas = y(chem_species%nSpecies+1)
      call chem_cal_rates
    end if
  end if
  !
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 21, 64) ! A + B -> C ! 53
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(2, i))
      case (1, 2, 3, 13, 61, 20) ! A -> B
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
      case (62)
        tmp = y(chem_net%reac(1, i)) / &
          (chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain)
        if (tmp .le. 1D-9) then
          rtmp = chem_net%rates(i) * tmp
        else
          rtmp = chem_net%rates(i) * (1D0 - exp(-tmp))
        end if
      case (75)
        tmp = y(chem_net%reac(1, i)) / &
          (chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain &
           * chem_net%ABC(3, i))
        if (tmp .le. 1D-9) then
          rtmp = chem_net%rates(i) * tmp
        else
          rtmp = chem_net%rates(i) * (1D0 - exp(-tmp))
        end if
      case (63) ! gA + gA -> gB
        ! dt(N(H2)) = k_HH * <H(H-1)>
        ! Moment equation:
        ! dt(N(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * N(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * n(H) * X(gH)
        ! dt(X(H2)) = k_HH / (k_HH + k_desorb) * sigma * v * X(H) * X(gH) * n_dust / D2G
        ! Rate equation:
        ! dt(X(H2)) = k_HH * X(H)**2 / D2G
        if (chem_net%reac_names(1, i) .eq. 'gH') then
          if (chemsol_params%H2_form_use_moeq) then
            i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
            rtmp = chem_net%rates(i) * y(i1) * y(chem_net%reac(1, i))
            ydot(i1) = ydot(i1) - rtmp ! It's like H + gH -> gH2. So dt(H) -= rtmp, dt(gH) += rtmp
            ydot(chem_net%reac(1, i)) = ydot(chem_net%reac(1, i)) + rtmp
          else
            rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(1, i))
          end if
        else
          rtmp = chem_net%rates(i) * y(chem_net%reac(1, i)) * y(chem_net%reac(1, i))
        end if
      case (0)
        rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
      case default
        cycle
    end select
    !
    do j=1, chem_net%n_reac(i)
      ydot(chem_net%reac(j, i)) = ydot(chem_net%reac(j, i)) - rtmp
    end do
    do j=1, chem_net%n_prod(i)
      ydot(chem_net%prod(j, i)) = ydot(chem_net%prod(j, i)) + rtmp
    end do
  end do
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    call realtime_heating_cooling_rate(ydot(chem_species%nSpecies+1), NEQ, y)
  else
    ydot(chem_species%nSpecies+1) = 0D0
  end if
  !
end subroutine chem_ode_f




subroutine realtime_heating_cooling_rate(r, NEQ, y)
  use chemistry
  use heating_cooling
  use disk
  double precision, intent(out) :: r
  integer, intent(in) :: NEQ
  double precision, dimension(NEQ), intent(in) :: y
  hc_params%Tgas    = y(chem_species%nSpecies+1)
  hc_params%X_H2    = y(chem_idx_some_spe%i_H2)
  hc_params%X_HI    = y(chem_idx_some_spe%i_HI)
  hc_params%X_CI    = y(chem_idx_some_spe%i_CI)
  hc_params%X_Cplus = y(chem_idx_some_spe%i_Cplus)
  hc_params%X_OI    = y(chem_idx_some_spe%i_OI)
  hc_params%X_CO    = y(chem_idx_some_spe%i_CO)
  hc_params%X_H2O   = y(chem_idx_some_spe%i_H2O)
  hc_params%X_OH    = y(chem_idx_some_spe%i_OH)
  hc_params%X_E     = y(chem_idx_some_spe%i_E)
  hc_params%X_Hplus = y(chem_idx_some_spe%i_Hplus)
  hc_params%X_gH    = y(chem_idx_some_spe%i_gH)
  hc_params%R_H2_form_rate_coeff = chem_params%R_H2_form_rate_coeff
  hc_params%R_H2_form_rate = &
    get_H2_form_rate( &
      hc_params%R_H2_form_rate_coeff, &
      hc_params%X_gH, &
      hc_params%X_HI, &
      hc_params%n_gas)
  hc_Tgas = y(chem_species%nSpecies+1)
  hc_Tdust = hc_params%Tdust
  r = &
    heating_minus_cooling() * phy_SecondsPerYear / &
    (chem_params%n_gas * phy_kBoltzmann_CGS)
  !call disp_h_c_rates
end subroutine realtime_heating_cooling_rate




subroutine chem_ode_jac(NEQ, t, y, j, ian, jan, pdj)
  use chemistry
  use heating_cooling
  use trivials
  implicit none
  double precision t, rtmp, tmp, tmp1
  double precision, dimension(NEQ) :: y, pdj
  double precision, dimension(:) :: ian, jan
  integer NEQ, i, j, k, i1
  double precision dT_dt_1, dT_dt_2, del_ratio, del_0, delta_y
  double precision, dimension(NEQ) :: ydot1, ydot2
  del_ratio = 1D-3
  del_0 = 1D-12
  pdj = 0D0
  do i=1, chem_net%nReactions
    select case (chem_net%itype(i))
      case (5, 21, 64) ! A + B -> C
        if (j .EQ. chem_net%reac(1, i)) then
          if (chem_net%reac(1, i) .ne. chem_net%reac(2, i)) then
            rtmp = chem_net%rates(i) * y(chem_net%reac(2, i))
          else
            rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(2, i))
          end if
        else if (j .EQ. chem_net%reac(2, i)) then
          if (chem_net%reac(1, i) .ne. chem_net%reac(2, i)) then
            rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
          else
            rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(1, i))
          end if
        else
          rtmp = 0D0
        end if
      case (1, 2, 3, 13, 61, 20) ! A -> B
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          rtmp = chem_net%rates(i)
        end if
      case (62)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          tmp1 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain
          tmp = y(chem_net%reac(1, i)) / tmp1
          if (tmp .le. 1D-9) then
            rtmp = chem_net%rates(i) / tmp1
          else
            rtmp = chem_net%rates(i) / tmp1 * exp(-tmp)
          end if
        end if
      case (75)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          tmp1 = chem_params%ratioDust2HnucNum * chem_params%SitesPerGrain &
                 * chem_net%ABC(3, i)
          tmp = y(chem_net%reac(1, i)) / tmp1
          if (tmp .le. 1D-9) then
            rtmp = chem_net%rates(i) / tmp1
          else
            rtmp = chem_net%rates(i) / tmp1 * exp(-tmp)
          end if
        end if
      case (63) ! gA + gA -> gB
        if (chem_net%reac_names(1, i) .eq. 'gH') then
          if (chemsol_params%H2_form_use_moeq) then
            i1 = chem_species%idx_gasgrain_counterpart(chem_net%reac(1, i))
            if (j .eq. chem_net%reac(1, i)) then
              rtmp = chem_net%rates(i) * y(i1)
              pdj(i1) = pdj(i1) - rtmp
              pdj(chem_net%reac(1, i)) = pdj(chem_net%reac(1, i)) + rtmp
            else if (j .eq. i1) then
              rtmp = chem_net%rates(i) * y(chem_net%reac(1, i))
              pdj(i1) = pdj(i1) - rtmp
              pdj(chem_net%reac(1, i)) = pdj(chem_net%reac(1, i)) + rtmp
            else
              rtmp = 0D0
            end if
          else
            if (j .eq. chem_net%reac(1, i)) then
              rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(1, i))
            else
              rtmp = 0D0
            end if
          end if
        else
          if (j .eq. chem_net%reac(1, i)) then
            rtmp = 2D0 * chem_net%rates(i) * y(chem_net%reac(1, i))
          else
            rtmp = 0D0
          end if
        end if
      case (0)
        if (j .ne. chem_net%reac(1, i)) then
          rtmp = 0D0
        else
          rtmp = chem_net%rates(i)
        end if
      case default
        cycle
    end select
    !
    if (rtmp .NE. 0D0) then
      do k=1, chem_net%n_reac(i)
        pdj(chem_net%reac(k, i)) = pdj(chem_net%reac(k, i)) - rtmp
      end do
      do k=1, chem_net%n_prod(i)
        pdj(chem_net%prod(k, i)) = pdj(chem_net%prod(k, i)) + rtmp
      end do
    end if
  end do
  !
  if (chemsol_params%evolT .and. (NEQ .ge. chem_species%nSpecies+1)) then
    if (is_in_list_int(j, chem_idx_some_spe%nItem, chem_idx_some_spe%idx)) then
      call realtime_heating_cooling_rate(dT_dt_1, NEQ, y)
      delta_y = y(j) * del_ratio + del_0
      rtmp = y(j)
      y(j) = y(j) + delta_y
      call realtime_heating_cooling_rate(dT_dt_2, NEQ, y)
      pdj(chem_species%nSpecies+1) = (dT_dt_2 - dT_dt_1) / delta_y
      y(j) = rtmp
    else if (j .eq. (chem_species%nSpecies+1)) then
      call chem_ode_f(NEQ, t, y, ydot1)
      delta_y = y(j) * del_ratio + del_0
      rtmp = y(j)
      y(j) = y(j) + delta_y
      call chem_ode_f(NEQ, t, y, ydot2)
      pdj = (ydot2 - ydot1) / delta_y
      y(j) = rtmp
      chem_params%Tgas = rtmp
    end if
  else
    pdj(chem_species%nSpecies+1) = 0D0
  end if
end subroutine chem_ode_jac
