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


module configure

use trivials
use my_radex

implicit none

character(len=128) :: filename_exe, filename_config = ''


contains


subroutine config_do
  use my_timer
  integer fU
  !
  if (.NOT. getFileUnit(fU)) then
    write(*,*) 'Cannot get a file unit!'
    stop
  end if
  !
  ! First set to a nonharmful value.
  rdxx_cfg%ndens    = 0D0
  rdxx_cfg%n_H2     = 0D0
  rdxx_cfg%n_HI     = 0D0
  rdxx_cfg%n_oH2    = 0D0
  rdxx_cfg%n_pH2    = 0D0
  rdxx_cfg%n_Hplus  = 0D0
  rdxx_cfg%n_E      = 0D0
  rdxx_cfg%n_He     = 0D0
  !
  call openFileSequentialRead(fU, filename_config, 99999)
  !
  read(fU, nml=rdxx_configure)
  !
  close(fU, status='KEEP')
  !
  !! Make a backup of the configure file.
  !if (file_exist(trim(combine_dir_filename(a_book_keeping%dir, a_book_keeping%filename_log)))) then
  !  write(*,*) trim(a_disk_iter_params%iter_files_dir), ' is not empty!'
  !  write(*,*) 'I would rather not overwrite it.'
  !  stop
  !else
  !  call my_cp_to_dir(filename_config, a_book_keeping%dir)
  !end if
  !
  if (.NOT. dir_exist(rdxx_cfg%dir_save)) then
    call my_mkdir(rdxx_cfg%dir_save)
  end if
  !
  !if (.NOT. getFileUnit(a_book_keeping%fU)) then
  !  write(*,*) 'Cannot get a file unit for logging.'
  !end if
  !call openFileSequentialWrite(a_book_keeping%fU, &
  !  trim(combine_dir_filename(a_book_keeping%dir, a_book_keeping%filename_log)), 9999)
  !write(a_book_keeping%fU, '(A)') '! Current time: ' // trim(a_date_time%date_time_str())
  !flush(a_book_keeping%fU)
  !!
  !if (a_disk%backup_src) then
  !  write(*,*) 'Backing up your source code...'
  !  call my_cp_to_dir(a_disk%filename_exe, a_book_keeping%dir)
  !  call system(trim(a_disk%backup_src_cmd) // ' ' // trim(a_book_keeping%dir))
  !  write(*,*) 'Source code backup finished.'
  !end if
end subroutine config_do

end module configure
