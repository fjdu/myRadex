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


program main

use configure
use my_timer
use trivials

implicit none

integer i, j

type(atimer) timer
!type(date_time) a_date_time

call get_command_argument(0, filename_exe,    i, j)
call get_command_argument(1, filename_config, i, j)
if (i .EQ. 0) then
  filename_config = 'configure.dat'
end if

! Load the configure file
call config_do

call timer%init('')

! Do the work
call do_my_radex

call timer%elapse

end program main
