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
