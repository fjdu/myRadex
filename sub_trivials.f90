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


module trivials

implicit none

contains

!  Open a file for sequential read.
subroutine openFileSequentialRead (fU, filename, maxRowLen, getu)
integer, intent(inout) :: fU
character(len=*), intent(in) :: filename
integer, intent(in), optional :: maxRowLen, getu
integer maxRecLen, ios
if (present(getu)) then
  if (getu .ne. 0) then
    if (.not. getFileUnit(fU)) then
      write(*,*) 'Cannot get a file unit.'
      return
    end if
  end if
end if
if (present(maxRowLen)) then
  maxRecLen = maxRowLen
else
  maxRecLen = 9999
end if
open (UNIT=fU, FILE=trim(filename), IOSTAT=ios, &
     STATUS='OLD', ACCESS='SEQUENTIAL', FORM='FORMATTED', &
     RECL=maxRecLen, BLANK='NULL', POSITION='REWIND', &
     ACTION='READ', DELIM='NONE', PAD='YES')
if (ios .NE. 0) then
  write (*, '(A)') 'In openFileSequentialRead:'
  write (*, '(/A, I8, /A, A)') 'Open File Error: IOSTAT=', ios, &
    'Filename: ', filename
  ! Some people say that a subroutine should not
  ! terminate the whole program. But, ...
  stop
end if
end subroutine openFileSequentialRead


!  Open a file for sequential write.
subroutine openFileSequentialWrite (fU, filename, maxRowLen, getu)
integer, intent(inout) :: fU
character(len=*), intent(in) :: filename
integer, intent(in), optional :: maxRowLen, getu
integer maxRecLen, ios
if (present(getu)) then
  if (getu .ne. 0) then
    if (.not. getFileUnit(fU)) then
      write(*,*) 'Cannot get a file unit.'
      return
    end if
  end if
end if
if (present(maxRowLen)) then
  maxRecLen = maxRowLen
else
  maxRecLen = 9999
end if
open (UNIT=fU, FILE=trim(filename), IOSTAT=ios, &
     STATUS='REPLACE', ACCESS='SEQUENTIAL', FORM='FORMATTED', &
     RECL=maxRecLen, BLANK='NULL', POSITION='REWIND', &
     ACTION='WRITE', DELIM='NONE', PAD='YES')
if (ios .NE. 0) then
  write (*, '(A)') 'In openFileSequentialWrite'
  write (*, '(A, I8, /A, A)') 'Open File Error: IOSTAT=', ios, &
    'Filename: ', filename
  stop
end if
end subroutine openFileSequentialWrite


!  Get a file unit to output.
function getFileUnit (fU)
integer fU
logical getFileUnit
logical uExist, uOpened
! 0: standard error
! 5: standard input
! 6: standard output
! 101:
! 102:
do fU=10, 99
  inquire (unit=fU, exist=uExist, opened=uOpened)
  if (uExist .AND. .NOT. (uOpened)) then
    getFileUnit = .TRUE.
    return
  end if
end do
getFileUnit = .FALSE.
fU = -1
return
end function getFileUnit


!  Check if a file unit is opened.
function FileUnitOpened (fU)
logical FileUnitOpened
integer fU
inquire (unit=fU, opened=FileUnitOpened)
return
end function FileUnitOpened


function is_in_list_int(n, d, list)
  logical is_in_list_int
  integer, intent(in) :: n, d
  integer, dimension(d), intent(in) :: list
  integer i
  if (d .lt. 1) then
    is_in_list_int = .false.
    return
  end if
  do i=1, d
    if (n .eq. list(i)) then
      is_in_list_int = .true.
      return
    end if
  end do
  is_in_list_int = .false.
  return
end function is_in_list_int


subroutine my_mkdir(dir)
  character(len=*) dir
  call system('mkdir ' // trim(dir))
end subroutine my_mkdir


subroutine my_cp_to_dir(filename, dir)
  character(len=*) filename, dir
  call system('cp -p ' // trim(filename) // ' ' // trim(dir))
end subroutine my_cp_to_dir


function file_exist(filename)
  logical file_exist
  character(len=*) filename
  INQUIRE(FILE=filename, EXIST=file_exist)
end function file_exist


function dir_exist(dirname)
  logical dir_exist
  character(len=*) dirname
  character(len=128) filename_tmp
  integer ios, fU
  write(filename_tmp, '("rXq", I16, "AjYnWd")') get_simple_rand_integer()
  !filename_tmp = trim(dirname) // trim(filename_tmp)
  filename_tmp = trim(combine_dir_filename(trim(adjustl(dirname)), trim(adjustl(filename_tmp))))
  if (file_exist(filename_tmp)) then
    dir_exist = .true.
    return
  else
    if (.NOT. getFileUnit(fU)) then
      write(*,*) 'Cannot get a file unit in dir_exist!'
      stop
    end if
    open(unit=fU, file=filename_tmp, iostat=ios, action='WRITE')
    if (ios .NE. 0) then
      dir_exist = .false.
    else
      dir_exist = .true.
      close(fU, status='DELETE')
      !call system('rm -f ' // trim(filename_tmp))
    end if
  end if
end function dir_exist


function combine_dir_filename(dir, filename)
  character(len=*) dir, filename
  character(len=256) combine_dir_filename
  integer i, i1
  logical found
  found = .false.
  do i=len_trim(dir), 1, -1
    if (dir(i:i) .eq. '/') then
      if (i .eq. 1) then
        found = .true.
        i1 = i
        exit
      else if (dir((i-1):(i-1)) .ne. '/') then
        found = .true.
        i1 = i
        exit
      end if
    end if
  end do
  if (.not. found) then
    combine_dir_filename = trim(dir) // '/' // filename
  else
    combine_dir_filename = trim(dir(1:i1)) // filename
  end if
end function combine_dir_filename


subroutine dropout_char(str, ch)
  character(len=*), intent(inout) :: str
  character, intent(in) :: ch
  character, dimension(:), allocatable :: stmp
  integer i, i0, n
  n = len(str)
  allocate(stmp(n))
  i0 = 0
  do i=1, n
    if (str(i:i) .ne. ch) then
      i0 = i0 + 1
      stmp(i0) = str(i:i)
    end if
  end do
  do i=1, i0
    str(i:i) = stmp(i)
  end do
  do i=i0+1, n
    str(i:i) = ''
  end do
  deallocate(stmp)
end subroutine dropout_char


function integer2char(n)
  integer, intent(in) :: n
  character(len=8) integer2char
  write(integer2char, '(I8)') n
end function integer2char


function get_simple_rand_integer()
  integer get_simple_rand_integer, i
  associate(j => get_simple_rand_integer)
    call system_clock(count=j)
    i = int(abs(sin(real(j))*1E6))
    i = iand(i, 65535) * 36969 + ishft(i, -16)
    j = iand(j, 65535) * 18000 + ishft(j, -16)
    get_simple_rand_integer = ishft(i, 16) + j
  end associate
end function get_simple_rand_integer


!  Get the number of lines in a file.
subroutine GetFileLen_ (fU, FileName, nFileLen)
integer fU, nFileLen, ios
character(len=*) FileName
character strtmp
CALL openFileSequentialRead &
  (fU, FileName, 999)
nFileLen = 0
do
  read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
  if (ios .LT. 0) exit
  nFileLen = nFileLen + 1
end do
close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
end subroutine GetFileLen_


!  Get the number of lines in a file.
function GetFileLen(FileName)
  integer GetFileLen
  integer fU, ios
  character(len=*) FileName
  character strtmp
  GetFileLen = 0
  if (.NOT. getFileUnit(fU)) then
    write(*,*) 'No freee file unit!'
    return
  end if
  CALL openFileSequentialRead(fU, FileName, 999)
  do
    read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
    if (ios .LT. 0) exit
    GetFileLen = GetFileLen + 1
  end do
  close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
end function GetFileLen


!  Get the number of lines in a file.
function GetFileLen_comment(FileName, commentchar)
integer GetFileLen_comment
integer fU, ios
character(len=*) FileName
character commentchar
character(len=32) strtmp
GetFileLen_comment = 0
if (.NOT. getFileUnit(fU)) then
  write(*,*) 'No freee file unit!'
  return
end if
CALL openFileSequentialRead(fU, FileName, 999)
do
  read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
  if (ios .LT. 0) exit
  strtmp = adjustl(strtmp)
  if (strtmp(1:1) .NE. commentchar) then
    GetFileLen_comment = GetFileLen_comment + 1
  end if
end do
close(UNIT=fU, IOSTAT=ios, STATUS='KEEP')
end function GetFileLen_comment


!  Get the number of lines in a file.
function GetFileLen_comment_blank(FileName, commentchar)
integer GetFileLen_comment_blank
integer fU, ios
character(len=*) FileName
character commentchar
character(len=32) strtmp
GetFileLen_comment_blank = 0
if (.NOT. getFileUnit(fU)) then
  write(*,*) 'No freee file unit!'
  return
end if
CALL openFileSequentialRead(fU, FileName, 999)
do
  read (UNIT=fU, FMT='(A)', IOSTAT=ios) strtmp
  if (ios .LT. 0) exit
  strtmp = adjustl(strtmp)
  if ((strtmp(1:1) .NE. commentchar) .AND. &
      (len_trim(strtmp) .NE. 0)) then
      !(strtmp(1:1) .NE. ' ')) then
    GetFileLen_comment_blank = GetFileLen_comment_blank + 1
  end if
end do
close(UNIT=fU, IOSTAT=ios, STATUS='KEEP')
end function GetFileLen_comment_blank


function getFilePreName(strFileName)
integer ntrim, i
character(len=*) strFileName
character(len=128) getFilePreName
ntrim = len_trim(strFileName)
do i=ntrim, 1, -1
  if (strFileName(i:i) .EQ. '.') exit
end do
if (i .EQ. 1) then
  getFilePreName = strFileName(1:ntrim)
else
  getFilePreName = strFileName(1:(i-1))
end if
end function getFilePreName


function str_pad_to_len(str, len)
  character(len=*) str
  integer l, len
  character(len=len) :: str_pad_to_len
  l = len_trim(str)
  str_pad_to_len = ''
  str_pad_to_len(len-l+1 : len) = trim(str)
end function str_pad_to_len


function IsWordChar(ch)
logical IsWordChar
character ch
IsWordChar = &
  (LGE(ch, '0') .AND. LLE(ch, '9')) .OR. &
  (LGE(ch, 'A') .AND. LLE(ch, 'Z')) .OR. &
  (LGE(ch, 'a') .AND. LLE(ch, 'z')) .OR. &
  (ch .EQ. '_')
end function IsWordChar


function IsDigitChar(ch)
logical IsDigitChar
character ch
IsDigitChar = &
  (LGE(ch, '0') .AND. LLE(ch, '9'))
end function IsDigitChar


function CharCountsInStr(Str, Ch)
integer CharCountsInStr
character(len=*) Str
character Ch
integer :: nTmp, i
CharCountsInStr = 0
nTmp = len(Str)
do i=1, nTmp
  if (Str(i:i) .EQ. Ch) &
    CharCountsInStr = CharCountsInStr + 1
end do
end function CharCountsInStr


subroutine write_no_advance(str)
  character(len=*) str
  character(len=*), parameter :: pre = char(27)//'[A'
  write(*, '(A)') pre//str
end subroutine write_no_advance



function binary_search(list, ndim, key, step)
! Assuming list is sorted in either ascending or declining order
  integer binary_search
  integer, intent(in) :: ndim
  double precision, dimension(ndim), intent(in) :: list
  double precision, intent(in) :: key
  integer, intent(in) :: step
  integer i, n, imid, itmp, itmp1, itmp2, imin, imax
  logical asc
  n = ndim/step
  imin = 1
  imax = n
  if (list(ndim) .GT. list(1)) then
    asc = .TRUE.
  else
    asc = .FALSE.
  end if
  do i=1, n
    imid = (imin + imax) / 2
    itmp = (imid - 1) * step + 1
    if (key .LT. list(itmp)) then
      if (asc) then
        imax = imid - 1
      else
        imin = imid + 1
      end if
    else if (key .GT. list(itmp)) then
      if (asc) then
        imin = imid + 1
      else
        imax = imid - 1
      end if
    else
      binary_search = itmp
      return
    end if
    if (imin .GT. imax) then
      itmp1 = max(itmp - step, 1)
      itmp2 = min(itmp + step, ndim)
      if (abs(list(itmp1) - key) .LT. abs(list(itmp) - key)) then
        binary_search = itmp1
      else if (abs(list(itmp2) - key) .LT. abs(list(itmp) - key)) then
        binary_search = itmp2
      else
        binary_search = itmp
      end if
      return
    else if (imax .LT. 1) then
      binary_search = 0
      return
    else if (imin .GT. ndim) then
      binary_search = ndim+1
      return
    end if
  end do
end function binary_search



function sign_dble(x)
  integer sign_dble
  double precision x
  if (x .GT. 0D0) then
    sign_dble = 1
  else if (x .LT. 0D0) then
    sign_dble = -1
  else
    sign_dble = 0
  end if
end function sign_dble


function calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
  double precision :: calc_four_point_linear_interpol
  double precision, intent(in) :: x, y, x1, x2, y1, y2, z11, z12, z21, z22
  double precision k1, k2
  k1 = (z12-z11) / (y2-y1)
  k2 = (z22-z21) / (y2-y1)
  associate( &
    z0_1 => z11, &
    z0_2 => z21)
    associate( &
      k_k  => (k2-k1)/(x2-x1), &
      k_0  => k1, &
      k_z0 => (z0_2-z0_1)/(x2-x1), &
      dx => x-x1, &
      dy => y-y1)
      calc_four_point_linear_interpol = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
    end associate
  end associate
end function calc_four_point_linear_interpol


subroutine load_array_from_txt(filename, array, ncol, nrow, nx, ny, commentstr)
  character(len=*), intent(in) :: filename
  double precision, dimension(:,:), allocatable, intent(out) :: array
  character(len=64) fmtstr
  character(len=1024) tmpstr
  character commentchar
  integer i, fU, ios, nc, nr, nx_, ny_
  integer, intent(out), optional :: ncol, nrow, nx, ny
  character(len=*), intent(out), optional :: commentstr
  if (.not. getFileUnit(fU)) then
    write(*,*) 'Cannot get a free file unit.  In load_array_from_txt.'
    stop
  end if
  call openFileSequentialRead(fU, filename, 99999)
  ! Search for the format part
  do
    call read_a_nonempty_row(fU, fmtstr, '(A128)', ios)
    if (ios .ne. 0) then
      exit
    end if
    ios = index(fmtstr, 'Format:')
    if (ios .ne. 0) then
      exit
    end if
  end do
  !
  read(fU, '(A1, 4I10)') commentchar, nc, nr, nx_, ny_
  if (present(ncol)) then
    ncol = nc
  end if
  if (present(nrow)) then
    nrow = nr
  end if
  if (present(nx)) then
    nx = nx_
  end if
  if (present(ny)) then
    ny = ny_
  end if
  if (.not. allocated(array)) then
    allocate(array(nc, nr))
  end if
  read(fU, '(1X, A64)') fmtstr
  i = 0
  do
    read(fU, '(A)', IOSTAT=ios) tmpstr
    if (ios .lt. 0) then
      exit
    end if
    if ((tmpstr(1:1) .eq. commentchar) .or. (len_trim(tmpstr) .eq. 0)) then
      if (present(commentstr)) then
        commentstr = trim(commentstr) // trim(tmpstr)
      end if
      cycle
    end if
    i = i + 1
    read(tmpstr, fmtstr) array(:, i)
  end do
  close(fU)
end subroutine load_array_from_txt


subroutine split_str_by_space(str, str_split, n, nout)
  character(len=*), intent(in) :: str
  integer, intent(in) :: n
  character(len=*), dimension(n), intent(out) :: str_split
  integer, intent(out), optional :: nout
  integer nlen, i, istart, iend, nfound
  logical flag
  character TABCHAR
  TABCHAR = char(9)
  !
  nfound = 0
  !
  nlen = len_trim(str)
  ! Get the first non-blank position
  istart = 1
  do
    ! Find the starting position of the current substring
    flag = .true.
    do i=istart, nlen
      if ((len_trim(str(i:i)) .ne. 0) .and. (str(i:i) .ne. TABCHAR)) then
        flag = .false.
        exit
      end if
    end do
    if (flag) then
      exit
    end if
    ! Find the end position of the current substring
    istart = i
    iend = i + 1
    flag = .true.
    do i=iend, nlen
      if ((len_trim(str(i:i)) .eq. 0) .or. (str(i:i) .eq. TABCHAR)) then
        flag = .false.
        exit
      end if
    end do
    if (flag) then
      iend = nlen
    else
      iend = i - 1
    end if
    !
    nfound = nfound + 1
    if (nfound .gt. n) then
      nfound = n
      exit
    end if
    str_split(nfound) = str(istart : iend)
    istart = iend + 1
    if (istart .gt. nlen) then
      exit
    end if
  end do
  if (present(nout)) then
    nout = nfound
  end if 
end subroutine split_str_by_space


subroutine join_str_arr(sarr, n, lenS, sSep, sout)
  integer, intent(in) :: n, lenS
  character(len=*), intent(in) :: sSep
  character(len=*), dimension(n), intent(in) :: sarr
  character(len=lenS), intent(out) :: sout
  integer i
  sout = trim(adjustl(sarr(1)))
  do i=2,n-1
    sout = trim(sout) // sSep // trim(adjustl(sarr(i)))
  end do
  if (n .gt. 1) then
    sout = trim(sout) // sSep // trim(adjustl(sarr(n)))
  end if
end subroutine join_str_arr


subroutine read_a_nonempty_row(fU, str, fmtstr, ios)
  integer, intent(in) :: fU
  character(len=*), intent(out) :: str
  character(len=*), intent(in) :: fmtstr
  integer, intent(out) :: ios
  do
    read(fU, fmt=fmtstr, iostat=ios) str
    if ((ios .ne. 0) .or. (len_trim(str) .gt. 0)) then
      exit
    end if
  end do
end subroutine read_a_nonempty_row



function get_num_of_interval_log(xmin, xmax, dx0, del_ratio)
  integer get_num_of_interval_log
  double precision, intent(in) :: xmin, xmax, dx0, del_ratio
  get_num_of_interval_log = ceiling(log( &
         (xmax-xmin)/dx0 * (del_ratio - 1D0) + 1D0) / log(del_ratio))
end function get_num_of_interval_log


function get_ratio_of_interval_log(xmin, xmax, dx0, n, iout)
  double precision get_ratio_of_interval_log
  double precision, intent(in) :: xmin, xmax, dx0
  integer, intent(in) :: n
  integer, intent(out), optional :: iout
  integer :: i, imax = 100
  double precision, parameter :: frac_precision = 1D-3
  double precision k, p
  double precision tmp
  k = (xmax - xmin) / dx0
  p = 1D0 / dble(n)
  get_ratio_of_interval_log = &
    exp(p*log(k + 1D0)) - 1D0
  do i=1, imax
    tmp = get_ratio_of_interval_log
    get_ratio_of_interval_log = &
      exp(p * log(get_ratio_of_interval_log*k + 1D0)) - 1D0
    if (abs(tmp - get_ratio_of_interval_log) .le. get_ratio_of_interval_log * frac_precision) then
      exit
    end if
  end do
  get_ratio_of_interval_log = get_ratio_of_interval_log + 1D0
  if (present(iout)) then
    iout = i
  end if
end function get_ratio_of_interval_log


function discrete_integral(n, x, y, xmin, xmax) result(s)
  double precision s
  integer, intent(in) :: n
  double precision, dimension(:), intent(in) :: x, y
  double precision, intent(in) :: xmin, xmax
  double precision k, b, x1, x2
  integer i, i1, i2
  !
  if ((x(1) .ge. xmax) .or. (x(n) .le. xmin)) then
    s = 0D0
    return
  end if
  i1 = 0
  i2 = 0
  do i=1, n-1
    if ((x(i) .le. xmin) .and. (x(i+1) .ge. xmin)) then
      i1 = i
    end if
    if ((x(i) .le. xmax) .and. (x(i+1) .ge. xmax)) then
      i2 = i
    end if
  end do
  x1 = xmin
  x2 = xmax
  if (i1 .eq. 0) then
    i1 = 1
    x1 = x(1)
  end if
  if (i2 .eq. 0) then
    i2 = n-1
    x2 = x(n)
  end if
  s = 0D0
  do i=i1, i2
    k = (y(i+1) - y(i)) / (x(i+1) - x(i))
    b = y(i) - k * x(i)
    if (i1 .eq. i2) then
      s = (x2 - x1) * (0.5D0*k*(x2+x1) + b)
      return
    end if
    if (i .eq. i1) then
      s = s + (x(i+1) - x1) * (0.5D0*k*(x(i+1)+x1) + b)
    else if (i .eq. i2) then
      s = s + (x2 - x(i)) * (0.5D0*k*(x2+x(i)) + b)
    else
      s = s + (x(i+1) - x(i)) * (0.5D0*k*(x(i+1)+x(i)) + b)
    end if
  end do
end function discrete_integral



function tau2beta(tau, factor)
  double precision tau2beta
  double precision, intent(in) :: tau
  double precision, intent(in), optional :: factor
  double precision fac
  double precision, parameter :: const_small_num = 1D-8
  if (.not. present(factor)) then
    fac = 3D0
  else
    fac = factor
  end if
  if (tau .le. const_small_num) then
    tau2beta = 1D0
  else
    tau2beta = (1D0 - exp(-fac * tau)) / (fac * tau)
  end if
end function tau2beta


subroutine display_string_both(str, fU)
  character(len=*), intent(in) :: str
  integer, intent(in) :: fU
  write(*, '(A)') trim(str)
  if (FileUnitOpened(fU)) then
    write(fU, '(A)') trim(str)
  end if
end subroutine display_string_both


Pure Function to_upper (str) Result (string)
! From: https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string
    Integer :: ic, i
    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do
End Function to_upper


end module trivials

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module
  public :: QsortC
  private :: Partition

  contains

  recursive subroutine QsortC(A)
    double precision, intent(in out), dimension(:) :: A
    integer :: iq

    if(size(A) > 1) then
       call Partition(A, iq)
       call QsortC(A(:iq-1))
       call QsortC(A(iq:))
    endif
  end subroutine QsortC

  subroutine Partition(A, marker)
    double precision, intent(in out), dimension(:) :: A
    integer, intent(out) :: marker
    integer :: i, j
    double precision :: temp
    double precision :: x      ! pivot point
    x = A(1)
    i= 0
    j= size(A) + 1

    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine Partition

end module qsort_c_module



module my_timer
  private
  !
  type, public :: atimer
    character(len=8) :: name_
    real start_cpu_time, current_cpu_time, elapsed_cpu_time
    contains
      procedure :: init => init_atimer
      procedure :: elapse => elapse_atimer
      procedure :: elapsed_time => elapsed_time_atimer
  end type atimer
  !
  type, public :: date_time
    character(len=16) date_str, time_str
    contains
      procedure :: show => print_date_time_str
      procedure :: date_time_str => get_date_time_str
  end type date_time
  !
  contains
  !
    subroutine init_atimer(this, name_)
      class(atimer) this
      character(len=*), intent(in), optional :: name_
      call cpu_time(this%start_cpu_time)
      if (present(name_)) then
        this%name_ = name_
      else
        this%name_ = ''
      end if
    end subroutine init_atimer
    !
    function elapsed_time_atimer(this)
      class(atimer) this
      real elapsed_time_atimer
      call cpu_time(this%current_cpu_time)
      this%elapsed_cpu_time = this%current_cpu_time - this%start_cpu_time
      elapsed_time_atimer = this%elapsed_cpu_time
    end function elapsed_time_atimer
    !
    subroutine elapse_atimer(this)
      class(atimer) this
      call cpu_time(this%current_cpu_time)
      this%elapsed_cpu_time = this%current_cpu_time - this%start_cpu_time
      write(*,  '(/"Timer ", A, " Seconds elapsed: ", F10.3/)') &
        this%name_, this%elapsed_cpu_time
    end subroutine elapse_atimer
    !
    subroutine print_date_time_str(this)
      class(date_time) this
      call date_and_time(date=this%date_str, time=this%time_str)
      write(*, '("Current date and time: ", A4, "-", A2, "-", A2, 2X, A2, ":", A2, ":", A6)') &
        this%date_str(1:4), this%date_str(5:6), this%date_str(7:8), &
        this%time_str(1:2), this%time_str(3:4), this%time_str(5:10)
    end subroutine print_date_time_str
    !
    function get_date_time_str(this)
      class(date_time) this
      character(len=32) get_date_time_str
      call date_and_time(date=this%date_str, time=this%time_str)
      write(get_date_time_str, '(A4, "-", A2, "-", A2, 2X, A2, ":", A2, ":", A6)') &
        this%date_str(1:4), this%date_str(5:6), this%date_str(7:8), &
        this%time_str(1:2), this%time_str(3:4), this%time_str(5:10)
    end function get_date_time_str
end module my_timer
