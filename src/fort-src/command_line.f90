! Copyright (c) 2015 Simone Marsili
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module command_line
  implicit none
  private
  public :: command_line_read
  character(2) :: opts(4) = ['-i','-w','-l','-g']

contains

  subroutine command_line_read(nerrs,data_file,wid,regu,lambda,fix_gaps)
    use nrtype
    integer, intent(out) :: nerrs
    character(1000), intent(out) :: data_file
    real(SP), intent(out) :: wid
    integer(I4B), intent(out) :: regu
    real(SP), intent(out) :: lambda
    character(1000) :: cmd
    integer :: nargs
    character(100) :: arg
    integer :: iarg
    logical :: fix_gaps

    call get_command(cmd)
    nargs = command_argument_count()
    write(0,*) trim(cmd)

    iarg = 1
    nerrs = 0
    wid = 0.0_SP
    regu = 0
    fix_gaps = .false.
    do while(iarg <= nargs)
       call get_command_argument(iarg,arg)
       if (any(opts == trim(arg))) then
          select case(trim(arg))
          case('-i','--input')
             ! input file
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             data_file = arg
             if(len_trim(data_file) == 0) then
                nerrs = nerrs + 1
                write(0,*) 'error: no data file name'
             end if
             if(data_file(1:1) == '-') then
                nerrs = nerrs + 1
                write(0,*) 'error: no data file name'
             end if
          case('-w','--weigths')
             ! weights file
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             read(arg,*) wid
          case('-l')
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             regu = 2
             read(arg,*) lambda
          case('-g')
             ! gaps are independent variables
             fix_gaps = .true.
          end select
       else
          write(0,'(a,1x,a)') 'error: unknown option',trim(arg)
          nerrs = nerrs + 1
       end if
       iarg = iarg + 1
    end do

  end subroutine command_line_read



end module command_line
