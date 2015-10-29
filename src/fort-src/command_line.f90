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
