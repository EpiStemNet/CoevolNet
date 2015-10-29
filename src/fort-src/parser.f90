module parser
  implicit none
  private
  public :: parser_nfields
  ! space, tab, comma, colon, semicolon
  character(1) :: delimiters(5)=[" ", achar(9), ",", ":", ";"]
  
contains

  subroutine parser_nfields(line,parsed,nfields)
    character,intent(in):: line*(*)
    integer, intent(out) :: nfields
    character,intent(out):: parsed*(*)
    integer i, n, toks
    
    i = 1
    n = len_trim(line)
    toks = 0
    nfields = 0
    parsed = ""

    outer: do while(i <= n)
       do while(any(delimiters == line(i:i)))
          i = i + 1
          if (n < i) exit outer
       enddo
       toks = toks + 1
       nfields = toks
       if(nfields == 1) then 
          parsed=trim(parsed)//line(i:i)
       else
          parsed=trim(parsed)//" "//line(i:i)
       end if
       do
          i = i + 1
          if (n < i) exit outer
          if (any(delimiters == line(i:i))) exit
          parsed=trim(parsed)//line(i:i)
       enddo
    enddo outer
    
  end subroutine parser_nfields

end module parser
