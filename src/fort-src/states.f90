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

program data
  use parser, only: parser_nfields
  integer :: ng=0
  integer :: np=0
  integer :: nfields
  real, allocatable :: mat(:,:)
  integer :: err
  integer :: ip,ig,k
  integer :: ngt0,ind
  integer, allocatable :: states(:)
  real, allocatable :: wlist(:),wsort(:)
  real :: qt0,qt1,qt2,qt3,qt4,md
  integer,parameter :: nin=5
  integer :: ns(nin+1)
  real :: qt(nin-1)
  character(1280) :: longstring
  character(100000) :: line, newline

  ! read the first string
  read(5,'(a)',iostat=err) line
  call parser_nfields(line,newline,nfields)
  rewind(5)
  ng = nfields

  np=0
  do 
     read(5,'(a)',iostat=err) line
     if(err /= 0) exit
     np = np + 1
  end do
  rewind(5)

  allocate(mat(np,ng),stat=err)
  allocate(states(ng),stat=err)

  do ip=1,np
     read(*,*) mat(ip,:)
  end do

  do ip = 1,np
     ! number of dists > 0
     ngt0 = count(mat(ip,:) > 0.0)
     ! mean dist
     md = sum(mat(ip,:),mat(ip,:) > 0.0)/ real(ngt0)
     ! deviation
     m2d = sum(mat(ip,:)**2,mat(ip,:) > 0.0)/ real(ngt0)
     m2d = sqrt(m2d - md**2)
     allocate(wlist(ngt0),stat=err)
     allocate(wsort(ngt0),stat=err)
     k = 0
     do ig = 1,ng
        if(mat(ip,ig) > 0.0) then
           k = k + 1
           wlist(k) = mat(ip,ig)
        end if
     end do
     wsort = wlist
     call sort(wsort,ngt0)

     ! look for quartiles
     do i = 1,nin-1
        ind = int((real(i)/real(nin))*ngt0)
        qt(i) = 0.5 * (wsort(ind) + wsort(ind+1))
     end do

     ! assign states     
     states = 0
     where(mat(ip,:) < 0.0) states = nin+1
     do i = 1,nin-1
        where(mat(ip,:) < qt(i) .and. states == 0) states = i
     end do
     where(states == 0) states = nin

     write(*,'(58i3)') states 

     deallocate(wlist)
     deallocate(wsort)

  end do

CONTAINS

  ! --------------------------------------------------------------------
  ! INTEGER FUNCTION  FindMinimum():
  !    This function returns the location of the minimum in the section
  ! between Start and End.
  ! --------------------------------------------------------------------

  INTEGER FUNCTION  FindMinimum(x, Start, End)
    IMPLICIT  NONE
    real, DIMENSION(1:), INTENT(IN) :: x
    INTEGER, INTENT(IN)                :: Start, End
    real                            :: Minimum
    INTEGER                            :: Location
    INTEGER                            :: i

    Minimum  = x(Start)         ! assume the first is the min
    Location = Start                    ! record its position
    DO i = Start+1, End         ! start with next elements
       IF (x(i) < Minimum) THEN !   if x(i) less than the min?
          Minimum  = x(i)               !      Yes, a new minimum found
          Location = i                !      record its position
       END IF
    END DO
    FindMinimum = Location              ! return the position
  END FUNCTION  FindMinimum

  ! --------------------------------------------------------------------
  ! SUBROUTINE  Swap():
  !    This subroutine swaps the values of its two formal arguments.
  ! --------------------------------------------------------------------

  SUBROUTINE  Swap(a, b)
    IMPLICIT  NONE
    real, INTENT(INOUT) :: a, b
    real                :: Temp

    Temp = a
    a    = b
    b    = Temp
  END SUBROUTINE  Swap

  ! --------------------------------------------------------------------
  ! SUBROUTINE  Sort():
  !    This subroutine receives an array x() and sorts it into ascending
  ! order.
  ! --------------------------------------------------------------------

  SUBROUTINE  Sort(x, Size)
    IMPLICIT  NONE
    real, DIMENSION(1:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN)                   :: Size
    INTEGER                               :: i
    INTEGER                               :: Location

    DO i = 1, Size-1                    ! except for the last
       Location = FindMinimum(x, i, Size)       ! find min from this to last
       CALL  Swap(x(i), x(Location))    ! swap this and the minimum
    END DO
  END SUBROUTINE  Sort



end program data
