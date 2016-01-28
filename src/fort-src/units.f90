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

module units
  implicit none
  private 
  public :: units_initialize
  public :: units_open
  
  integer :: nunits

contains

  subroutine units_initialize
    implicit none
    
   nunits = 10
    
  end subroutine units_initialize
  
  subroutine units_open(filename,fileunit,flag,err)
    implicit none
    character(1000), intent(in) :: filename
    character(1), intent(in) :: flag
    integer, intent(out) :: fileunit
    integer, intent(out) :: err

    
    nunits = nunits + 1
    fileunit = nunits
    
    select case(flag)
       case('O')
          open(unit=fileunit,file=filename,status='OLD',iostat=err)
       case('N')
          open(unit=fileunit,file=filename,status='NEW',iostat=err)
       case('R')
          open(unit=fileunit,file=filename,status='REPLACE',iostat=err)
       case('U')
          open(unit=fileunit,file=filename,status='UNKNOWN',iostat=err)
       case default 
          open(unit=fileunit,file=filename,status='UNKNOWN',iostat=err)
    end select

    if(err /= 0) then 
       write(0,*) "error opening file ", filename, err,fileunit
    end if
    
  end subroutine units_open
  
end module units
