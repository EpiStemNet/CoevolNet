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
