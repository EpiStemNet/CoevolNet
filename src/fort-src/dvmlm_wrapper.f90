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

module dvmlm_wrapper
  use nrtype
  use model, only: prm, grd, etot
  ! wrapper to dvmlm subroutine 
  implicit none
  private 

  public :: dvmlm_minimize

contains

  subroutine dvmlm_min(ndim,mstep,task,wa)
    integer(I4B), intent(in) :: ndim,mstep
    character(60), intent(inout) :: task
    real(SP), intent(in) :: wa(:)
    integer(I4B) :: isave(5)
    real(SP) :: dsave(24)
    real(SP) :: f
    real(SP) :: frtol
    real(SP) :: fatol
    real(SP) :: fmin
    external dvmlm

    ! set prms for minimization
    ! frtol: desired relative error
    ! fatol: desired absolute error
    !    frtol = 1.0e-16_SP
    !    fatol = 1.0e-20_SP
    !    frtol = 1.0e-8_SP
    !    fatol = 1.0e-10_SP
    frtol = 1.0e-4_SP
    fatol = 1.0e-5_SP
    fmin = -1.9e30_SP
    
    f = - etot
    call dvmlm(ndim,prm,f,grd,frtol,fatol,fmin,task,mstep,&
         wa(1),wa(ndim*mstep+1),wa(2*ndim*mstep+1),&
         isave,dsave,wa(2*ndim*mstep+mstep+1),wa(2*ndim*mstep+mstep+ndim+1))

  end subroutine dvmlm_min

  subroutine dvmlm_minimize(iter,totiter)
    use model, only: compute_pseudo_likelihood, model_put_myv
    integer(I4B),intent(out) :: iter,totiter
    integer(I4B) :: err
    integer(I4B) :: ndim,mstep
    character(60) :: task
    real(SP), allocatable :: wa(:)
    integer :: lwa
    
    iter = 0
    totiter = 0
    task = 'START'

    ndim = size(prm)
    ! this is the number of steps for hessian approximation
    mstep = 100
    lwa = 2*ndim*mstep + 2*ndim + mstep
    allocate(wa(lwa),stat=err)
    

    call compute_pseudo_likelihood(iter)
    do 
       if(totiter > 100) then 
          write(0,*) 'warning: totiter > 100'
          flush(0)
       end if
       call dvmlm_min(ndim,mstep,task,wa)
       if(task(1:2) == 'FG') then 
          ! update etot and gradient for line search
          totiter = totiter + 1
          call compute_pseudo_likelihood(iter)
       elseif(task(1:4) == 'NEWX') then
          ! start new line search
          iter = iter + 1
       elseif(task(1:4) == 'WARN') then 
          write(0,*) 'warning ', iter
          flush(0)
       elseif(task(1:4) == 'CONV') then 
          ! compute final values for likelihood 
          call compute_pseudo_likelihood(iter)
          ! put my prms back in fields and couplings arrays 
          call model_put_myv()
          exit
       end if
    end do

    deallocate(wa)

  end subroutine dvmlm_minimize

end module dvmlm_wrapper
