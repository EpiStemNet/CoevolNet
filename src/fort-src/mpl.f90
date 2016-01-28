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

program mpl
  use nrtype
  use units
  use command_line, only: &
       command_line_read
  use data, only: &
       nv, &
       data_read
  use model, only: &
       model_initialize, &
       model_set_myv, &
       model_collect_prm, &
       fix_gaps
  use scrs, only: &
       compute_scores, &
       print_scores
  use dvmlm_wrapper, only: &
       dvmlm_minimize
  implicit none
  integer(I4B) :: err,toterr
  integer(I4B) :: iv
  integer(I4B) :: niter,neval
  real :: finish,start,start_min,end_min
  integer(I4B) :: udata,uscrs
  character(1000) :: data_file,scores_file
  real(SP) :: w_id
  integer(I4B) :: regu
  real(SP) :: lambda

  call units_initialize()

  toterr=0
  ! get command line args 
  call command_line_read(err,data_file,w_id,regu,lambda,fix_gaps)
  if(len_trim(data_file) == 0) then 
     toterr = toterr + 1
     write(0,*) 'error: data filename is needed'
  end if

  ! open data file
  call units_open(data_file,udata,'O',err)
  if(err /= 0) then 
     toterr = toterr + 1
     write(0,*) 'error: cant find data file'
  end if
  scores_file = trim(data_file)//'.scores'
  call units_open(scores_file,uscrs,'U',err)

  if (toterr > 0) then 
     write(0,*) 'program will exit'
     stop
  end if

  write(0,*) 'reading data...'
  call data_read(udata,w_id)

  write(0,*) 'initiale model...'
  call model_initialize(regu,lambda)

  call cpu_time(start_min)

  ! compute averages for the model 
  write(0,*) 'minimization...'

  do iv = 1,nv
     call model_set_myv(iv,err)
     niter = 0
     call cpu_time(start)
     call dvmlm_minimize(niter,neval)
     call cpu_time(finish)
     write(0,'(a,i5,a,2i5,a,f8.3,a)') ' variable ', iv, &
          '  converged (niter,neval) ', niter, neval, ' in ', finish-start, ' secs'
  end do
  
  flush(0)

  call cpu_time(end_min)
  write(0,*) 'minimization total (secs): ', end_min-start_min
  flush(0)
  call model_collect_prm()
  call compute_scores()
  call print_scores(uscrs)

end program mpl
