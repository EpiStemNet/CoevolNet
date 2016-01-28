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

module model
  use nrtype
  use data, only: nv,ns,nd
  implicit none
  private

  public :: prm,grd

  public :: fix_gaps

  public :: model_initialize
  public :: model_set_myv
  public :: model_put_myv
  public :: model_collect_prm
  public :: compute_pseudo_likelihood

  public :: mypl
  public :: etot
  public :: fields
  public :: couplings

  ! working variable for the processor
  integer(I4B) :: myv 
  real(SP), allocatable :: fields(:,:) ! NS x NV
  real(SP), allocatable :: couplings(:,:,:) ! NS x NS x NV(NV-1)/2 
  real(SP), allocatable :: p1(:,:) ! NS x NV
  real(SP), allocatable :: p2(:,:,:) ! NS x NS x NV(NV-1)/2 

  integer(I4B), allocatable :: shuffled_data_samples(:,:)

  ! arrays for fixed residue/sequence
  real(SP), allocatable, save :: my_fields(:) ! NS
  real(SP), allocatable, save :: my_couplings(:,:,:) ! NS x NV x NS
  real(SP), allocatable, save :: my_p1(:) ! NS 
  real(SP), allocatable, save :: my_p2(:,:,:) ! NS x NV x NS
  real(SP), allocatable, save :: my_f1(:) ! NS 
  real(SP), allocatable, save :: my_f2(:,:,:) ! NS x NV x NS

  ! gradients arrays
  real(SP) :: lambda=0.01_SP ! regularization strength; the default is l2 with lambda=0.01
  integer(I4B) :: regu=2 ! the default is l2; regu=0 => no regularization 
  real(SP) :: mypl,ereg,etot
  real(SP), allocatable, save :: grd(:) ! NS + NS x NV x NS
  real(SP), allocatable, save :: prm(:) ! NS + NS x NV x NS

  logical :: fix_gaps

contains

  subroutine model_initialize(regux,lambdax)
    integer(I4B) :: err
    integer(I4B) :: iv,jv
    integer(I4B) :: ind
    integer(I4B), intent(in) :: regux
    real(SP), intent(in) :: lambdax


    regu = regux
    lambda = lambdax

!    allocate(shuffled_data_samples(nv,nd),stat=err)
    allocate(p1(ns,nv),stat=err)
    allocate(p2(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(fields(ns,nv),stat=err)
    allocate(couplings(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(my_f1(ns),stat=err)
    allocate(my_f2(ns,nv,ns),stat=err)
    allocate(my_fields(ns),stat=err)
    allocate(my_couplings(ns,nv,ns),stat=err)
    allocate(my_p1(ns),stat=err)
    allocate(my_p2(ns,nv,ns),stat=err)
    allocate(grd(ns + ns*ns*nv),stat=err)
    allocate(prm(ns + ns*ns*nv),stat=err)


    p1 = 0.0_SP
    p2 = 0.0_SP
    fields = 0.0_SP
    couplings = 0.0_SP       

  end subroutine model_initialize

  subroutine model_zero()
    real(SP), parameter :: small_number=1.e-2_SP
    my_p1 = 0.0_SP
    my_p2 = 0.0_SP
    mypl = 0.0_SP
    etot = 0.0_SP
    select case(regu)
    case(0)
       ereg = 0.0_SP
    case(1)
!       ereg = - lambda * (sum(abs(prm(:ns))) + 0.5_SP * sum(abs(prm(ns+1:))))
       ereg = - lambda * small_number * (sum(log(cosh(prm(:ns)/small_number))) + 0.5_SP * sum(log(cosh(prm(ns+1:)/small_number))))
    case(2)
       ereg = - lambda * (sum(prm(:ns)**2) + 0.5_SP * sum(prm(ns+1:)**2))
    end select
  end subroutine model_zero

  subroutine model_set_myv(iv,err) ! my_couplings
    use data, only: nd,data_samples,w
    integer(I4B), intent(in) :: iv
    ! make my_couplings given myv 
    ! must be called before looping on data
    integer(I4B) :: id,jd,jv,k,ind,mys
    integer(I4B) :: err
    integer, allocatable :: list(:)
    real(SP) :: sum,rnd

    myv = iv
    my_p1 = 0.0_SP
    my_p2 = 0.0_SP
    my_f1 = 0.0_SP
    my_f2 = 0.0_SP
    my_fields = 0.0_SP
    my_couplings = 0.0_SP       
    mypl = 0.0_SP
    etot = 0.0_SP
    ereg = 0.0_SP
    prm = 0.0_SP
    grd = 0.0_SP

    allocate(list(nv),stat=err)
    
    my_f1 = 0.0_SP
    my_f2 = 0.0_SP
    do id = 1,nd
       list = data_samples(:,id)
       mys = list(myv)
       my_f1(mys) = my_f1(mys) + w(id)
       do jv = 1,nv
          if(jv /= myv) then 
             my_f2(list(jv),jv,mys) = my_f2(list(jv),jv,mys) + w(id)
          end if
       end do
    end do

    deallocate(list)

  end subroutine model_set_myv

  subroutine model_put_myv ! my_couplings
    ! make my_couplings given myv 
    ! must be called before looping on data
    integer(I4B) :: jv,k,ind

    fields(:,myv) = my_fields
    
    !      --myv--
    !      x x x x
    !   |  1 x x x
    !  jv  2 4 x x
    !   |  3 5 6 x
    !
    ! lower triangle packing: jv > myv

    ! remove gauge before adding to couplings
    call model_gauge()

    fields(:,myv) = my_fields
    do jv = myv+1,nv ! jv > myv
       k = (myv - 1) * nv - myv * (myv + 1) / 2 + jv 
       couplings(:,:,k) = couplings(:,:,k) + my_couplings(:,jv,:)
    end do

    do jv = 1,myv-1 ! myv > jv: submatrices must be transposed 
       k = (jv - 1) * nv - jv * (jv + 1) / 2 + myv
       couplings(:,:,k) =  couplings(:,:,k) + transpose(my_couplings(:,jv,:))
    end do

  end subroutine model_put_myv

  subroutine model_gauge
    integer(I4B) :: jv,is,js
    real(SP) :: mat(ns,ns),arr(ns),marr
    real(SP) :: rsum(ns),csum(ns),totsum

    arr = my_fields
    marr = sum(arr) / real(ns)
    arr = arr - marr
    do jv = 1,nv
       mat = my_couplings(:,jv,:)
       totsum = sum(mat)
       totsum = totsum / real(ns*ns)
       do is = 1,ns
          rsum(is) = sum(mat(is,:))
          csum(is) = sum(mat(:,is))
       end do
       rsum = rsum / real(ns)
       csum = csum / real(ns)
       if(jv /= myv) arr = arr + csum - totsum
       do js = 1,ns
          do is = 1,ns
             mat(is,js) = mat(is,js) - rsum(is) - csum(js) + totsum
          end do
       end do
       my_couplings(:,jv,:) = mat
    end do
    my_fields = arr
    
  end subroutine model_gauge

  subroutine compute_pseudo_conp()
    use data, only: data_samples,w,nd
    integer(I4B) :: list(nv)
    real(SP) :: conp(ns)
    integer(I4B) :: is,jv,kv
    real(SP) :: r,rsum
    real :: start,finish
    real(SP) :: pp,pp0,zz,zz0
    real(SP) :: tmp(nv)
    integer(I4B):: mys
    integer :: ind
    integer :: id
    integer(I4B) :: dd(nv)
    real(SP) :: ww,rnd
    integer(I4B) :: rd

    ! loop over data
    do id = 1,nd
       list = data_samples(:,id)
       ww = w(id)
       mys = list(myv)
       
       ! loop over the states of myv 
       do is = 1,ns
          r = my_fields(is) 
          do jv = 1,nv
             if(myv /= jv) then 
                r = r + my_couplings(list(jv),jv,is) 
             end if
          end do
          conp(is) = exp(r)
       end do
       
       rsum = sum(conp)
       conp = conp / rsum
       
       mypl = mypl + ww * log(conp(mys))
       
       ! update histograms 
       ! loop over the states of myv 
       do is = 1,ns
          pp = conp(is) * ww
          my_p1(is) = my_p1(is) + pp 
          do jv = 1,nv
             if(myv /= jv) then 
                my_p2(list(jv),jv,is) = my_p2(list(jv),jv,is) + pp
             end if
          end do
       end do
       
    end do
    
  end subroutine compute_pseudo_conp

  subroutine compute_pseudo_likelihood(it)
    integer(I4B), intent(in) :: it
    real(SP) :: etot0,de

    call model_parameters_unpack()

    etot0 = etot
    call model_zero()
    call compute_pseudo_conp()
    etot = mypl + ereg
    de = etot-etot0

    call model_parameters_pack()
    
  end subroutine compute_pseudo_likelihood

  subroutine model_parameters_pack
    integer(I4B) :: dim
    real(SP), parameter :: small_number=1.e-2_SP
    integer(I4B) :: is,jv,js

    dim = ns*ns*nv
    ! compute the gradient
    my_p1 = my_p1 - my_f1 
    my_p2 = my_p2 - my_f2 

    select case(regu)
    case(0)
       ! do nothing; no regularization
    case(1) 
       ! l1 regularization
       my_p1 = my_p1 + lambda * tanh(my_fields/small_number)
       my_p2 = my_p2 + 0.5_SP * lambda * tanh(my_couplings/small_number)
    case(2) 
       ! l2 regularization
       my_p1 = my_p1 + 2.0_SP * lambda * my_fields
       my_p2 = my_p2 + 2.0_SP * 0.5_SP * lambda * my_couplings
    end select

    prm(1:ns) = my_fields
    prm(ns+1:) = reshape(my_couplings,(/dim/))

    grd(1:ns) = my_p1
    grd(ns+1:) = reshape(my_p2,(/dim/))
    
  end subroutine model_parameters_pack

  subroutine model_parameters_unpack
    integer :: k

    my_fields = prm(1:ns) 
    my_couplings = reshape(prm(ns+1:),(/ns,nv,ns/))
    
  end subroutine model_parameters_unpack

  subroutine model_collect_prm

    couplings = 0.5_SP * couplings

  end subroutine model_collect_prm
  

end module model
  
