module selection_types
  type selection_buffer
    integer :: N, cur
    integer(8)      , pointer :: det(:,:,:)
    double precision, pointer :: val(:)
    double precision          :: mini
  endtype

  type pt2_type
    double precision, allocatable :: pt2(:)
    double precision, allocatable :: pt2_im(:)
    double precision, allocatable :: pt2_norm(:)
    double precision, allocatable :: rpt2(:)
    double precision, allocatable :: variance(:)
    double precision, allocatable :: overlap(:,:)
    double precision, allocatable :: diag(:,:)   ! (pt2_ndiag, N) diagnostics of the CS PT2 (pt2_diag)
  endtype

  ! Diagnostics of the CS PT2 (keyword pt2_diag, debug only), accumulated in
  ! pt2_type % diag(:,state) by selection_cs and estimated like the other sums.
  ! pt2_ndiag is set by pt2_alloc: pt2_diag_nfull if pt2_diag, else 0 (no cost).
  ! Layout of diag(:,state), e_alpha the contribution of an external determinant:
  !   1                 sum |Re e_alpha|
  !   2                 S = sum |e_alpha|
  !   3                 sum |e_alpha|**2           (n_eff = S**2 / this)
  !   4  .. 19          4 classes x (sum Re e, sum Im e, sum |e|, count);
  !                     class 1: single, Re(Delta)<0 ; 2: single, Re(Delta)>0 ;
  !                     class 3: double, Re(Delta)<0 ; 4: double, Re(Delta)>0 ;
  !                     Delta = E - H_alpha,alpha, degree with respect to the generator
  !   20 .. 91          histogram of arg(e_alpha) in pt2_diag_nbin sectors of
  !                     ]-pi,pi], weighted by |e_alpha| : Re(Delta)<0 (20..55),
  !                     then Re(Delta)>0 (56..91)
  integer, parameter :: pt2_diag_nbin    = 36
  integer, parameter :: pt2_diag_iclass  = 4
  integer, parameter :: pt2_diag_ihist   = 20
  integer, parameter :: pt2_diag_nfull   = 19 + 2*pt2_diag_nbin
  integer            :: pt2_ndiag        = 0

  contains

  integer function pt2_type_size(N)
    implicit none
    integer, intent(in) :: N
    pt2_type_size = (5*n + n*n + pt2_ndiag*n)
  end function

end module

