subroutine pt2_alloc(pt2_data,N)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  integer, intent(in)           :: N
  integer :: k

  ! Size of the diagnostics (debug keyword pt2_diag): the same on all the
  ! processes, as they read the same keyword, so the serialized size agrees.
  if (pt2_diag) then
    pt2_ndiag = pt2_diag_nfull
  else
    pt2_ndiag = 0
  endif

  allocate(pt2_data % pt2(N)           &
          ,pt2_data % pt2_im(N)        &
          ,pt2_data % pt2_norm(N)      &
          ,pt2_data % variance(N)      &
          ,pt2_data % rpt2(N)          &
          ,pt2_data % overlap(N,N)     &
          ,pt2_data % diag(pt2_ndiag,N) &
          )

  pt2_data % pt2(:)           = 0.d0
  pt2_data % pt2_im(:)        = 0.d0
  pt2_data % pt2_norm(:)      = 0.d0
  pt2_data % variance(:)      = 0.d0
  pt2_data % rpt2(:)          = 0.d0
  pt2_data % overlap(:,:)     = 0.d0
  pt2_data % diag(:,:)        = 0.d0

end subroutine

subroutine pt2_dealloc(pt2_data)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  deallocate(pt2_data % pt2         &
            ,pt2_data % pt2_im      &
            ,pt2_data % pt2_norm    &
            ,pt2_data % variance    &
            ,pt2_data % rpt2        &
            ,pt2_data % overlap     &
            ,pt2_data % diag        &
            )
end subroutine

subroutine pt2_add(p1, w, p2)
  implicit none
  use selection_types
  BEGIN_DOC
! p1 += w * p2
  END_DOC
  type(pt2_type), intent(inout) :: p1
  double precision, intent(in)  :: w
  type(pt2_type), intent(in)    :: p2

  if (w == 1.d0) then

    p1 % pt2(:)            = p1 % pt2(:)           + p2 % pt2(:)
    p1 % pt2_im(:)         = p1 % pt2_im(:)        + p2 % pt2_im(:)
    p1 % pt2_norm(:)       = p1 % pt2_norm(:)      + p2 % pt2_norm(:)
    p1 % rpt2(:)           = p1 % rpt2(:)          + p2 % rpt2(:)
    p1 % variance(:)       = p1 % variance(:)      + p2 % variance(:)
    p1 % overlap(:,:)      = p1 % overlap(:,:)     + p2 % overlap(:,:)
    p1 % diag(:,:)         = p1 % diag(:,:)        + p2 % diag(:,:)

  else

    p1 % pt2(:)            = p1 % pt2(:)           + w * p2 % pt2(:)
    p1 % pt2_im(:)         = p1 % pt2_im(:)        + w * p2 % pt2_im(:)
    p1 % pt2_norm(:)       = p1 % pt2_norm(:)      + w * p2 % pt2_norm(:)
    p1 % rpt2(:)           = p1 % rpt2(:)          + w * p2 % rpt2(:)
    p1 % variance(:)       = p1 % variance(:)      + w * p2 % variance(:)
    p1 % overlap(:,:)      = p1 % overlap(:,:)     + w * p2 % overlap(:,:)
    p1 % diag(:,:)         = p1 % diag(:,:)        + w * p2 % diag(:,:)

  endif

end subroutine


subroutine pt2_add2(p1, w, p2)
  implicit none
  use selection_types
  BEGIN_DOC
! p1 += w * p2**2
  END_DOC
  type(pt2_type), intent(inout) :: p1
  double precision, intent(in)  :: w
  type(pt2_type), intent(in)    :: p2

  if (w == 1.d0) then

    p1 % pt2(:)           = p1 % pt2(:)           + p2 % pt2(:)           * p2 % pt2(:)
    p1 % pt2_im(:)        = p1 % pt2_im(:)        + p2 % pt2_im(:)        * p2 % pt2_im(:)
    p1 % pt2_norm(:)      = p1 % pt2_norm(:)      + p2 % pt2_norm(:)      * p2 % pt2_norm(:)
    p1 % rpt2(:)          = p1 % rpt2(:)          + p2 % rpt2(:)          * p2 % rpt2(:)
    p1 % variance(:)      = p1 % variance(:)      + p2 % variance(:)      * p2 % variance(:)
    p1 % overlap(:,:)     = p1 % overlap(:,:)     + p2 % overlap(:,:)     * p2 % overlap(:,:)
    p1 % diag(:,:)        = p1 % diag(:,:)        + p2 % diag(:,:)        * p2 % diag(:,:)

  else

    p1 % pt2(:)           = p1 % pt2(:)           + w * p2 % pt2(:)           * p2 % pt2(:)
    p1 % pt2_im(:)        = p1 % pt2_im(:)        + w * p2 % pt2_im(:)        * p2 % pt2_im(:)
    p1 % pt2_norm(:)      = p1 % pt2_norm(:)      + w * p2 % pt2_norm(:)      * p2 % pt2_norm(:)
    p1 % rpt2(:)          = p1 % rpt2(:)          + w * p2 % rpt2(:)          * p2 % rpt2(:)
    p1 % variance(:)      = p1 % variance(:)      + w * p2 % variance(:)      * p2 % variance(:)
    p1 % overlap(:,:)     = p1 % overlap(:,:)     + w * p2 % overlap(:,:)     * p2 % overlap(:,:)
    p1 % diag(:,:)        = p1 % diag(:,:)        + w * p2 % diag(:,:)        * p2 % diag(:,:)

  endif

end subroutine


subroutine pt2_serialize(pt2_data, n, x)
  implicit none
  use selection_types
  type(pt2_type), intent(in)    :: pt2_data
  integer, intent(in)           :: n
  double precision, intent(out) :: x(*)

  integer :: i,k,n2

  n2 = n*n
  x(1:n)         =  pt2_data % pt2(1:n)
  k=n
  x(k+1:k+n)     =  pt2_data % pt2_im(1:n)
  k=k+n
  x(k+1:k+n)     =  pt2_data % pt2_norm(1:n)
  k=k+n
  x(k+1:k+n)     =  pt2_data % rpt2(1:n)
  k=k+n
  x(k+1:k+n)     =  pt2_data % variance(1:n)
  k=k+n
  x(k+1:k+n2)  =  reshape(pt2_data % overlap(1:n,1:n), (/ n2 /))
  k=k+n2
  if (pt2_ndiag > 0) then
    x(k+1:k+pt2_ndiag*n) = reshape(pt2_data % diag(1:pt2_ndiag,1:n), (/ pt2_ndiag*n /))
  endif

end

subroutine pt2_deserialize(pt2_data, n, x)
  implicit none
  use selection_types
  type(pt2_type), intent(inout) :: pt2_data
  integer, intent(in)           :: n
  double precision, intent(in)  :: x(*)

  integer :: i,k,n2

  n2 = n*n
  pt2_data % pt2(1:n)           =   x(1:n)
  k=n
  pt2_data % pt2_im(1:n)        =   x(k+1:k+n)
  k=k+n
  pt2_data % pt2_norm(1:n)      =   x(k+1:k+n)
  k=k+n
  pt2_data % rpt2(1:n)          =   x(k+1:k+n)
  k=k+n
  pt2_data % variance(1:n)      =   x(k+1:k+n)
  k=k+n
  pt2_data % overlap(1:n,1:n) = reshape(x(k+1:k+n2), (/ n, n /))
  k=k+n2
  if (pt2_ndiag > 0) then
    pt2_data % diag(1:pt2_ndiag,1:n) = reshape(x(k+1:k+pt2_ndiag*n), (/ pt2_ndiag, n /))
  endif

end
