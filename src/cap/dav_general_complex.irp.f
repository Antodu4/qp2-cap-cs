
subroutine davidson_general_complex(u_in,H_jj,energies,dim_in,sze,N_st,N_st_diag_in,converged,h_mat)
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! u_in : guess coefficients on the various states. Overwritten on exit
  !
  ! dim_in : leftmost dimension of u_in
  !
  ! sze : Number of determinants
  !
  ! N_st : Number of eigenstates
  !
  ! N_st_diag_in : Number of states in which H is diagonalized. Assumed > sze
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)             :: dim_in, sze, N_st, N_st_diag_in
  complex*16,  intent(in)   :: H_jj(sze),h_mat(sze,sze)
  complex*16, intent(inout) :: u_in(dim_in,N_st_diag_in)
  complex*16, intent(out)   :: energies(N_st)

  integer                        :: iter, N_st_diag
  integer                        :: i,j,k,l,m
  logical, intent(inout)         :: converged

  double precision, external     :: u_dot_v, u_dot_u

  integer                        :: k_pairs, kl

  integer                        :: iter2, itertot
  complex*16, allocatable        :: y(:,:), h(:,:), lambda(:), tmp(:,:)
  complex*16, allocatable        :: residual_norm(:), h_cp(:,:)
  character*(16384)              :: write_buffer
  complex*16                     :: to_print(2,N_st), res
  double precision               :: cpu, wall
  integer                        :: shift, shift2, itermax, istate
  double precision               :: r1, r2, alpha, val
  integer                        :: nproc_target, info
  integer                        :: order(N_st_diag_in)
  double precision               :: cmax
  complex*16      , allocatable  :: U(:,:), overlap(:,:)!, S_d(:,:)
  complex*16      , pointer      :: W(:,:)
  logical                        :: disk_based
  double precision, allocatable  :: max_U(:)
  integer, allocatable           :: pos_U(:)
  complex*16                     :: energy_shift(N_st_diag_in*davidson_sze_max)

  include 'constants.include.F'

  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, y, h, lambda
  if (N_st_diag*3 > sze) then
    print *,  'error in Davidson :'
    print *,  'Increase n_det_max_full to ', N_st_diag*3
    stop -1
  endif

  itermax = max(2,min(davidson_sze_max, sze/N_st_diag))+1
  itertot = 0

  if (state_following) then
    allocate(overlap(N_st_diag*itermax, N_st_diag*itermax))
  else
    allocate(overlap(1,1))  ! avoid 'if' for deallocate
  endif
  overlap = 0.d0

  provide threshold_davidson !nthreads_davidson
  call write_time(6)
  write(6,'(A)') ''
  write(6,'(A)') 'Davidson Diagonalization'
  write(6,'(A)') '------------------------'
  write(6,'(A)') ''

  
  ! Find max number of cores to fit in memory
  ! -----------------------------------------

  nproc_target = nproc
  double precision :: rss
  integer :: maxab
  maxab = sze 

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 2d0 * 8.d0 *                                   &!complex  bytes
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.d0*dble(sze*m)*(N_st_diag*itermax)  &! W
         + 2.0d0*(N_st_diag*itermax)**2           &! h,y
         + 2.d0*(N_st_diag*itermax)               &! s2,lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
                                                   ! In H_S2_u_0_nstates_zmq
         + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on collector
         + 3.d0*(N_st_diag*N_det)                 &! u_t, v_t, s_t on slave
         + 0.5d0*maxab                            &! idx0 in H_S2_u_0_nstates_openmp_work_*
         + nproc_target *                         &! In OMP section
           ( 1.d0*(N_int*maxab)                   &! buffer
           + 3.5d0*(maxab) )                      &! singles_a, singles_b, doubles, idx
         ) / 1024.d0**3

    if (nproc_target == 0) then
      call check_mem(r1,irp_here)
      nproc_target = 1
      exit
    endif

    if (r1+rss < qp_max_mem) then
      exit
    endif

    if (itermax > 4) then
      itermax = itermax - 1
    else if (m==1.and.disk_based_davidson) then
      m=0
      disk_based = .True.
      itermax = 6
    else
      nproc_target = nproc_target - 1
    endif

  enddo
  nthreads_davidson = nproc_target
  TOUCH nthreads_davidson
  call write_int(6,N_st,'Number of states')
  call write_int(6,N_st_diag,'Number of states in diagonalization')
  call write_int(6,sze,'Number of basis functions')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif

  !---------------

  write(6,'(A)') ''
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================  ===========  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = 'Iter'
  do i=1,N_st
    if (i==1) then
    write_buffer = trim(write_buffer)//'       Energy            Im            Residual       Im     '
    else 
    write_buffer = trim(write_buffer)//'            Energy            Im            Residual       Im     '
    endif
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================  ===========  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)


!  if (disk_based) then
!    ! Create memory-mapped files for W and S
!    type(c_ptr) :: ptr_w, ptr_s
!    integer :: fd_s, fd_w
!    call mmap(trim(ezfio_work_dir)//'davidson_w', (/int(sze,8),int(N_st_diag*itermax,8)/),&
!        8, fd_w, .False., ptr_w)
!    call mmap(trim(ezfio_work_dir)//'davidson_s', (/int(sze,8),int(N_st_diag*itermax,8)/),&
!        4, fd_s, .False., ptr_s)
!    call c_f_pointer(ptr_w, w, (/sze,N_st_diag*itermax/))
!    call c_f_pointer(ptr_s, s, (/sze,N_st_diag*itermax/))
!  else
    allocate(W(sze,N_st_diag*itermax))
!  endif

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax),                                      &
      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      h_cp(N_st_diag*itermax,N_st_diag*itermax),                     &
      tmp(N_st_diag*itermax,N_st_diag*itermax),                      &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      residual_norm(N_st_diag),                                      &
      lambda(N_st_diag*itermax),                                     &
      max_U(N_st_diag*itermax),pos_U(N_st_diag*itermax))

  h = (0.d0, 0d0)
  U = (0.d0, 0d0)
  y = (0.d0, 0d0)

  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)

  ! Davidson iterations
  ! ===================

  converged = .False.

  ! Initialize from N_st to N_st_diat with gaussian random numbers
  ! to be sure to have overlap with any eigenvectors
  do k=N_st+1,N_st_diag
    do i=1,sze
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      u_in(i,k) = dcmplx(r1*dcos(r2), 0d0) * u_in(i,k-N_st)
    enddo
    u_in(k,k) = u_in(k,k) + (10.d0,0d0)
  enddo
  ! Normalize all states 
  complex*16 :: norm

  ! Copy from the guess input "u_in" to the working vectors "U"
  U=dcmplx(1d300,1d300)
  do k=1,N_st_diag
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
  enddo

  do while (.not.converged)
    itertot = itertot+1
    if (itertot == 8) then
      exit
    endif

    do iter=1,itermax-1

      shift  = N_st_diag*(iter-1)
      shift2 = N_st_diag*iter
   
      ! Compute |W_k> = \sum_i |i><i|H|u_k>
      ! -----------------------------------

      ! Gram-Schmidt to orthogonalize all new guess with the previous vectors
      call qr_decomposition_c(U,size(U,1),sze,shift2)

      ! Change the sign of the guess vectors to match with the ones of the previous
      ! iterations
      if (iter > 1) then
        do j = 1, shift
          if (sign(1d0,max_U(j)) * sign(1d0,dble(U(pos_U(j),j))) < 0d0) then
            do i = 1, sze
              U(i,j) = - U(i,j)
            enddo
            !print*,'sign change!'
          endif
        enddo
      endif

      ! Sign of each vector
      do j = shift+1, shift2
        val = 0d0
        do i = 1, sze
          if (dabs(dble(U(i,j))) > dabs(val)) then
            val = dble(U(i,j))
            max_U(j) = dble(U(i,j))
            pos_U(j) = i
          endif
        enddo
      enddo

      call hpsi_complex(W(:,shift+1),U(:,shift+1),N_st_diag,sze,h_mat)

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call zgemm('T','N', shift2, shift2, sze,                       &
          (1.d0,0d0), U, size(U,1), W, size(W,1),                          &
          (0.d0,0d0), h, size(h,1))

      ! Diagonalize h y = lambda y
      ! ---------------

      h_cp = h
      call diag_general_complex(lambda,y,h_cp,size(h,1),shift2,info)
      call qr_decomposition_c(y,size(y,1),shift2,shift2)

      ! Compute Energy for each eigenvector
      ! -----------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), h, size(h,1), y, size(y,1),                          &
          (0d0,0.d0), tmp, size(tmp,1))

      call zgemm('T','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), tmp, size(tmp,1),                  &
          (0.d0,0d0), h, size(h,1))

      do k=1,shift2
        lambda(k) = h(k,k)
      enddo

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------

      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), U, size(U,1), y, size(y,1), (0d0,0.d0), U(1,shift2+1), size(U,1))

      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), W, size(W,1), y, size(y,1), (0.d0,0d0), W(1,shift2+1), size(W,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      do k=1,N_st_diag
        do i=1,sze
          if (dble(H_jj(i) - lambda (k)) >= 1d-2) then
            U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k) ) /(H_jj(i) - lambda (k))
          else
            U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
             / dcmplx(1d-2, dimag(H_jj(i) - lambda (k)))
          endif
        enddo

        if (k <= N_st) then
          call inner_product_complex(U(1,shift2+k),U(1,shift2+k),sze, residual_norm(k))
          to_print(1,k) = lambda(k) + nuclear_repulsion
          to_print(2,k) = residual_norm(k)
        endif
      enddo

      if ((itertot>1).and.(iter == 1)) then
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,F16.10,1X,ES12.3,1X,ES12.3,1X))') iter-1, to_print(1:2,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
          converged = dabs(maxval(dble(residual_norm(1:N_st)))) < threshold_davidson
      endif   
      

      do k=1,N_st
        if (dble(residual_norm(k)) > 1.e8) then
          print *, 'Davidson failed'
          stop -1
        endif
      enddo
      if (converged) then
        exit
      endif

      logical, external :: qp_stop
      if (qp_stop()) then
        converged = .True.
        exit
      endif

    enddo

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        W, size(W,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        U, size(U,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo

    do j = 1, N_st_diag
      val = 0d0
      do i = 1, sze
        if (dabs(dble(U(i,j))) > dabs(val)) then
          val = dble(U(i,j))
          max_U(j) = dble(U(i,j))
          pos_U(j) = i
        endif
      enddo
    enddo

    do j = 1, N_st_diag
      if (sign(1d0,max_U(j)) * sign(1d0,dble(U(pos_U(j),j))) < 0d0) then
        do i = 1, sze
          u_in(i,j) = - U(i,j)
        enddo
      else
        do i = 1, sze
          u_in(i,j) = U(i,j)
        enddo
      endif
    enddo
  enddo

  do k=1,N_st
    energies(k) = lambda(k)
  enddo
  write_buffer = '====='
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================  ================  ===========  ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''

  call write_time(6)
  write(6,'(A)') ''

    deallocate(W)

  deallocate (                                                       &
      residual_norm,                                                 &
      U, h,                                                 &
      y,                                                             &
      lambda                                                         &
      )
  deallocate(overlap)
  FREE nthreads_davidson
end

subroutine hpsi_complex(v,u,N_st,sze,h_mat)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v = H_mat * u for all N_st states simultaneously.
  END_DOC
  integer, intent(in)              :: N_st,sze
  complex*16, intent(in)     :: u(sze,N_st),h_mat(sze,sze)
  complex*16, intent(inout)  :: v(sze,N_st)
  integer :: i,j,istate
  v = (0.d0, 0d0)
  do istate = 1, N_st
   do i = 1, sze
    do j = 1, sze
      v(i,istate) += h_mat(j,i) * u(j,istate)
    enddo
   enddo
  enddo
end

