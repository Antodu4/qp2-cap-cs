subroutine davidson_diag_hs2_complex_cs(dets_in,u_in,s2_out,dim_in,energies,sze,N_st,N_st_diag,Nint,converged)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization for Complex Scaling (CS)
  !
  ! Diagonalizes H_CS = e^{-i*theta} * (H + (e^{-i*theta} - 1) * T)
  !
  ! dets_in   : bitmasks corresponding to determinants
  !
  ! u_in      : guess coefficients on the various states. Overwritten on exit
  !             with the converged eigenvectors (normalized with the Hermitian
  !             product; diagonalize_ci_cs c-normalizes them afterwards)
  !
  ! dim_in    : leftmost dimension of u_in
  !
  ! sze       : number of determinants
  !
  ! N_st      : number of eigenstates to compute
  !
  ! N_st_diag : number of states used in the diagonalization (>= N_st)
  !
  ! energies  : (output) complex eigenvalues of H_CS for each state
  !
  ! s2_out    : (output) expectation value of S^2 for each state
  !
  ! converged : (output) .true. if the Davidson procedure converged
  !
  ! Initial guess vectors are not necessarily orthonormal
  END_DOC
  integer, intent(in)            :: dim_in, sze, N_st, N_st_diag, Nint
  integer(bit_kind), intent(in)  :: dets_in(Nint,2,sze)
  complex*16, intent(inout) :: u_in(dim_in,N_st_diag)
  complex*16, intent(out)  :: energies(N_st_diag), s2_out(N_st_diag)
  logical, intent(out)           :: converged
  complex*16, allocatable  :: H_jj(:)

  double precision, external     :: diag_H_mat_elem, diag_S_mat_elem
  complex*16, external           :: diag_H_mat_elem_cs
  integer                        :: i,k,l
  complex*16                     :: e_itheta_local
  ASSERT (N_st > 0)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  PROVIDE mo_two_e_integrals_in_map theta_cs
  allocate(H_jj(sze))

  ! H^{CS}_{jj} = e^{-i*theta} * H^{CAP/CS}_{jj}
  ! = e^{-i*theta} * (H_{jj} + (e^{-i*theta}-1)*T_{jj})
  e_itheta_local = dcmplx(dcos(theta_cs), -dsin(theta_cs))

  H_jj(1) = e_itheta_local * (dcmplx(diag_h_mat_elem(dets_in(1,1,1),Nint), 0d0) &
           + diag_H_mat_elem_cs(dets_in(1,1,1),Nint))
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(sze,H_jj, dets_in,Nint,e_itheta_local)          &
      !$OMP  PRIVATE(i)
  !$OMP DO SCHEDULE(static)
  do i=2,sze
    H_jj(i) = e_itheta_local * (dcmplx(diag_H_mat_elem(dets_in(1,1,i),Nint), 0d0) &
             + diag_H_mat_elem_cs(dets_in(1,1,i),Nint))
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call davidson_diag_hjj_sjj_complex(dets_in,u_in,H_jj,S2_out,energies,dim_in,sze,N_st,N_st_diag,Nint,converged)
  deallocate (H_jj)
end


subroutine davidson_diag_hjj_sjj_complex(dets_in,u_in,H_jj,s2_out,energies,dim_in,sze,N_st,N_st_diag_in,Nint,converged)
  use bitmasks
  use mmap_module
  implicit none
  BEGIN_DOC
  ! Davidson diagonalization with specific diagonal elements of the H matrix
  !
  ! H_jj : specific diagonal H matrix elements to diagonalize de Davidson
  !
  ! S2_out : Output : s^2
  !
  ! dets_in : bitmasks corresponding to determinants
  !
  ! u_in : guess coefficients on the various states. Overwritten
  !   on exit
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
  integer, intent(in)           :: dim_in, sze, N_st, N_st_diag_in, Nint
  integer(bit_kind), intent(in) :: dets_in(Nint,2,sze)
  complex*16,  intent(in)       :: H_jj(sze)
  complex*16,  intent(inout)    :: s2_out(N_st_diag_in)
  complex*16, intent(inout)     :: u_in(dim_in,N_st_diag_in)
  complex*16, intent(out)       :: energies(N_st_diag_in)

  integer                       :: iter, N_st_diag
  integer                       :: i,j,k,l,m
  logical, intent(inout)        :: converged

  integer                       :: k_pairs, kl

  integer                       :: iter2, itertot
  complex*16, allocatable       :: y(:,:), h(:,:), h_p(:,:), lambda(:), s2(:), h_cp(:,:), su(:,:)
  complex*8, allocatable        :: y_s(:,:)
  complex*16, allocatable       :: s_(:,:), s_tmp(:,:), y_left(:,:), s_cp(:,:), prev_y(:,:)
  double precision              :: diag_h_mat_elem
  complex*16, allocatable       :: residual_norm(:)
  character*(16384)             :: write_buffer
  double precision              :: to_print(5,N_st)
  double precision              :: cpu, wall
  integer                       :: shift, shift2, itermax, istate
  double precision              :: r1, r2, alpha, t1, t2
  logical                       :: state_ok(N_st_diag_in*davidson_sze_max)
  integer                       :: nproc_target, info
  integer                       :: order(N_st_diag_in)
  double precision              :: cmax, val
  complex*16, allocatable       :: U(:,:), overlap(:,:), S_d(:,:)
  double precision, allocatable :: max_U(:)
  integer, allocatable          :: pos_U(:)
  complex*16, pointer           :: W(:,:)
  complex*8, pointer            :: S(:,:)
  logical                       :: disk_based
  complex*16                    :: energy_shift(N_st_diag_in*davidson_sze_max), res

  ! Debug of the state following (cs_debug_state_following)
  complex*16, allocatable       :: dbg_u0(:,:), dbg_uold(:,:)
  double precision, allocatable :: dbg_z(:), dbg_hn(:), dbg_oh(:,:), dbg_o0(:,:)
  logical, allocatable          :: dbg_ok(:)
  complex*16                    :: dbg_c
  double precision, external    :: dbg_sf_ovh

  ! State following and restarts (cs_sf_fixed_ref, cs_sf_window, cs_dav_*)
  complex*16, allocatable       :: u_ref(:,:), u_thick(:,:), y_keep(:,:), y_keep_prev(:,:)
  complex*16                    :: e_win_ref(N_st), e_tot
  logical                       :: window_on, win_ok, thick_next
  integer                       :: off, idx0, n_rej_s2, n_rej_im, n_rej_win, n_thick
  character*(16)                :: dv_why

  include 'constants.include.F'

  N_st_diag = N_st_diag_in
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: U, W, S, y, y_s, S_d, h, lambda
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
  overlap = (0.d0,0d0)

  PROVIDE nuclear_repulsion expected_s2 psi_bilinear_matrix_order psi_bilinear_matrix_order_reverse threshold_davidson_pt2 threshold_davidson_from_pt2 theta_cs

  ! CS phase factor e^{-i*theta} for converting electronic energy to total energy
  complex*16 :: e_itheta_cs
  e_itheta_cs = dcmplx(dcos(theta_cs), -dsin(theta_cs))

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
  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1

  m=1
  disk_based = .False.
  call resident_memory(rss)
  do
    r1 = 2d0 * 8.d0 *                                   &! complex bytes
         ( dble(sze)*(N_st_diag*itermax)          &! U
         + 1.5d0*dble(sze*m)*(N_st_diag*itermax)  &! W,S
         + 1.d0*dble(sze)*(N_st_diag)             &! S_d
         + 4.5d0*(N_st_diag*itermax)**2           &! h,y,y_s,s_,s_tmp
         + 2.d0*(N_st_diag*itermax)               &! s2,lambda
         + 1.d0*(N_st_diag)                       &! residual_norm
         + 3d0*(N_st_diag*itermax)**2             &! h_cp,y_left,s_cp
         + 2d0*(N_st_diag*itermax)                &! max_U,pos_U
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
  call write_int(6,sze,'Number of determinants')
  call write_int(6,nproc_target,'Number of threads for diagonalization')
  call write_double(6, r1, 'Memory(Gb)')
  if (disk_based) then
    print *, 'Using swap space to reduce RAM'
  endif


  write(6,'(A)') ''
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================ ========== ========== =========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = ' Iter'
  do i=1,N_st
    if (i==1) then
    write_buffer = trim(write_buffer)//'       Energy                           S^2                Residual             '
    else
    write_buffer = trim(write_buffer)//'       Energy                           S^2                Residual             '
    endif
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================ ========== ========== =========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*N_st*2)


  allocate(W(sze,N_st_diag*itermax), S(sze,N_st_diag*itermax))

  allocate(                                                          &
      ! Large
      U(sze,N_st_diag*itermax), S_d(sze,N_st_diag),                  &

      ! Small
      h(N_st_diag*itermax,N_st_diag*itermax),                        &
      y(N_st_diag*itermax,N_st_diag*itermax),                        &
      prev_y(N_st_diag*itermax,N_st_diag*itermax),                   &
      s_(N_st_diag*itermax,N_st_diag*itermax),                       &
      s_tmp(N_st_diag*itermax,N_st_diag*itermax),                    &
      residual_norm(N_st_diag),                                      &
      s2(N_st_diag*itermax),                                         &
      y_s(N_st_diag*itermax,N_st_diag*itermax),                      &
      lambda(N_st_diag*itermax),                                     &
      h_cp(N_st_diag*itermax,N_st_diag*itermax),                     &
      y_left(N_st_diag*itermax,N_st_diag*itermax),                   &
      s_cp(N_st_diag*itermax,N_st_diag*itermax),                     &
      max_U(N_st_diag*itermax), pos_U(N_st_diag*itermax))

  h = (0.d0,0d0)
  U = (0.d0,0d0)
  y = (0.d0,0d0)
  s_ = (0.d0,0d0)
  s_tmp = (0.d0,0d0)
  pos_U = 0
  max_U = 0d0

  prev_y = (0d0,0d0)
  do i = 1, N_st_diag*itermax
    prev_y(i,i) = (1d0,0d0)
  enddo

  ASSERT (N_st > 0)
  ASSERT (N_st_diag >= N_st)
  ASSERT (sze > 0)
  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)

  ! Davidson iterations
  ! ===================

  converged = .False.

  if (cs_debug_state_following) then
    ! u_0 : input vectors, kept frozen to measure the drift of the reference u_in
    allocate(dbg_u0(sze,N_st), dbg_uold(sze,N_st))
    allocate(dbg_z(N_st_diag*itermax), dbg_hn(N_st_diag*itermax), dbg_ok(N_st_diag*itermax))
    allocate(dbg_oh(N_st_diag*itermax,N_st), dbg_o0(N_st_diag*itermax,N_st))
    do l = 1, N_st
      dbg_u0(1:sze,l) = u_in(1:sze,l)
    enddo
    call dbg_sf_legend
  endif

  ! u_ref : input vectors of the call, frozen (reference of cs_sf_fixed_ref,
  ! drift measured by cs_dav_debug)
  if (cs_sf_fixed_ref .or. cs_dav_debug) then
    allocate(u_ref(sze,N_st))
    do l = 1, N_st
      u_ref(1:sze,l) = u_in(1:sze,l)
    enddo
  endif
  if (cs_dav_thick_restart) then
    allocate(u_thick(sze,N_st_diag))
    allocate(y_keep(N_st_diag*itermax,N_st_diag), y_keep_prev(N_st_diag*itermax,N_st_diag))
  endif
  window_on = (cs_sf_window > 0.d0) .and. state_following
  e_win_ref = (0.d0,0.d0)
  off       = 0
  n_rej_s2  = 0
  n_rej_im  = 0
  n_rej_win = 0
  n_thick   = 0

  ! Guess
  do k=N_st+1,N_st_diag
    if (cs_dav_retry_ritz) then
      ! Retry of diagonalize_ci_cs: the extra vectors are the Ritz vectors of
      ! the previous call, keep them
      if (maxval(cdabs(u_in(1:sze,k))) > 0.d0) cycle
    endif
    do i=1,sze
        call random_number(r1)
        call random_number(r2)
        r1 = dsqrt(-2.d0*dlog(r1))
        r2 = dtwo_pi*r2
        u_in(i,k) = dcmplx(r1*dcos(r2),0d0) * u_in(i,k-N_st)
    enddo
    u_in(k,k) = u_in(k,k) + (10.d0,0d0)
  enddo
  
  do k=1,N_st_diag
    do i=1,sze
      U(i,k) = u_in(i,k)
    enddo
  enddo

  ! Orthogonalization
  if (do_qr_dav) then
    ! Sign of each vector
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,val)
    do j = 1, N_st_diag
      val = 0d0
      do i = 1, sze
        if (dabs(dble(U(i,j))) > dabs(val)) then
          val = dble(U(i,j))
          max_U(j) = val
          pos_U(j) = i
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO

    call ortho_qr_complex(U,size(U,1),sze,N_st_diag) 
    !call qr_decomposition_c(U,size(U,1),sze,N_st_diag)

    ! Change the sign of the guess vectors to match with the ones of the guess
    ! vectors
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
    do j = 1, N_st_diag
      if (sign(1d0,max_U(j)) * sign(1d0,dble(U(pos_U(j),j))) < 0d0) then
        do i = 1, sze
          U(i,j) = - U(i,j)
        enddo
      endif
    enddo
    !$OMP END PARALLEL DO
  endif

  do while (.not.converged)
    itertot = itertot+1
    if (itertot > max(1,cs_dav_max_cycles)) then
      exit
    endif

    if (cs_dav_thick_restart) then
      y_keep      = (0.d0,0.d0)
      y_keep_prev = (0.d0,0.d0)
    endif

    ! off = N_st_diag after a thick restart: the cycle starts with the current
    ! and the previous Ritz vectors (2 N_st_diag vectors). Not possible if
    ! itermax was reduced too much (zggev failure below).
    if (itermax-1 - off/N_st_diag < 1) off = 0
    iter = 0
    do while (iter < itermax-1 - off/N_st_diag)
      iter += 1

      if (iter == 1) then
        shift  = 0
      else
        shift  = N_st_diag*(iter-1) + off
      endif
      shift2 = N_st_diag*iter + off

        if (do_qr_dav) then
          ! Orthogonalization of the guess vectors
          call ortho_qr_complex(U,size(U,1),sze,shift2)
          !call qr_decomposition_c(U,size(U,1),sze,shift2)

          ! Change the sign of the guess vectors to match with the ones of the previous 
          ! iterations
          if (iter > 1) then
            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
            do j = 1, shift
              if (sign(1d0,max_U(j)) * sign(1d0,dble(U(pos_U(j),j))) < 0d0) then
                do i = 1, sze
                  U(i,j) = - U(i,j)
                enddo
              endif
            enddo
            !$OMP END PARALLEL DO
          endif

          ! Sign of each vector
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,val)
          do j = shift+1, shift2
            val = 0d0
            do i = 1, sze
              if (dabs(dble(U(i,j))) > dabs(val)) then
                val = dble(U(i,j))
                max_U(j) = val
                pos_U(j) = i
              endif
            enddo
          enddo
          !$OMP END PARALLEL DO
        endif

        ! N_st_diag new vectors, or 2 N_st_diag at the first iteration after a thick restart
        do j = shift+1, shift2, N_st_diag
          call H_S2_u_0_nstates_openmp_complex_cs(W(1,j),S_d,U(1,j),N_st_diag,sze)
          S(1:sze,j:j+N_st_diag-1) = cmplx(real(dble(S_d(1:sze,1:N_st_diag))), real(dimag(S_d(1:sze,1:N_st_diag))))
        enddo

      ! Compute s_kl = <u_k | S_l> = <u_k| S2 |u_l>
      ! -------------------------------------------
      ! The subspace is projected with the Hermitian product, as in the CAP
      ! Davidson: U is orthonormalized by a Hermitian QR, so U^H U ~ Id stays
      ! well conditioned, whereas U^T U can become nearly singular for genuinely
      ! complex vectors (spurious eigenvalues, zggev info=1, divergence).
      ! diagonalize_ci_cs c-normalizes the converged vectors and recomputes the
      ! energy as the c-product quotient U^T H U / U^T U afterwards.

       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k) COLLAPSE(2)
       do j=1,shift2
         do i=1,shift2
           s_(i,j) = 0.d0
           do k=1,sze
             s_(i,j) = s_(i,j) + DCONJG(U(k,i)) * dcmplx(dble(S(k,j)),dble(aimag(S(k,j))))
           enddo
          enddo
        enddo
        !$OMP END PARALLEL DO

      ! Compute h_kl = <u_k | W_l> = <u_k| H |u_l>
      ! -------------------------------------------

      call zgemm('C','N', shift2, shift2, sze,                       &
          (1.d0,0d0), U, size(U,1), W, size(W,1),                          &
          (0.d0,0d0), h, size(h,1))

      call zgemm('C','N', shift2, shift2, sze,                       &
          (1.d0,0d0), U, size(U,1), U, size(U,1),                          &
          (0.d0,0d0), s_tmp, size(s_tmp,1))

      ! cs_sf_window: reference energy of each state, the Hermitian Rayleigh
      ! quotient of its input vector at the first iteration of the call
      if (window_on .and. (itertot == 1) .and. (iter == 1)) then
        do l = 1, N_st
          e_win_ref(l) = h(l,l) / s_tmp(l,l) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0)
          if (cs_dav_debug) then
            write(6,'(A,I3,A,2F18.10,A,ES10.2)') ' DV window state', l, '  E_ref =', e_win_ref(l), &
                '  half-width =', cs_sf_window
          endif
        enddo
      endif

      ! Diagonalize the projected Hamiltonian using the generalized eigenvalue problem
      ! h y = lambda s_tmp y
      h_cp = h
      s_cp = s_tmp
      call lapack_zggev(h_cp,s_tmp,size(h,1),shift2,lambda,y_left,y,info)

       ! Normalization to have Y^* s_tmp Y = Id
       ! => y = 1/sqrt(y^* s_tmp y)
       ! --------------------------------------

       call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), s_cp, size(h,1), y, size(y,1),                          &
          (0d0,0.d0), s_tmp, size(s_tmp,1))

       call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), s_tmp, size(s_tmp,1),                  &
          (0.d0,0d0), s_cp, size(h,1))

       do i = 1, shift2
         y(:,i) = y(:,i) / cdsqrt(s_cp(i,i))
       enddo

       if (info > 0) then
         ! Numerical errors propagate. We need to reduce the number of iterations
         print*,'info:',info
         itermax = iter-1
         shift2 = shift2 - N_st_diag
         y = prev_y
         exit
       endif

      ! Compute Energy for each eigenvector
      ! -----------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), h, size(h,1), y, size(y,1),                          &
          (0d0,0.d0), s_tmp, size(s_tmp,1))

      call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), s_tmp, size(s_tmp,1),                  &
          (0.d0,0d0), h, size(h,1))

      do k=1,shift2
        lambda(k) = h(k,k)
      enddo

      ! Compute S2 for each eigenvector
      ! -------------------------------

      call zgemm('N','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), s_, size(s_,1), y, size(y,1),                        &
          (0.d0,0d0), s_tmp, size(s_tmp,1))

      call zgemm('C','N',shift2,shift2,shift2,                       &
          (1.d0,0d0), y, size(y,1), s_tmp, size(s_tmp,1),                  &
          (0.d0,0d0), s_, size(s_,1))

      do k=1,shift2
        s2(k) = s_(k,k)
      enddo

      if (only_expected_s2) then
          do k=1,shift2
            state_ok(k) = ((dabs(dble(s2(k))-expected_s2) < 0.6d0) .and. (dabs(dimag(s2(k))) < 0.1d0))
          enddo
      else
        do k=1,size(state_ok)
          state_ok(k) = .True.
        enddo
      endif

      if (state_following) then

        integer ::  state(N_st), idx
        double precision :: omax
        logical :: used
        logical, allocatable :: ok(:)
        complex*16, allocatable :: overlp(:,:)

        allocate(overlp(shift2,N_st),ok(shift2))

        overlp = (0d0,0d0)
        do j = 1, shift2-1, N_st_diag

          ! Computes some states from the guess vectors
          ! Psi(:,j:j+N_st_diag) = U y(:,j:j+N_st_diag) and put them
          ! in U(1,shift2+1:shift2+1+N_st_diag) as temporary array
          call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), U, size(U,1), y(1,j), size(y,1), (0.d0,0d0), U(1,shift2+1), size(U,1))
          do k = 1, N_st_diag
            if (cs_debug_state_following) then
              ! Norms and Hermitian overlaps of the Ritz vector before its c-normalization
              call inner_prod_c(U(1,shift2+k), U(1,shift2+k), sze, dbg_c)
              dbg_z(k+j-1) = cdabs(dbg_c)
              call inner_product_complex(U(1,shift2+k), U(1,shift2+k), sze, dbg_c)
              dbg_hn(k+j-1) = dble(dbg_c)
              do l = 1, N_st
                dbg_oh(k+j-1,l) = dbg_sf_ovh(U(1,shift2+k), u_in(1,l), sze)
                dbg_o0(k+j-1,l) = dbg_sf_ovh(U(1,shift2+k), dbg_u0(1,l), sze)
              enddo
            endif
            ! The c-normalization divides the overlap by sqrt(|psi^T psi|) and favours
            ! spurious vectors with dispersed phases. With cs_sf_hermitian_overlap, the
            ! Ritz vectors are normalized with the Hermitian norm, consistent with the
            ! Hermitian overlap below (a no-op up to round-off: y^H (U^H U) y = 1 above).
            if (cs_sf_hermitian_overlap) then
              call normalize_complex(U(:,shift2+k),sze)
            else
              call normalize_c(U(:,shift2+k),sze)
            endif
          enddo

          ! Overlap using the Hermitian product for state following. Reference:
          ! the input vectors of the call (cs_sf_fixed_ref) or the vectors of
          ! the last restart (u_in)
          if (cs_sf_fixed_ref) then
            do l = 1, N_st
              do k = 1, N_st_diag
                do i = 1, sze
                  overlp(k+j-1,l) += u_ref(i,l) * DCONJG(U(i,shift2+k))
                enddo
              enddo
            enddo
          else
            do l = 1, N_st
              do k = 1, N_st_diag
                do i = 1, sze
                  overlp(k+j-1,l) += u_in(i,l) * DCONJG(U(i,shift2+k))
                enddo
              enddo
            enddo
          endif

        enddo

        ! state_ok is overwritten by the selection below
        if (cs_debug_state_following) dbg_ok(1:shift2) = state_ok(1:shift2)

        ! cs_sf_guard: Ritz vectors with Im(E_total) > cs_sf_im_tol are not physical
        integer :: pass
        logical :: accept
        logical, allocatable :: im_bad(:)
        allocate(im_bad(shift2))
        do k = 1, shift2
          im_bad(k) = cs_sf_guard .and.                                              &
              (dimag(lambda(k) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0)) > cs_sf_im_tol)
        enddo

        state = 0
        do l = 1, N_st

          ! Pass 1: S2, Im(E) (cs_sf_guard) and energy window (cs_sf_window)
          ! filters. If no candidate is acceptable, pass 2 drops the Im(E)
          ! filter, pass 3 the window and pass 4 the S2 filter, instead of
          ! leaving idx = 0. Without cs_sf_guard and cs_sf_window, only pass 1.
          do pass = 1, 4
            omax = 0d0
            idx = 0
            do k = 1, shift2

              ! Already used ?
              used = .False.
              do i = 1, N_st
                if (state(i) == k) then
                  used = .True.
                endif
              enddo

              win_ok = .True.
              if (window_on) then
                e_tot  = lambda(k) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0)
                win_ok = cdabs(e_tot - e_win_ref(l)) <= cs_sf_window
              endif

              if (pass == 1) then
                accept = state_ok(k) .and. (.not. im_bad(k)) .and. win_ok
              else if (pass == 2) then
                accept = state_ok(k) .and. win_ok
              else if (pass == 3) then
                accept = state_ok(k)
              else
                accept = .True.
              endif

              ! Maximum overlap
              if ((cdabs(overlp(k,l)) > omax) .and. (.not. used) .and. accept) then
                omax = cdabs(overlp(k,l))
                idx = k
              endif
            enddo
            if ((idx /= 0) .or. (.not. (cs_sf_guard .or. window_on))) exit
          enddo

          if (pass > 1) then
            if (pass == 2) then
              print *, 'cs_sf_guard: no candidate with Im(E) <= cs_sf_im_tol for state', l, &
                  '; Im(E) filter relaxed'
            else if (pass == 3) then
              print *, 'cs_sf_window: no candidate in the energy window for state', l,      &
                  '; window and Im(E) filters relaxed'
            else
              print *, 'cs_sf: no candidate with the expected S2 for state', l,              &
                  '; all filters relaxed'
            endif
          endif

          ! cs_dav_debug: the candidate of maximum overlap, without any filter,
          ! and the filters that rejected it when it is not the one chosen
          if (cs_dav_debug) then
            omax = 0d0
            idx0 = 0
            do k = 1, shift2
              used = .False.
              do i = 1, l-1
                if (state(i) == k) used = .True.
              enddo
              if ((cdabs(overlp(k,l)) > omax) .and. (.not. used)) then
                omax = cdabs(overlp(k,l))
                idx0 = k
              endif
            enddo
            if ((idx0 /= idx) .and. (idx0 > 0)) then
              dv_why = ''
              if (.not. state_ok(idx0)) then
                dv_why = trim(dv_why)//' S2'
                n_rej_s2 += 1
              endif
              if (im_bad(idx0)) then
                dv_why = trim(dv_why)//' Im'
                n_rej_im += 1
              endif
              if (window_on) then
                e_tot = lambda(idx0) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0)
                if (cdabs(e_tot - e_win_ref(l)) > cs_sf_window) then
                  dv_why = trim(dv_why)//' window'
                  n_rej_win += 1
                endif
              endif
              e_tot = lambda(idx0) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0)
              write(6,'(A,I4,A,I3,A,I2,A,2F16.10,A,F8.5,A,A)') ' DV reject cyc', itertot,  &
                  ' it', iter, ' st', l, '  max-overlap E =', e_tot,                        &
                  '  |ov| =', cdabs(overlp(idx0,l)), '  filters:', trim(dv_why)
              if (idx > 0) then
                e_tot = lambda(idx) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0)
                write(6,'(A,I4,A,I3,A,I2,A,2F16.10,A,F8.5)') ' DV chosen cyc', itertot,   &
                    ' it', iter, ' st', l, '  E =', e_tot, '  |ov| =', cdabs(overlp(idx,l))
              endif
            endif
          endif

          state(l) = idx
        enddo
        deallocate(im_bad)

        if (cs_debug_state_following) then
          do l = 1, N_st
            call dbg_sf_report(l, itertot, iter, sze, shift2, U, size(U,1), W, size(W,1), &
                y, size(y,1), lambda, s2, H_jj, overlp(1,l), dbg_z, dbg_hn,                 &
                dbg_oh(1,l), dbg_o0(1,l), dbg_ok, state(l), expected_s2,                   &
                e_itheta_cs * dcmplx(nuclear_repulsion, 0d0))
          enddo
        endif

        ! tmp array before setting state_ok
        ok = .False.
        do l = 1, N_st
          ok(state(l)) = .True.
        enddo

        do k = 1, shift2
          if (.not. ok(k)) then
            state_ok(k) = .False.
          endif
        enddo

        ! A candidate chosen after relaxing the filters must be moved to the front too
        if (cs_sf_guard) then
          do l = 1, N_st
            if (state(l) > 0) state_ok(state(l)) = .True.
          enddo
        endif

        deallocate(overlp,ok)
      endif

      do k=1,shift2
        if (.not. state_ok(k)) then
          do l=k+1,shift2
            if (state_ok(l)) then
              call zswap(shift2, y(1,k), 1, y(1,l), 1)
              call zswap(1, s2(k), 1, s2(l), 1)
              call zswap(1, lambda(k), 1, lambda(l), 1)
              state_ok(k) = .True.
              state_ok(l) = .False.
              exit
            endif
          enddo
        endif
      enddo

      ! Swapped eigenvectors
      prev_y = y

      ! cs_dav_thick_restart: Ritz vectors (in the basis U) of this iteration
      ! and of the previous one
      if (cs_dav_thick_restart) then
        y_keep_prev = y_keep
        y_keep = (0.d0,0.d0)
        y_keep(1:shift2,1:N_st_diag) = y(1:shift2,1:N_st_diag)
      endif

      ! Express eigenvectors of h in the determinant basis
      ! --------------------------------------------------
      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), U, size(U,1), y, size(y,1), (0.d0,0d0), U(1,shift2+1), size(U,1))

      call zgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.d0,0d0), W, size(W,1), y, size(y,1), (0.d0,0d0), W(1,shift2+1), size(W,1))

      y_s(:,:) = cmplx(real(dble(y(:,:))), real(dimag(y(:,:))))

      call cgemm('N','N', sze, N_st_diag, shift2,                    &
          (1.,0.), S, size(S,1), y_s, size(y_s,1), (0.,0.), S(1,shift2+1), size(S,1))

      ! Compute residual vector and davidson step
      ! -----------------------------------------

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag
        do i=1,sze
           if (dabs(dble(H_jj(i) - lambda (k))) >= 1d-2) then 
             U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k) ) /(H_jj(i) - lambda (k))
           else
             U(i,shift2+k) = (lambda(k) * U(i,shift2+k) - W(i,shift2+k) )      &
                / dcmplx(dsign(1d-2, dble(H_jj(i) - lambda (k))), dimag(H_jj(i) - lambda (k)))
           endif
        enddo
        if (k <= N_st) then
          ! Residual norm uses Hermitian inner product: always real and positive, reliable convergence check
          ! Total energy adds nuclear repulsion scaled by e^{-i*theta}
          call inner_product_complex(U(1,shift2+k),U(1,shift2+k),sze,residual_norm(k))
          to_print(1,k) = dble(lambda(k) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0))
          to_print(2,k) = dimag(lambda(k) + e_itheta_cs * dcmplx(nuclear_repulsion, 0d0))
          to_print(3,k) = dble(s2(k))
          to_print(4,k) = dimag(s2(k))
          to_print(5,k) = dble(residual_norm(k))
        endif
      enddo
      !$OMP END PARALLEL DO

      if ((itertot>1).and.(iter == 1)) then
        !don't print
        continue
      else
        write(*,'(1X,I3,1X,100(1X,F16.10,1X,F16.10,1X,ES10.2,1X,ES10.2,1X,ES10.2))') iter-1, to_print(1:5,1:N_st)
      endif

      ! Check convergence
      if (iter > 1) then
        if (threshold_davidson_from_pt2) then
          converged = dabs(maxval(abs(residual_norm(1:N_st)))) < threshold_davidson_pt2
        else
          converged = dabs(maxval(abs(residual_norm(1:N_st)))) < threshold_davidson
        endif
      endif

      do k=1,N_st
        if (dble(residual_norm(k)) > 1.d8) then
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

    ! cs_dav_thick_restart: previous Ritz vectors, kept in the next cycle
    ! (needs two iterations in this cycle; not after convergence)
    thick_next = cs_dav_thick_restart .and. (iter >= 2) .and. (.not. converged)
    if (thick_next) then
      call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),    &
          U, size(U,1), y_keep_prev, size(y_keep_prev,1), (0.d0,0d0), u_thick, size(u_thick,1))
    endif

    ! Re-contract U and update S and W
    ! --------------------------------

    call cgemm('N','N', sze, N_st_diag, shift2, (1.,0.),      &
        S, size(S,1), y_s, size(y_s,1), (0.,0.), S(1,shift2+1), size(S,1))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        S(i,k) = S(i,shift2+k)
      enddo
    enddo
    !$OMP END PARALLEL DO

    if (cs_debug_state_following) then
      do l = 1, N_st
        dbg_uold(1:sze,l) = u_in(1:sze,l)
      enddo
    endif

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        W, size(W,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        W(i,k) = u_in(i,k)
      enddo
    enddo
    !$OMP END PARALLEL DO

    call zgemm('N','N', sze, N_st_diag, shift2, (1.d0,0d0),      &
        U, size(U,1), y, size(y,1), (0.d0,0d0), u_in, size(u_in,1))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    !$OMP END PARALLEL DO

    call nullify_small_elements_complex(sze,N_st_diag,U,size(U,1),threshold_davidson_pt2)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
    do k=1,N_st_diag
      do i=1,sze
        U(i,k) = u_in(i,k)
      enddo
    enddo
    !$OMP END PARALLEL DO

    ! Thick restart: the next cycle starts with the current Ritz vectors
    ! (1:N_st_diag) and the previous ones (N_st_diag+1:2 N_st_diag); H and S2
    ! are applied to both at its first iteration
    if (thick_next) then
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k)
      do k=1,N_st_diag
        do i=1,sze
          U(i,N_st_diag+k) = u_thick(i,k)
        enddo
      enddo
      !$OMP END PARALLEL DO
      off = N_st_diag
      n_thick += 1
    else
      off = 0
    endif

    if (cs_dav_debug) then
      write(6,'(A,I4,A,I3,A,I5,A,L2,A,100ES10.2)') ' DV restart cyc', itertot, '  iterations', iter, &
          '  subspace', shift2, '  thick', thick_next, '  residuals', (dble(residual_norm(l)), l=1,N_st)
      do l = 1, N_st
        write(6,'(A,I4,A,I2,A,F10.7)') ' DV restart cyc', itertot, ' st', l,                         &
            '  |<u|u_input>| =', dbg_sf_ovh(u_in(1,l), u_ref(1,l), sze)
      enddo
    endif

    if (cs_debug_state_following) then
      ! Drift of the reference vector of the state following at this restart
      do l = 1, N_st
        write(6,'(A,I3,A,I2,2(A,F8.5))') ' SFR cyc', itertot, ' st', l,          &
            '  |<u_new|u_old>| =', dbg_sf_ovh(u_in(1,l), dbg_uold(1,l), sze),       &
            '  |<u_new|u_0>| =', dbg_sf_ovh(u_in(1,l), dbg_u0(1,l), sze)
      enddo
    endif

  enddo

  if (cs_debug_state_following) then
    deallocate(dbg_u0, dbg_uold, dbg_z, dbg_hn, dbg_ok, dbg_oh, dbg_o0)
  endif

  if (cs_dav_debug) then
    write(6,'(A,I4,A,L2,A,I4,A,3I6)') ' DV end  cycles', min(itertot, max(1,cs_dav_max_cycles)), &
        '  converged', converged, '  thick restarts', n_thick,                                &
        '  max-overlap rejected by S2 / Im / window', n_rej_s2, n_rej_im, n_rej_win
  endif
  if (allocated(u_ref)) deallocate(u_ref)
  if (allocated(u_thick)) deallocate(u_thick, y_keep, y_keep_prev)

  do k=1,N_st_diag
    energies(k) = lambda(k)
    s2_out(k) = s2(k)
  enddo
  write_buffer = '======'
  do i=1,N_st
    write_buffer = trim(write_buffer)//' ================ ================ ========== ========== =========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''
  call write_time(6)

  deallocate(W,S)

  deallocate (                                                       &
      residual_norm,                                                 &
      U, overlap,                                                    &
      h, y_s, S_d,                                                   &
      y, s_, s_tmp,                                                  &
      lambda, prev_y                                                 &
      )
  FREE nthreads_davidson
end
