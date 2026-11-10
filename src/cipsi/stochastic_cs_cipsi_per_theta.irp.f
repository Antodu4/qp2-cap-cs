
subroutine run_stochastic_cs_cipsi_per_theta(Ev,PT2)
  use selection_types
  implicit none
  BEGIN_DOC
! PROOF OF CONCEPT -- inverted Complex-Scaling CIPSI loop.
!
! Outer loop over the theta grid (theta_cs, theta_cs+theta_step_size, ...,
! n_steps_cs points). For EACH theta a complete CIPSI is converged: the
! determinant selection is driven by the CS PT2 at that theta, and the space
! is grown until the usual stopping criteria (N_det_max / pt2_max / ...).
! Once converged, the complex energy is recorded and we move to the next theta.
!
! This is the mirror image of run_stochastic_cs_cipsi, where a single
! determinant space (optimized at the base theta) is reused for every theta.
!
! For a faithful "complete CIPSI per theta", every theta restarts from the
! same bootstrap determinant set (see the RESET block below), with a real
! Davidson guess.
!
! WARM START (cs_warm_start_det / cs_warm_start_cfg): every theta after the
! first starts from the converged complex wave function of the previous theta,
! truncated to 1-cs_warm_start_tol of its weight, by determinants (det) or by
! whole configurations (cfg). The Davidson guess is the truncated complex
! wave function (already rotated), and the discarded determinants can be
! selected again at the new theta. cs_warm_start_debug prints diagnostics
! (lines starting with WS).
  END_DOC

  integer                        :: i,j,k,l
  double precision, intent(out)  :: Ev(N_states), PT2(N_states)
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop
  double precision, external     :: pt2_stop_value

  double precision               :: rss
  double precision, external     :: memory_of_double

  ! theta grid
  integer                        :: n_theta, k_theta
  double precision               :: theta0, dtheta
  complex*16, allocatable        :: e_theta(:,:)
  double precision               :: theta_val

  ! CS working arrays
  complex*16, allocatable        :: psi_cs(:,:), energy(:), corr(:)
  complex*16, allocatable        :: tmp_cs(:,:)
  integer(bit_kind), allocatable :: tmp_det(:,:,:)

  ! bootstrap (starting) determinant set, restored before each theta
  integer                        :: n_det_boot
  integer(bit_kind), allocatable :: det_boot(:,:,:)
  double precision, allocatable  :: coef_boot(:,:)

  ! remap bookkeeping
  integer(bit_kind), allocatable :: l_det(:,:,:)
  integer, allocatable           :: key(:), prev_key(:)
  integer                        :: idx
  logical                        :: ok

  ! warm start: truncated wave function of the previous theta
  logical                        :: warm
  integer                        :: n_keep, n_drop, n_resel
  logical, allocatable           :: keep(:)
  integer(bit_kind), allocatable :: det_keep(:,:,:), det_drop(:,:,:)
  complex*16, allocatable        :: cs_keep(:,:), guess(:,:)
  double precision               :: t_theta, t0, t1
  double precision, allocatable  :: ovl(:)
  integer, external              :: get_index_in_psi_det_sorted_bit

  complex*16                     :: e_itheta
  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error
  double precision               :: correlation_energy_ratio

  ! Trajectory collection: one record (row) per (state,theta,iter), printed at the end.
  ! rec columns: 1=iter 2=n_det 3=theta 4=state 5=ReE 6=ImE 7=ReE+PT2 8=ImE+PT2 9=RePT2 10=ImPT2 11=variance
  integer                        :: inner_iter, s, nrec, cap
  double precision, allocatable  :: rec(:,:), rec_tmp(:,:)

  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  relative_error = PT2_relative_error

  ! ---------------------------------------------------------------------------
  ! theta grid validation
  ! ---------------------------------------------------------------------------
  if (theta_cs < 0d0) then
    print*,'theta_cs must be >= 0d0 for a per-theta CS-CIPSI scan'
    call abort()
  endif
  if (n_steps_cs < 1) then
    print*,'n_steps_cs must be >= 1'
    call abort()
  endif
  if (n_steps_cs > 1 .and. theta_step_size <= 0d0) then
    print*,'Please set theta_step_size > 0d0 to screen a range of theta values'
    call abort()
  endif
  if (cs_warm_start_det .and. cs_warm_start_cfg) then
    print*,'cs_warm_start_det and cs_warm_start_cfg are exclusive'
    call abort()
  endif
  n_theta = n_steps_cs
  theta0  = theta_cs
  dtheta  = theta_step_size
  allocate(e_theta(N_states,n_theta))

  nrec = 0
  cap  = 0

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  ! ---------------------------------------------------------------------------
  ! Bootstrap: build an initial real |Psi> on the starting determinant set,
  ! then stash that determinant set so each theta can start a fresh CIPSI.
  ! ---------------------------------------------------------------------------
  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  call diagonalize_CI
  call save_wavefunction

  n_det_boot = N_det
  allocate(det_boot(N_int,2,n_det_boot), coef_boot(n_det_boot,N_states))
  det_boot(:,:,1:n_det_boot)        = psi_det (:,:,1:n_det_boot)
  coef_boot(1:n_det_boot,1:N_states) = psi_coef(1:n_det_boot,1:N_states)

  allocate (zeros(N_states))
  zeros = 0.d0
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  ! ===========================================================================
  ! OUTER LOOP OVER THETA
  ! ===========================================================================
  do k_theta = 1, n_theta

    theta_val = theta0 + dble(k_theta-1)*dtheta
    theta_cs  = theta_val
    touch theta_cs
    e_itheta = dcmplx(dcos(theta_cs), -dsin(theta_cs))

    write(*,'(A)') ''
    write(*,'(A)') '################################################################################'
    write(*,'(A,I6,A,I6,A,ES16.8)') ' ### CS-CIPSI per-theta -- step ', k_theta, ' /', n_theta, &
                                     ' -- theta = ', theta_cs
    write(*,'(A)') '################################################################################'

    call wall_time(t_theta)
    warm = (cs_warm_start_det .or. cs_warm_start_cfg) .and. (k_theta > 1)

    if (warm) then
      ! -----------------------------------------------------------------------
      ! WARM START from the previous theta: its converged space (tmp_det,
      ! tmp_cs), truncated to 1-cs_warm_start_tol of the weight, by
      ! determinants (cs_warm_start_det) or by configurations (cs_warm_start_cfg)
      ! -----------------------------------------------------------------------
      if (allocated(det_keep)) deallocate(det_keep,cs_keep,keep)
      allocate(det_keep(N_int,2,size(tmp_cs,1)), cs_keep(size(tmp_cs,1),N_states))
      allocate(keep(size(tmp_cs,1)))
      call cs_truncate_wf(size(tmp_cs,1), tmp_det, tmp_cs, cs_warm_start_tol, &
                          cs_warm_start_cfg, n_keep, det_keep, cs_keep, keep)
      write(*,'(A,A,A,I10,A,I10,A,ES10.2,A)') ' Warm start (', &
                          trim(merge('cfg','det',cs_warm_start_cfg)), '): ', n_keep, ' /', &
                          size(tmp_cs,1), ' determinants kept (discarded weight < ', &
                          cs_warm_start_tol, ')'
      if (cs_warm_start_debug) then
        call cs_ws_print_truncation(size(tmp_cs,1), tmp_det, tmp_cs, keep)
        ! discarded determinants, to count those selected again at convergence
        n_drop = size(tmp_cs,1) - n_keep
        if (allocated(det_drop)) deallocate(det_drop)
        allocate(det_drop(N_int,2,max(n_drop,1)))
        j = 0
        do i = 1, size(tmp_cs,1)
          if (.not. keep(i)) then
            j = j+1
            det_drop(:,:,j) = tmp_det(:,:,i)
          endif
        enddo
      endif
      N_det = n_keep
      psi_det (:,:,1:n_keep)        = det_keep(:,:,1:n_keep)
      psi_coef(1:n_keep,1:N_states) = dble(cs_keep(1:n_keep,1:N_states))
      SOFT_TOUCH N_det psi_det psi_coef
    else
      ! -----------------------------------------------------------------------
      ! RESET to the bootstrap determinant set => a complete CIPSI for each theta.
      ! -----------------------------------------------------------------------
      N_det = n_det_boot
      psi_det (:,:,1:n_det_boot)         = det_boot(:,:,1:n_det_boot)
      psi_coef(1:n_det_boot,1:N_states)  = coef_boot(1:n_det_boot,1:N_states)
      SOFT_TOUCH N_det psi_det psi_coef
    endif
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    call diagonalize_CI
    if (warm .and. cs_warm_start_debug) then
      write(*,'(A,I10,A,I10)') ' WS spin completion: N_det = ', N_det, &
                               '   added determinants = ', N_det - n_keep
    endif

    ! Initial CS wavefunction at this theta: the previous theta's complex
    ! coefficients (warm start) or the real coefficients
    if (allocated(psi_cs)) deallocate(psi_cs,energy,corr)
    allocate(psi_cs(N_det,N_states), energy(N_states), corr(N_states))
    if (warm) then
      call cs_remap_guess(n_keep, det_keep, cs_keep, psi_cs)
    else
      psi_cs = dcmplx(psi_coef,0d0)
    endif
    if (cs_warm_start_debug) then
      if (allocated(ovl)) deallocate(ovl)
      allocate(guess(N_det,N_states), ovl(N_states))
      guess = psi_cs
      call wall_time(t0)
    endif
    call diagonalize_ci_cs(psi_cs,energy,corr)   ! sets psi_cs_coef, psi_energy_cs
    if (cs_warm_start_debug) then
      call wall_time(t1)
      call cs_hermitian_overlap(N_det, guess, psi_cs, ovl)
      write(*,'(A,A,A,F12.2,A)') ' WS first CS diagonalization (', &
                                 trim(merge('warm','cold',warm)), ' guess): ', t1-t0, ' s'
      do s = 1, N_states
        write(*,'(A,I3,A,F12.8,A,2F18.10)') ' WS state ', s, '  |<guess|psi>| = ', ovl(s), &
                                            '   E = ', dble(energy(s)), dimag(energy(s))
      enddo
      deallocate(guess)
    endif

    if (allocated(tmp_cs)) deallocate(tmp_cs,tmp_det)
    allocate(tmp_cs(N_det,N_states), tmp_det(N_int,2,N_det))
    tmp_cs  = psi_cs_coef
    tmp_det = psi_det

    ! Reset PT2 bookkeeping so the inner loop starts fresh for this theta
    correlation_energy_ratio = 0.d0
    pt2_data % pt2      = -huge(1.e0)
    pt2_data % pt2_im   = -huge(1.e0)
    pt2_data % rpt2     = -huge(1.e0)
    pt2_data % overlap  = 0.d0
    pt2_data % variance = huge(1.e0)

    inner_iter = 0

    ! -------------------------------------------------------------------------
    ! INNER CIPSI LOOP at fixed theta
    ! -------------------------------------------------------------------------
    do while (                                                         &
          (N_det < N_det_max) .and.                                    &
          (pt2_stop_value(pt2_data) > pt2_max) .and.                   &
          (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and. &
          (correlation_energy_ratio <= correlation_energy_ratio_max)   &
          )
      write(*,'(A)')  '--------------------------------------------------------------------------------'

      inner_iter = inner_iter + 1

      to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
      to_select = max(N_states_diag, to_select)

      Ev(1:N_states) = psi_energy_with_nucl_rep(1:N_states)

      ! Selection + stochastic PT2, driven by the CS wavefunction at this theta
      call build_psi_coef_cs_sorted(psi_cs_coef)

      call pt2_dealloc(pt2_data)
      call pt2_dealloc(pt2_data_err)
      call pt2_alloc(pt2_data, N_states)
      call pt2_alloc(pt2_data_err, N_states)
      call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,to_select)

      PT2(1:N_states) = pt2_data % pt2(1:N_states)
      correlation_energy_ratio = (psi_energy_with_nucl_rep(1) - hf_energy_ref)  / &
                      (psi_energy_with_nucl_rep(1) + pt2_data % rpt2(1) - hf_energy_ref)
      correlation_energy_ratio = min(1.d0,correlation_energy_ratio)

      call write_double(6,correlation_energy_ratio, 'Correlation ratio')
      call print_summary_cs(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)
      call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)
      call increment_n_iter(psi_energy_with_nucl_rep, pt2_data)
      call print_mol_properties()
      call write_cipsi_json(pt2_data,pt2_data_err)

      ! Collect one record per state: E_var and PT2 of the CURRENT determinant set
      if (nrec + N_states > cap) then
        cap = max(2*cap, nrec + N_states, 1024)
        allocate(rec_tmp(11,cap))
        if (nrec > 0) rec_tmp(:,1:nrec) = rec(:,1:nrec)
        call move_alloc(rec_tmp, rec)
      endif
      do s = 1, N_states
        nrec = nrec + 1
        rec(1,nrec)  = dble(inner_iter)
        rec(2,nrec)  = dble(N_det)
        rec(3,nrec)  = theta_cs
        rec(4,nrec)  = dble(s)
        rec(5,nrec)  = dble(energy(s))
        rec(6,nrec)  = dimag(energy(s))
        rec(7,nrec)  = dble(energy(s))  + pt2_data % pt2(s)
        rec(8,nrec)  = dimag(energy(s)) + pt2_data % pt2_im(s)
        rec(9,nrec)  = pt2_data % pt2(s)
        rec(10,nrec) = pt2_data % pt2_im(s)
        rec(11,nrec) = pt2_data % variance(s)
      enddo

      ! build_psi_coef_cs_sorted overwrote the real psi_coef; restore it
      psi_coef = psi_coef_tmpsave
      touch psi_coef

      if (qp_stop()) exit

      ! Add selected determinants
      call copy_H_apply_buffer_to_wf()
      if (save_wf_after_selection) then
        call save_wavefunction
      endif

      PROVIDE  psi_coef
      PROVIDE  psi_det
      PROVIDE  psi_det_sorted

      call diagonalize_CI
      call save_wavefunction
      call save_energy(psi_energy_with_nucl_rep, zeros)

      ! Remap the previous CS coefficients onto the new determinant ordering
      ! (new determinants are zero-padded); reuse as the Davidson guess.
      allocate(l_det(N_int,2,N_det),key(N_det))
      allocate(prev_key(size(tmp_cs,1)))

      l_det = psi_det
      do i = 1, N_det
        key(i) = i
      enddo
      idx = 1
      call recursive_int_sort(l_det,key,N_det,N_int,idx)

      do i = 1, size(tmp_cs,1)
        prev_key(i) = i
      enddo
      idx = 1
      call recursive_int_sort(tmp_det,prev_key,size(tmp_cs,1),N_int,idx)

      psi_cs_coef = (0d0,0d0)
      idx = 1
      do i = 1, size(tmp_cs,1)
        do l = idx, N_det
          ok = .True.
          do j = 1, 2
            do k = 1, N_int
              if (l_det(k,j,l) /= tmp_det(k,j,i)) then
                ok = .False.
                exit
              endif
            enddo
          enddo
          if (ok) then
            psi_cs_coef(key(l),:) = tmp_cs(prev_key(i),:)
            idx = l+1
            exit
          endif
        enddo
      enddo

      deallocate(key,prev_key,l_det)
      touch psi_cs_coef

      ! Diagonalize H_CS(theta) on the enlarged space
      deallocate(psi_cs)
      allocate(psi_cs(N_det,N_states))
      psi_cs = psi_cs_coef
      call diagonalize_ci_cs(psi_cs,energy,corr)

      deallocate(tmp_cs,tmp_det)
      allocate(tmp_cs(N_det,N_states), tmp_det(N_int,2,N_det))
      tmp_cs  = psi_cs_coef
      tmp_det = psi_det

      if (qp_stop()) exit
    enddo

    ! -------------------------------------------------------------------------
    ! Final CS PT2 on the converged determinant set (no selection), so the
    ! converged point carries a PT2 like every other iteration. Same CS path
    ! as the final-PT2 block of run_stochastic_cs_cipsi.
    ! -------------------------------------------------------------------------
    if (.not. qp_stop()) then
      call pt2_dealloc(pt2_data)
      call pt2_dealloc(pt2_data_err)
      call pt2_alloc(pt2_data, N_states)
      call pt2_alloc(pt2_data_err, N_states)
      call build_psi_coef_cs_sorted(psi_cs_coef)
      call ZMQ_pt2(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0)
      psi_coef = psi_coef_tmpsave
      touch psi_coef
    endif

    ! Record the converged complex energy at this theta
    do k = 1, N_states
      e_theta(k,k_theta) = psi_energy_cs(k) + e_itheta * dcmplx(nuclear_repulsion,0d0)
    enddo

    ! Collect the converged point (E_var + final CS PT2) for each state
    inner_iter = inner_iter + 1
    if (nrec + N_states > cap) then
      cap = max(2*cap, nrec + N_states, 1024)
      allocate(rec_tmp(11,cap))
      if (nrec > 0) rec_tmp(:,1:nrec) = rec(:,1:nrec)
      call move_alloc(rec_tmp, rec)
    endif
    do s = 1, N_states
      nrec = nrec + 1
      rec(1,nrec)  = dble(inner_iter)
      rec(2,nrec)  = dble(N_det)
      rec(3,nrec)  = theta_cs
      rec(4,nrec)  = dble(s)
      rec(5,nrec)  = dble(e_theta(s,k_theta))
      rec(6,nrec)  = dimag(e_theta(s,k_theta))
      rec(7,nrec)  = dble(e_theta(s,k_theta))  + pt2_data % pt2(s)
      rec(8,nrec)  = dimag(e_theta(s,k_theta)) + pt2_data % pt2_im(s)
      rec(9,nrec)  = pt2_data % pt2(s)
      rec(10,nrec) = pt2_data % pt2_im(s)
      rec(11,nrec) = pt2_data % variance(s)
    enddo

    write(*,'(A,ES16.8,A,I10)') ' Converged theta = ', theta_cs, '   N_det = ', N_det

    if (cs_warm_start_debug) then
      call wall_time(t1)
      write(*,'(A,ES12.4,A,I4,A,I10,A,F12.2,A)') ' WS converged theta = ', theta_cs, &
            '   CIPSI iterations = ', inner_iter-1, '   N_det = ', N_det, &
            '   time = ', t1-t_theta, ' s'
      if (warm) then
        PROVIDE psi_det_sorted_bit
        n_resel = 0
        do i = 1, n_drop
          if (get_index_in_psi_det_sorted_bit(det_drop(1,1,i), N_int) > 0) then
            n_resel = n_resel + 1
          endif
        enddo
        write(*,'(A,I10,A,I10,A,I10)') ' WS start space = ', n_keep, &
              '   discarded = ', n_drop, '   discarded and selected again = ', n_resel
      endif
    endif

    if (qp_stop()) exit
  enddo

  ! ===========================================================================
  ! Final unified output (shared routine, identical layout for both algorithms).
  ! ===========================================================================
  if (.not. allocated(rec)) allocate(rec(11,1))
  call cs_print_traj('per_theta', n_theta, theta0, dtheta, N_states, N_det_max, pt2_max, nrec, rec, e_theta)

  ! Restore the base theta
  theta_cs = theta0
  touch theta_cs

  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  deallocate(e_theta, det_boot, coef_boot, zeros)
  if (allocated(psi_cs)) deallocate(psi_cs,energy,corr)
  if (allocated(tmp_cs)) deallocate(tmp_cs,tmp_det)
  if (allocated(det_keep)) deallocate(det_keep,cs_keep,keep)
  if (allocated(det_drop)) deallocate(det_drop)
  if (allocated(ovl)) deallocate(ovl)

end


subroutine cs_det_weights(n, coef, w)
  implicit none
  BEGIN_DOC
! Weight of each determinant in a complex wave function, averaged over the
! states: w_i = 1/N_states sum_s |c_is|^2 / <c_s|c_s> (Hermitian norm).
  END_DOC
  integer, intent(in)            :: n
  complex*16, intent(in)         :: coef(n,N_states)
  double precision, intent(out)  :: w(n)
  double precision               :: norm
  integer                        :: i, s

  w = 0.d0
  do s = 1, N_states
    norm = 0.d0
    do i = 1, n
      norm += dble(dconjg(coef(i,s))*coef(i,s))
    enddo
    do i = 1, n
      w(i) += dble(dconjg(coef(i,s))*coef(i,s)) / (norm*dble(N_states))
    enddo
  enddo

end


subroutine cs_det_to_cfg(n, det, cfg, n_cfg)
  use bitmasks
  implicit none
  BEGIN_DOC
! Configuration index (1..n_cfg) of each determinant of det: determinants
! with the same spatial occupation (singly and doubly occupied orbitals)
! belong to the same configuration.
  END_DOC
  integer, intent(in)            :: n
  integer(bit_kind), intent(in)  :: det(N_int,2,n)
  integer, intent(out)           :: cfg(n), n_cfg
  integer(bit_kind), allocatable :: occ(:,:,:)
  integer, allocatable           :: key(:)
  integer                        :: i, k

  allocate(occ(N_int,2,n), key(n))
  do i = 1, n
    do k = 1, N_int
      occ(k,1,i) = ieor(det(k,1,i),det(k,2,i))
      occ(k,2,i) = iand(det(k,1,i),det(k,2,i))
    enddo
    key(i) = i
  enddo
  call recursive_int_sort(occ,key,n,N_int,1)

  n_cfg = 0
  do i = 1, n
    if (i == 1) then
      n_cfg = 1
    else if (any(occ(:,:,i) /= occ(:,:,i-1))) then
      n_cfg = n_cfg + 1
    endif
    cfg(key(i)) = n_cfg
  enddo
  deallocate(occ, key)

end


subroutine cs_truncate_wf(n_in, det_in, coef_in, tol, by_cfg, n_out, det_out, coef_out, keep)
  use bitmasks
  implicit none
  BEGIN_DOC
! Truncates a complex wave function to 1-tol of its weight (cs_det_weights),
! discarding the smallest weights first:
!   by_cfg = .False. : determinant by determinant;
!   by_cfg = .True.  : whole configurations (weight of a configuration = sum of
!                      the weights of its determinants), so the truncated space
!                      stays spin-complete.
! At least N_states_diag determinants are kept. keep(i) tells whether the
! determinant i is kept; the kept determinants keep their original order.
  END_DOC
  integer, intent(in)            :: n_in
  integer(bit_kind), intent(in)  :: det_in(N_int,2,n_in)
  complex*16, intent(in)         :: coef_in(n_in,N_states)
  double precision, intent(in)   :: tol
  logical, intent(in)            :: by_cfg
  integer, intent(out)           :: n_out
  integer(bit_kind), intent(out) :: det_out(N_int,2,n_in)
  complex*16, intent(out)        :: coef_out(n_in,N_states)
  logical, intent(out)           :: keep(n_in)

  double precision, allocatable  :: w(:), wu(:)
  integer, allocatable           :: cfg(:), iorder(:), n_det_u(:)
  logical, allocatable           :: keep_u(:)
  double precision               :: wsum
  integer                        :: i, n_u, n_det_kept

  allocate(w(n_in), cfg(n_in))
  call cs_det_weights(n_in, coef_in, w)

  ! Units of the truncation: determinants or configurations
  if (by_cfg) then
    call cs_det_to_cfg(n_in, det_in, cfg, n_u)
  else
    n_u = n_in
    do i = 1, n_in
      cfg(i) = i
    enddo
  endif
  allocate(wu(n_u), iorder(n_u), keep_u(n_u), n_det_u(n_u))
  wu = 0.d0
  n_det_u = 0
  do i = 1, n_in
    wu(cfg(i)) -= w(i)                  ! -weight: dsort gives decreasing weight
    n_det_u(cfg(i)) += 1
  enddo
  do i = 1, n_u
    iorder(i) = i
  enddo
  call dsort(wu, iorder, n_u)

  ! Keep the heaviest units until 1-tol of the weight (and N_states_diag dets)
  keep_u = .False.
  wsum = 0.d0
  n_det_kept = 0
  do i = 1, n_u
    if ( (wsum >= 1.d0 - tol) .and. (n_det_kept >= N_states_diag) ) exit
    keep_u(iorder(i)) = .True.
    wsum = wsum - wu(i)
    n_det_kept = n_det_kept + n_det_u(iorder(i))
  enddo

  n_out = 0
  do i = 1, n_in
    keep(i) = keep_u(cfg(i))
    if (keep(i)) then
      n_out = n_out + 1
      det_out(:,:,n_out) = det_in(:,:,i)
      coef_out(n_out,:)  = coef_in(i,:)
    endif
  enddo
  deallocate(w, wu, cfg, iorder, keep_u, n_det_u)

end


subroutine cs_ws_print_truncation(n, det, coef, keep)
  use bitmasks
  implicit none
  BEGIN_DOC
! Debug of the warm-start truncation (lines starting with WS): determinants and
! configurations kept, configurations cut (partly kept, only possible for a
! truncation by determinants), discarded weight per state.
  END_DOC
  integer, intent(in)            :: n
  integer(bit_kind), intent(in)  :: det(N_int,2,n)
  complex*16, intent(in)         :: coef(n,N_states)
  logical, intent(in)            :: keep(n)
  integer, allocatable           :: cfg(:), n_in_cfg(:), n_kept_cfg(:)
  double precision               :: norm, wdrop
  integer                        :: i, s, n_cfg

  allocate(cfg(n))
  call cs_det_to_cfg(n, det, cfg, n_cfg)
  allocate(n_in_cfg(n_cfg), n_kept_cfg(n_cfg))
  n_in_cfg = 0
  n_kept_cfg = 0
  do i = 1, n
    n_in_cfg(cfg(i)) += 1
    if (keep(i)) n_kept_cfg(cfg(i)) += 1
  enddo
  write(*,'(A,I10,A,I10)') ' WS truncation: determinants kept = ', count(keep), ' /', n
  write(*,'(A,I10,A,I10,A,I10)') ' WS truncation: configurations kept = ', &
        count(n_kept_cfg == n_in_cfg), ' /', n_cfg, '   cut = ', &
        count((n_kept_cfg > 0) .and. (n_kept_cfg < n_in_cfg))
  do s = 1, N_states
    norm  = 0.d0
    wdrop = 0.d0
    do i = 1, n
      norm = norm + dble(dconjg(coef(i,s))*coef(i,s))
      if (.not. keep(i)) wdrop = wdrop + dble(dconjg(coef(i,s))*coef(i,s))
    enddo
    write(*,'(A,I3,A,ES12.4)') ' WS truncation: state ', s, '   discarded weight = ', wdrop/norm
  enddo
  deallocate(cfg, n_in_cfg, n_kept_cfg)

end


subroutine cs_remap_guess(n_old, det_old, coef_old, u)
  use bitmasks
  implicit none
  BEGIN_DOC
! Maps the complex coefficients of the determinants det_old onto the current
! psi_det (determinants absent from det_old get 0, e.g. those added by
! make_s2_eigenfunction) and normalizes each state with the Hermitian norm.
! Used as the guess of the complex Davidson.
  END_DOC
  integer, intent(in)            :: n_old
  integer(bit_kind), intent(in)  :: det_old(N_int,2,n_old)
  complex*16, intent(in)         :: coef_old(n_old,N_states)
  complex*16, intent(out)        :: u(N_det,N_states)

  integer, external              :: get_index_in_psi_det_sorted_bit
  integer(bit_kind), allocatable :: det_srt(:,:,:)
  double precision, allocatable  :: pos(:,:), pos_srt(:,:)
  double precision               :: norm
  integer                        :: i, k, s

  ! psi_det index of each determinant of psi_det_sorted_bit: sort the indices
  ! with the same routine as the psi_det_sorted_bit provider
  allocate(det_srt(N_int,2,N_det), pos(N_det,1), pos_srt(N_det,1))
  do i = 1, N_det
    pos(i,1) = dble(i)
  enddo
  call sort_dets_by_det_search_key(N_det, psi_det, pos, N_det, det_srt, pos_srt, 1)
  PROVIDE psi_det_sorted_bit

  u = (0.d0,0.d0)
  do i = 1, n_old
    k = get_index_in_psi_det_sorted_bit(det_old(1,1,i), N_int)
    if (k > 0) then
      u(nint(pos_srt(k,1)),:) = coef_old(i,:)
    endif
  enddo

  do s = 1, N_states
    norm = 0.d0
    do i = 1, N_det
      norm += dble(dconjg(u(i,s))*u(i,s))
    enddo
    u(:,s) = u(:,s) / dsqrt(norm)
  enddo
  deallocate(det_srt, pos, pos_srt)

end


subroutine cs_hermitian_overlap(n, u, v, ovl)
  implicit none
  BEGIN_DOC
! |<u_s|v_s>| / (||u_s|| ||v_s||) for each state, with the Hermitian product.
  END_DOC
  integer, intent(in)            :: n
  complex*16, intent(in)         :: u(n,N_states), v(n,N_states)
  double precision, intent(out)  :: ovl(N_states)
  complex*16                     :: uv
  double precision               :: uu, vv
  integer                        :: i, s

  do s = 1, N_states
    uv = (0.d0,0.d0)
    uu = 0.d0
    vv = 0.d0
    do i = 1, n
      uv = uv + dconjg(u(i,s))*v(i,s)
      uu = uu + dble(dconjg(u(i,s))*u(i,s))
      vv = vv + dble(dconjg(v(i,s))*v(i,s))
    enddo
    ovl(s) = abs(uv) / dsqrt(uu*vv)
  enddo

end
