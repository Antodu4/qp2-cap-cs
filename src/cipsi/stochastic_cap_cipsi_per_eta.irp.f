
subroutine run_stochastic_cap_cipsi_per_eta(Ev,PT2)
  use selection_types
  implicit none
  BEGIN_DOC
! Inverted CAP-CIPSI loop (mirror of run_stochastic_cs_cipsi_per_theta).
!
! Outer loop over the eta grid (eta_cap, eta_cap+eta_step_size, ...,
! n_steps_cap points). For EACH eta a complete CIPSI is converged: the
! determinant selection is driven by the CAP PT2 at that eta, and the space
! is grown until the usual stopping criteria (N_det_max / pt2_max / ...).
! Once converged, the complex energy is recorded and we move to the next eta.
!
! This is the mirror image of run_stochastic_cap_cipsi, where a single
! determinant space (optimized at the base eta) is reused for every eta.
!
! Cold start: every eta restarts from the same bootstrap determinant set
! (see the RESET block below), so each eta is an independent CIPSI.
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

  ! eta grid
  integer                        :: n_eta, k_eta
  double precision               :: eta0, deta
  complex*16, allocatable        :: e_eta(:,:)
  double precision               :: eta_val

  ! CAP working arrays
  complex*16, allocatable        :: psi_cap(:,:), energy(:), corr(:)
  complex*16, allocatable        :: tmp_cap(:,:)
  integer(bit_kind), allocatable :: tmp_det(:,:,:)

  ! bootstrap (starting) determinant set, restored before each eta
  integer                        :: n_det_boot
  integer(bit_kind), allocatable :: det_boot(:,:,:)
  double precision, allocatable  :: coef_boot(:,:)

  ! remap bookkeeping
  integer(bit_kind), allocatable :: l_det(:,:,:)
  integer, allocatable           :: key(:), prev_key(:)
  integer                        :: idx
  logical                        :: ok

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error
  double precision               :: correlation_energy_ratio

  ! Trajectory collection: one record (row) per (state,eta,iter), printed at the end.
  ! rec columns: 1=iter 2=n_det 3=eta 4=state 5=ReE 6=ImE 7=ReE+PT2 8=ImE+PT2 9=RePT2 10=ImPT2 11=variance
  ! (cap_rec instead of cap: a local named cap would shadow the subroutine cap())
  integer                        :: inner_iter, s, nrec, cap_rec
  double precision, allocatable  :: rec(:,:), rec_tmp(:,:)
  integer                        :: n_cipsi
  logical                        :: final_pt2

  ! Debug diagnostics of the PT2 (pt2_diag): values of the space on which the
  ! stopping test is evaluated, and one record per eta and state (CAPREC)
  logical                        :: stop_next
  double precision, allocatable  :: var_stop(:), abspt2_stop(:), sumim_stop(:), sabs_stop(:)
  double precision, allocatable  :: caprec(:,:)
  integer                        :: ncaprec

  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  relative_error = PT2_relative_error

  ! ---------------------------------------------------------------------------
  ! CAP parameters / eta grid validation
  ! ---------------------------------------------------------------------------
  if (cap_onset_x < 0d0 .or. cap_onset_y < 0d0 .or. cap_onset_z < 0d0) then
    print*,'cap_onset_x, cap_onset_y and cap_onset_z must be >= 0d0'
    call abort()
  endif
  if (eta_cap < 0d0) then
    print*,'eta_cap must be >= 0d0 for a per-eta CAP-CIPSI scan'
    call abort()
  endif
  if (n_steps_cap < 1) then
    print*,'n_steps_cap must be >= 1'
    call abort()
  endif
  if (n_steps_cap > 1 .and. eta_step_size <= 0d0) then
    print*,'Please set eta_step_size > 0d0 to screen a range of eta values'
    call abort()
  endif
  n_eta = n_steps_cap
  eta0  = eta_cap
  deta  = eta_step_size
  allocate(e_eta(N_states,n_eta))

  nrec    = 0
  cap_rec = 0

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  ! ---------------------------------------------------------------------------
  ! Bootstrap: build an initial real |Psi> on the starting determinant set,
  ! then stash that determinant set so each eta can start a fresh CIPSI.
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
  ! OUTER LOOP OVER ETA
  ! ===========================================================================
  ncaprec = 0
  if (pt2_diag) then
    allocate(var_stop(N_states), abspt2_stop(N_states), sumim_stop(N_states), sabs_stop(N_states))
    allocate(caprec(15, n_eta*N_states))
  endif

  do k_eta = 1, n_eta

    eta_val = eta0 + dble(k_eta-1)*deta
    eta_cap = eta_val
    touch eta_cap

    write(*,'(A)') ''
    write(*,'(A)') '################################################################################'
    write(*,'(A,I6,A,I6,A,ES16.8)') ' ### CAP-CIPSI per-eta -- step ', k_eta, ' /', n_eta, &
                                     ' -- eta = ', eta_cap
    write(*,'(A)') '################################################################################'

    ! -------------------------------------------------------------------------
    ! RESET to the bootstrap determinant set => a complete CIPSI for each eta
    ! (cold start).
    ! -------------------------------------------------------------------------
    N_det = n_det_boot
    psi_det (:,:,1:n_det_boot)         = det_boot(:,:,1:n_det_boot)
    psi_coef(1:n_det_boot,1:N_states)  = coef_boot(1:n_det_boot,1:N_states)
    SOFT_TOUCH N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    call diagonalize_CI

    ! Initial CAP wavefunction at this eta (guess from the real coefficients)
    if (allocated(psi_cap)) deallocate(psi_cap,energy,corr)
    allocate(psi_cap(N_det,N_states), energy(N_states), corr(N_states))
    psi_cap = dcmplx(psi_coef,0d0)
    call diagonalize_ci_cap(psi_cap,energy,corr)   ! sets psi_cap_coef, psi_energy_cap

    if (allocated(tmp_cap)) deallocate(tmp_cap,tmp_det)
    allocate(tmp_cap(N_det,N_states), tmp_det(N_int,2,N_det))
    tmp_cap = psi_cap_coef
    tmp_det = psi_det

    ! Reset PT2 bookkeeping so the inner loop starts fresh for this eta
    correlation_energy_ratio = 0.d0
    pt2_data % pt2      = -huge(1.e0)
    pt2_data % pt2_im   = -huge(1.e0)
    pt2_data % rpt2     = -huge(1.e0)
    pt2_data % overlap  = 0.d0
    pt2_data % variance = huge(1.e0)

    inner_iter = 0
    if (pt2_diag) then
      ! -1 : no CIPSI iteration at this eta (the loop was not entered)
      var_stop(:) = -1.d0 ; abspt2_stop(:) = -1.d0 ; sumim_stop(:) = -1.d0 ; sabs_stop(:) = -1.d0
    endif

    ! -------------------------------------------------------------------------
    ! INNER CIPSI LOOP at fixed eta
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

      ! Selection + stochastic PT2, driven by the CAP wavefunction at this eta
      call build_psi_coef_cap_sorted(psi_cap_coef)

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
      call print_summary_cap(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)
      if (pt2_diag) then
        call print_pt2_diag_cs(pt2_data, 'iter')
        ! Values of this space: the stopping test of the next loop head uses them
        do s = 1, N_states
          var_stop(s)    = pt2_data % variance(s)
          abspt2_stop(s) = dsqrt(pt2_data % pt2(s)**2 + pt2_data % pt2_im(s)**2)
          sumim_stop(s)  = pt2_data % pt2_norm(s)
          sabs_stop(s)   = pt2_data % diag(2,s)
        enddo
      endif
      call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)
      call increment_n_iter(psi_energy_with_nucl_rep, pt2_data)
      call print_mol_properties()
      call write_cipsi_json(pt2_data,pt2_data_err)

      ! Collect one record per state: E_var and PT2 of the CURRENT determinant set
      if (nrec + N_states > cap_rec) then
        cap_rec = max(2*cap_rec, nrec + N_states, 1024)
        allocate(rec_tmp(11,cap_rec))
        if (nrec > 0) rec_tmp(:,1:nrec) = rec(:,1:nrec)
        call move_alloc(rec_tmp, rec)
      endif
      do s = 1, N_states
        nrec = nrec + 1
        rec(1,nrec)  = dble(inner_iter)
        rec(2,nrec)  = dble(N_det)
        rec(3,nrec)  = eta_cap
        rec(4,nrec)  = dble(s)
        rec(5,nrec)  = dble(energy(s))
        rec(6,nrec)  = dimag(energy(s))
        rec(7,nrec)  = dble(energy(s))  + pt2_data % pt2(s)
        rec(8,nrec)  = dimag(energy(s)) + pt2_data % pt2_im(s)
        rec(9,nrec)  = pt2_data % pt2(s)
        rec(10,nrec) = pt2_data % pt2_im(s)
        rec(11,nrec) = pt2_data % variance(s)
      enddo

      ! build_psi_coef_cap_sorted overwrote the real psi_coef; restore it
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

      ! Remap the previous CAP coefficients onto the new determinant ordering
      ! (new determinants are zero-padded); reuse as the Davidson guess.
      allocate(l_det(N_int,2,N_det),key(N_det))
      allocate(prev_key(size(tmp_cap,1)))

      l_det = psi_det
      do i = 1, N_det
        key(i) = i
      enddo
      idx = 1
      call recursive_int_sort(l_det,key,N_det,N_int,idx)

      do i = 1, size(tmp_cap,1)
        prev_key(i) = i
      enddo
      idx = 1
      call recursive_int_sort(tmp_det,prev_key,size(tmp_cap,1),N_int,idx)

      psi_cap_coef = (0d0,0d0)
      idx = 1
      do i = 1, size(tmp_cap,1)
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
            psi_cap_coef(key(l),:) = tmp_cap(prev_key(i),:)
            idx = l+1
            exit
          endif
        enddo
      enddo

      deallocate(key,prev_key,l_det)
      touch psi_cap_coef

      ! Diagonalize H_CAP(eta) on the enlarged space
      deallocate(psi_cap)
      allocate(psi_cap(N_det,N_states))
      psi_cap = psi_cap_coef
      call diagonalize_ci_cap(psi_cap,energy,corr)

      deallocate(tmp_cap,tmp_det)
      allocate(tmp_cap(N_det,N_states), tmp_det(N_int,2,N_det))
      tmp_cap = psi_cap_coef
      tmp_det = psi_det

      if (pt2_diag) then
        ! Same condition as the loop head: the decision taken on the PT2 of the
        ! space before the enlargement that has just been done
        stop_next = .not. ( (N_det < N_det_max) .and.                               &
            (pt2_stop_value(pt2_data) > pt2_max) .and.                                 &
            (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and.        &
            (correlation_energy_ratio <= correlation_energy_ratio_max) )
        do s = 1, N_states
          write(*,'(A,ES11.4,A,I3,A,I2,A,I10,A,7(1X,ES12.4),A,L2)') ' STOPTEST eta', eta_cap, ' it', inner_iter, &
              ' st', s, ' Ndet', N_det,                                                &
              ' RePT2 ImPT2 |PT2| var sumIm S stopval :', pt2_data % pt2(s), pt2_data % pt2_im(s), &
              abspt2_stop(s), var_stop(s), sumim_stop(s), sabs_stop(s), pt2_stop_value(pt2_data), &
              '   stop =', stop_next
        enddo
      endif

      if (qp_stop()) exit
    enddo

    ! -------------------------------------------------------------------------
    ! End of the eta, as in run_stochastic_cap_cipsi (and the real CIPSI):
    ! the converged energy of the final space is recorded without PT2. A final
    ! PT2 (no selection) is computed only when the loop stopped because
    ! N_det > N_det_max, i.e. when no space passed the stopping test.
    ! -------------------------------------------------------------------------
    n_cipsi = inner_iter
    final_pt2 = (.not.qp_stop()) .and.                                  &
        (N_det > N_det_max) .and.                                       &
        (pt2_stop_value(pt2_data) > pt2_max) .and.                      &
        (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and. &
        (correlation_energy_ratio <= correlation_energy_ratio_max)
    if (final_pt2) then
      call pt2_dealloc(pt2_data)
      call pt2_dealloc(pt2_data_err)
      call pt2_alloc(pt2_data, N_states)
      call pt2_alloc(pt2_data_err, N_states)
      call build_psi_coef_cap_sorted(psi_cap_coef)
      call ZMQ_pt2(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0)
      psi_coef = psi_coef_tmpsave
      touch psi_coef
      call print_summary_cap(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)
    endif

    ! Record the converged complex energy at this eta (final space, no PT2;
    ! V_nn is not rotated with a CAP)
    do k = 1, N_states
      e_eta(k,k_eta) = psi_energy_cap(k) + dcmplx(nuclear_repulsion,0d0)
    enddo

    ! Trajectory record of the final space: only with a final PT2 (N_det_max)
    if (final_pt2) then
      inner_iter = inner_iter + 1
      if (nrec + N_states > cap_rec) then
        cap_rec = max(2*cap_rec, nrec + N_states, 1024)
        allocate(rec_tmp(11,cap_rec))
        if (nrec > 0) rec_tmp(:,1:nrec) = rec(:,1:nrec)
        call move_alloc(rec_tmp, rec)
      endif
      do s = 1, N_states
        nrec = nrec + 1
        rec(1,nrec)  = dble(inner_iter)
        rec(2,nrec)  = dble(N_det)
        rec(3,nrec)  = eta_cap
        rec(4,nrec)  = dble(s)
        rec(5,nrec)  = dble(e_eta(s,k_eta))
        rec(6,nrec)  = dimag(e_eta(s,k_eta))
        rec(7,nrec)  = dble(e_eta(s,k_eta))  + pt2_data % pt2(s)
        rec(8,nrec)  = dimag(e_eta(s,k_eta)) + pt2_data % pt2_im(s)
        rec(9,nrec)  = pt2_data % pt2(s)
        rec(10,nrec) = pt2_data % pt2_im(s)
        rec(11,nrec) = pt2_data % variance(s)
      enddo
    endif

    write(*,'(A,ES16.8,A,I10)') ' Converged eta = ', eta_cap, '   N_det = ', N_det

    if (pt2_diag) then
      if (final_pt2) then
        call print_pt2_diag_cs(pt2_data, 'final')
      endif
      ! Record of this eta: final space (E, and its PT2 only if computed),
      ! and space of the stopping decision
      do s = 1, N_states
        ncaprec = ncaprec + 1
        caprec(1,ncaprec)  = eta_cap
        caprec(2,ncaprec)  = dble(s)
        caprec(3,ncaprec)  = dble(N_det)
        caprec(4,ncaprec)  = dble(n_cipsi)
        caprec(5,ncaprec)  = dble(e_eta(s,k_eta))
        caprec(6,ncaprec)  = dimag(e_eta(s,k_eta))
        if (final_pt2) then
          caprec(7,ncaprec)  = dble(e_eta(s,k_eta))  + pt2_data % pt2(s)
          caprec(8,ncaprec)  = dimag(e_eta(s,k_eta)) + pt2_data % pt2_im(s)
          caprec(9,ncaprec)  = pt2_data % variance(s)
          caprec(10,ncaprec) = pt2_data % pt2_norm(s)
          caprec(11,ncaprec) = pt2_data % diag(2,s)
        else
          caprec(7,ncaprec)  = dble(e_eta(s,k_eta))
          caprec(8,ncaprec)  = dimag(e_eta(s,k_eta))
          caprec(9,ncaprec)  = -1.d0
          caprec(10,ncaprec) = -1.d0
          caprec(11,ncaprec) = -1.d0
        endif
        caprec(12,ncaprec) = var_stop(s)
        caprec(13,ncaprec) = abspt2_stop(s)
        caprec(14,ncaprec) = sumim_stop(s)
        caprec(15,ncaprec) = sabs_stop(s)
        call print_caprec(caprec(1,ncaprec))
      enddo
    endif

    if (qp_stop()) exit
  enddo

  ! ===========================================================================
  ! Final unified output (shared routine, identical layout for both algorithms).
  ! ===========================================================================
  if (.not. allocated(rec)) allocate(rec(11,1))
  call cap_print_traj('per_eta', n_eta, eta0, deta, N_states, N_det_max, pt2_max, nrec, rec, e_eta)

  if (pt2_diag) then
    write(*,'(A)') ''
    write(*,'(A)') ' CAPREC : one line per eta and state (debug pt2_diag). Final space: E, E+PT2, variance, sum|Im e|,'
    write(*,'(A)') ' CAPREC   S = sum|e| (no final PT2 unless N_det_max: E+PT2 = E and -1) ; space of the stopping'
    write(*,'(A)') ' CAPREC   decision: variance, |PT2|, sum|Im e|, S ; flag ok if the variance (final space, or else'
    write(*,'(A)') ' CAPREC   space of the stopping decision) is below 1e-3, suspect otherwise.'
    do i = 1, ncaprec
      call print_caprec(caprec(1,i))
    enddo
    deallocate(var_stop, abspt2_stop, sumim_stop, sabs_stop, caprec)
  endif

  ! Restore the base eta
  eta_cap = eta0
  touch eta_cap

  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  deallocate(e_eta, det_boot, coef_boot, zeros)
  if (allocated(psi_cap)) deallocate(psi_cap,energy,corr)
  if (allocated(tmp_cap)) deallocate(tmp_cap,tmp_det)

end
