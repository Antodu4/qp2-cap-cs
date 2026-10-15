
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
! same bootstrap determinant set (see the RESET block below). Switching to a
! warm start (carry the previous theta's space over) is a one-line change
! documented there.
  END_DOC

  integer                        :: i,j,k,l
  double precision, intent(out)  :: Ev(N_states), PT2(N_states)
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop

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

  complex*16                     :: e_itheta
  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error
  double precision               :: correlation_energy_ratio

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
  n_theta = n_steps_cs
  theta0  = theta_cs
  dtheta  = theta_step_size
  allocate(e_theta(N_states,n_theta))

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

    ! -------------------------------------------------------------------------
    ! RESET to the bootstrap determinant set => a complete CIPSI for each theta.
    ! For a warm start instead (reuse the previous theta's space), comment the
    ! three assignments + SOFT_TOUCH below and the diagonalize_CI that follows.
    ! -------------------------------------------------------------------------
    N_det = n_det_boot
    psi_det (:,:,1:n_det_boot)         = det_boot(:,:,1:n_det_boot)
    psi_coef(1:n_det_boot,1:N_states)  = coef_boot(1:n_det_boot,1:N_states)
    SOFT_TOUCH N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    call diagonalize_CI

    ! Initial CS wavefunction at this theta (guess from the real coefficients)
    if (allocated(psi_cs)) deallocate(psi_cs,energy,corr)
    allocate(psi_cs(N_det,N_states), energy(N_states), corr(N_states))
    psi_cs = dcmplx(psi_coef,0d0)
    call diagonalize_ci_cs(psi_cs,energy,corr)   ! sets psi_cs_coef, psi_energy_cs

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

    ! -------------------------------------------------------------------------
    ! INNER CIPSI LOOP at fixed theta
    ! -------------------------------------------------------------------------
    do while (                                                         &
          (N_det < N_det_max) .and.                                    &
          (maxval(abs(pt2_data % pt2(1:N_states))) > pt2_max) .and.    &
          (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and. &
          (correlation_energy_ratio <= correlation_energy_ratio_max)   &
          )
      write(*,'(A)')  '--------------------------------------------------------------------------------'

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

    ! Record the converged complex energy at this theta
    do k = 1, N_states
      e_theta(k,k_theta) = psi_energy_cs(k) + e_itheta * dcmplx(nuclear_repulsion,0d0)
    enddo

    write(*,'(A,ES16.8,A,I10)') ' Converged theta = ', theta_cs, '   N_det = ', N_det

    if (qp_stop()) exit
  enddo

  ! ===========================================================================
  ! Summary table E(theta) for each state
  ! ===========================================================================
  write(*,*) ''
  write(*,*) ' #####################################'
  write(*,*) ' ###  CS per-theta CIPSI energies  ###'
  write(*,*) ' #####################################'
  write(*,*) ''
  write(*,'(A15)',advance='no') '   Theta     '
  do k = 1, N_states
    write(*,'(A33)',advance='no') '   Energy (Re)     Energy (Im)  '
  enddo
  write(*,*) ''
  do k_theta = 1, n_theta
    theta_val = theta0 + dble(k_theta-1)*dtheta
    write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') theta_val, &
      (dble(e_theta(k,k_theta)), dimag(e_theta(k,k_theta)), k=1,N_states)
  enddo
  write(*,*) ''

  ! Restore the base theta
  theta_cs = theta0
  touch theta_cs

  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)
  deallocate(e_theta, det_boot, coef_boot, zeros)
  if (allocated(psi_cs)) deallocate(psi_cs,energy,corr)
  if (allocated(tmp_cs)) deallocate(tmp_cs,tmp_det)

end
