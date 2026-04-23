
subroutine run_stochastic_cs_cipsi(Ev,PT2) 
  use selection_types
  implicit none
  BEGIN_DOC
! Selected Full Configuration Interaction with Stochastic selection and PT2.
! Uses Complex Scaling (CS) instead of CAP.
  END_DOC
  integer                        :: i,j,k,l
  double precision, intent(out)  :: Ev(N_states), PT2(N_states) 
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop


  double precision :: rss
  double precision, external :: memory_of_double
  complex*16, allocatable :: tmp_cs(:,:)
  integer(bit_kind), allocatable :: tmp_det(:,:,:)
  
  complex*16 :: e_itheta
  e_itheta = dcmplx(cos(theta_cs), -sin(theta_cs))

  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  allocate (zeros(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error

  relative_error=PT2_relative_error

  zeros = 0.d0
  pt2_data % pt2   = -huge(1.e0)
  pt2_data % pt2_im = -huge(1.e0)
  pt2_data % rpt2  = -huge(1.e0)
  pt2_data % overlap= 0.d0
  pt2_data % variance = huge(1.e0)

  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  call diagonalize_CI
  call save_wavefunction
  if (do_cs) then
    do k = 1, N_states
        do i = 1, N_det
            psi_cs_coef(i,k) = dcmplx(psi_coef(i,k),0d0)
        enddo
    enddo
    touch psi_cs_coef
    call cs()
    allocate(tmp_cs(N_det,N_states), tmp_det(N_int,2,N_det))
    tmp_cs = psi_cs_coef
    tmp_det = psi_det
  endif

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  if (N_det > N_det_max) then
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    call diagonalize_CI
    call save_wavefunction

    if (do_cs) then
      do i = 1, N_det
        psi_cs_coef(i,:) = tmp_cs(i,:)
      enddo
      touch psi_cs_coef

      call cs()

      deallocate(tmp_cs,tmp_det)
      allocate(tmp_cs(N_det,N_states), tmp_det(N_int,2,N_det))

      tmp_cs = psi_cs_coef
      tmp_det = psi_det
    endif

  endif

  double precision :: correlation_energy_ratio

  correlation_energy_ratio = 0.d0

  do while (                                                         &
        (N_det < N_det_max) .and.                                    &
        (maxval(abs(pt2_data % pt2(1:N_states))) > pt2_max) .and.               &
        (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and.     &
        (correlation_energy_ratio <= correlation_energy_ratio_max)      &
        )
      write(*,'(A)')  '--------------------------------------------------------------------------------'
    
    to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
    to_select = max(N_states_diag, to_select)

    Ev(1:N_states) = psi_energy_with_nucl_rep(1:N_states)

    if (do_cs .and. cs_pt2) then
        call build_psi_coef_cs_sorted(psi_cs_coef)
    endif
    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,to_select) ! Stochastic PT2 and selection

    PT2(1:N_states) = pt2_data % pt2(1:N_states)
    correlation_energy_ratio = (psi_energy_with_nucl_rep(1) - hf_energy_ref)  /     &
                    (psi_energy_with_nucl_rep(1) + pt2_data % rpt2(1) - hf_energy_ref)
    correlation_energy_ratio = min(1.d0,correlation_energy_ratio)

    call write_double(6,correlation_energy_ratio, 'Correlation ratio')

    if (do_cs .and. cs_pt2) then
      call print_summary_cs(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)
    else
      call print_summary(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)
    endif

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)

    call increment_n_iter(psi_energy_with_nucl_rep, pt2_data)
    if (.not. cs_pt2) then
        call print_extrapolated_energy()
    endif
    call print_mol_properties()
    call write_cipsi_json(pt2_data,pt2_data_err)

    if (do_cs) then
        psi_coef = psi_coef_tmpsave
        touch psi_coef
    endif

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
    if (do_cs) then
      integer(bit_kind), allocatable :: l_det(:,:,:), l_prev_det(:,:,:)
      integer, allocatable :: key(:), prev_key(:)
      integer :: idx
      logical :: ok

      allocate(l_det(N_int,2,N_det),key(N_det))
      allocate(prev_key(size(tmp_cs,1)))

      ! To match the new order of psi_det...
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
      call cs()

      deallocate(tmp_cs,tmp_det)
      allocate(tmp_cs(N_det,N_states), tmp_det(N_int,2,N_det))

      tmp_cs = psi_cs_coef
      tmp_det = psi_det
    endif
    if (qp_stop()) exit
  enddo

  ! If stopped because N_det > N_det_max, do an extra iteration to compute the PT2
  if ((.not.qp_stop()).and.                                          &
        (N_det > N_det_max) .and.                                    &
        (maxval(abs(pt2_data % pt2(1:N_states))) > pt2_max) .and.    &
        (maxval(abs(pt2_data % variance(1:N_states))) > variance_max) .and.&
        (correlation_energy_ratio <= correlation_energy_ratio_max)   &
        ) then
    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    if (do_cs .and. cs_pt2) then
        call build_psi_coef_cs_sorted(psi_cs_coef)
    endif
    if (do_cs .and. cs_pt2) then
      deallocate(tmp_cs)
    endif
    if (do_cs .and.  cs_pt2) then
      call ZMQ_pt2(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0) ! Stochastic PT2
    endif

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)
    if (do_cs .and. cs_pt2) then
      call print_summary_cs(psi_energy_with_nucl_rep, &
       pt2_data , pt2_data_err, N_det, N_configuration, N_states, psi_s2)
    else
      call print_summary(psi_energy_with_nucl_rep, &
       pt2_data , pt2_data_err, N_det, N_configuration, N_states, psi_s2)
    endif
    call increment_n_iter(psi_energy_with_nucl_rep, pt2_data)
    call print_extrapolated_energy()
    call print_mol_properties()
    call write_cipsi_json(pt2_data,pt2_data_err)
  endif
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)

end
