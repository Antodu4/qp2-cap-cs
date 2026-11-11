! Helpers of the warm start of the per-theta CS-CIPSI
! (run_stochastic_cs_cipsi_per_theta, keywords cs_warm_start and
! cs_warm_start_tol).
!
! At each theta after the first, the converged complex wave function of the
! previous theta is truncated (cs_truncate_wf) and its coefficients are mapped
! onto the new psi_det (cs_remap_guess) to be the guess of the complex
! Davidson.


subroutine cs_det_weights(n, coef, w)
  implicit none
  BEGIN_DOC
! Weight of each determinant in a complex wave function, averaged over the
! states:
!   w_i = 1/N_states sum_s |c_is|^2 / <c_s|c_s>
! with the Hermitian norm, so that the weights are real, positive and sum to
! one. The c-norm (sum c_is^2) is complex and cannot be used to rank weights.
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
! Configuration index (1..n_cfg) of each determinant of det. Determinants with
! the same spatial occupation (same singly and doubly occupied orbitals) belong
! to the same configuration: they differ only by their spin coupling and are
! all needed to build an eigenfunction of S^2.
  END_DOC
  integer, intent(in)            :: n
  integer(bit_kind), intent(in)  :: det(N_int,2,n)
  integer, intent(out)           :: cfg(n), n_cfg
  integer(bit_kind), allocatable :: occ(:,:,:)
  integer, allocatable           :: key(:)
  integer                        :: i, k

  ! Spatial occupation: singly occupied (xor) and doubly occupied (and) orbitals
  allocate(occ(N_int,2,n), key(n))
  do i = 1, n
    do k = 1, N_int
      occ(k,1,i) = ieor(det(k,1,i),det(k,2,i))
      occ(k,2,i) = iand(det(k,1,i),det(k,2,i))
    enddo
    key(i) = i
  enddo

  ! After sorting, the determinants of a configuration are contiguous
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


subroutine cs_truncate_wf(n_in, det_in, coef_in, tol, n_out, det_out, coef_out)
  use bitmasks
  implicit none
  BEGIN_DOC
! Truncates a complex wave function to 1-tol of its weight (cs_det_weights)
! by discarding whole configurations, the lightest first. The weight of a
! configuration is the sum of the weights of its determinants.
!
! Whole configurations are kept so that the truncated space stays
! spin-complete: if a configuration were cut, make_s2_eigenfunction would add
! its missing determinants back with a zero coefficient, and the guess of the
! complex Davidson would not be an eigenfunction of S^2. Such a spin
! contamination of the guess is what lets the state following of the Davidson
! drift to a spurious root.
!
! At least N_states_diag determinants are kept, so that the Davidson has
! enough vectors. The kept determinants keep their original order in det_out
! and coef_out.
  END_DOC
  integer, intent(in)            :: n_in
  integer(bit_kind), intent(in)  :: det_in(N_int,2,n_in)
  complex*16, intent(in)         :: coef_in(n_in,N_states)
  double precision, intent(in)   :: tol
  integer, intent(out)           :: n_out
  integer(bit_kind), intent(out) :: det_out(N_int,2,n_in)
  complex*16, intent(out)        :: coef_out(n_in,N_states)

  double precision, allocatable  :: w(:), w_cfg(:)
  integer, allocatable           :: cfg(:), iorder(:), n_det_cfg(:)
  logical, allocatable           :: keep_cfg(:)
  double precision               :: wsum
  integer                        :: i, n_cfg, n_det_kept

  allocate(w(n_in), cfg(n_in))
  call cs_det_weights(n_in, coef_in, w)
  call cs_det_to_cfg(n_in, det_in, cfg, n_cfg)

  ! Weight and number of determinants of each configuration. The weights are
  ! stored with a minus sign so that dsort (increasing order) ranks the
  ! configurations by decreasing weight.
  allocate(w_cfg(n_cfg), iorder(n_cfg), keep_cfg(n_cfg), n_det_cfg(n_cfg))
  w_cfg = 0.d0
  n_det_cfg = 0
  do i = 1, n_in
    w_cfg(cfg(i)) -= w(i)
    n_det_cfg(cfg(i)) += 1
  enddo
  do i = 1, n_cfg
    iorder(i) = i
  enddo
  call dsort(w_cfg, iorder, n_cfg)

  ! Keep the heaviest configurations until 1-tol of the weight is reached
  ! (and at least N_states_diag determinants)
  keep_cfg = .False.
  wsum = 0.d0
  n_det_kept = 0
  do i = 1, n_cfg
    if ( (wsum >= 1.d0 - tol) .and. (n_det_kept >= N_states_diag) ) exit
    keep_cfg(iorder(i)) = .True.
    wsum = wsum - w_cfg(i)
    n_det_kept = n_det_kept + n_det_cfg(iorder(i))
  enddo

  ! Determinants of the kept configurations, in their original order
  n_out = 0
  do i = 1, n_in
    if (keep_cfg(cfg(i))) then
      n_out = n_out + 1
      det_out(:,:,n_out) = det_in(:,:,i)
      coef_out(n_out,:)  = coef_in(i,:)
    endif
  enddo
  deallocate(w, w_cfg, cfg, iorder, keep_cfg, n_det_cfg)

end


subroutine cs_remap_guess(n_old, det_old, coef_old, u)
  use bitmasks
  implicit none
  BEGIN_DOC
! Maps the complex coefficients of the determinants det_old onto the current
! psi_det, and normalizes each state with the Hermitian norm. Determinants of
! psi_det absent from det_old (e.g. those added by make_s2_eigenfunction, if any) get
! a zero coefficient. Used as the guess of the complex Davidson.
!
! psi_det may have been reordered by make_s2_eigenfunction, so the
! determinants are found with get_index_in_psi_det_sorted_bit, and the
! position in psi_det of each entry of psi_det_sorted_bit is recovered by
! sorting the indices 1..N_det with the same routine as the provider of
! psi_det_sorted_bit.
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

  ! pos_srt(k,1): index in psi_det of the k-th determinant of psi_det_sorted_bit
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
