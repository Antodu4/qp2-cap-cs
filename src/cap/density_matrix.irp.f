subroutine get_one_e_rdm_mo_cap(psi_coef_cap, one_e_dm_mo_cap)

  implicit none
  BEGIN_DOC
  ! Computes the total (alpha + beta) one-body density matrix in the MO basis
  ! for the CAP wave function using the c-inner product (no complex conjugation).
  !
  ! dm(p,q,istate) = sum_{spin} one_e_dm_mo_{spin}(p,q,istate)
  !
  ! Inputs:
  !   psi_coef_cap(N_det,N_states) : complex CI coefficients of the CAP wavefunction
  !
  ! Output:
  !   one_e_dm_mo_cap(mo_num,mo_num,N_states) : total one-body density matrix
  END_DOC

  complex*16, intent(in) :: psi_coef_cap(N_det,N_int)
  complex*16, intent(out) :: one_e_dm_mo_cap(mo_num,mo_num,N_states)

  complex*16, allocatable :: one_e_dm_mo_alpha_cap(:,:,:), one_e_dm_mo_beta_cap(:,:,:)
  
  allocate(one_e_dm_mo_alpha_cap(mo_num,mo_num,N_states))
  allocate(one_e_dm_mo_beta_cap(mo_num,mo_num,N_states))

  call get_one_e_rdm_mo_ab_cap(psi_coef_cap, one_e_dm_mo_alpha_cap, one_e_dm_mo_beta_cap)

  one_e_dm_mo_cap = one_e_dm_mo_alpha_cap + one_e_dm_mo_beta_cap

  deallocate(one_e_dm_mo_alpha_cap,one_e_dm_mo_beta_cap)

end

subroutine get_one_e_rdm_mo_ab_cap(psi_coef_cap, one_e_dm_mo_alpha_cap, one_e_dm_mo_beta_cap)

  implicit none
 
  BEGIN_DOC
  ! $\alpha$ and $\beta$ one-body density matrix for each state
  END_DOC

  complex*16, intent(in) :: psi_coef_cap(N_det,N_int)
  complex*16, intent(out) :: one_e_dm_mo_alpha_cap(mo_num,mo_num,N_states)
  complex*16, intent(out) :: one_e_dm_mo_beta_cap(mo_num,mo_num,N_states)

  integer                        :: j,k,l,m,k_a,k_b
  integer                        :: occ(N_int*bit_kind_size,2)
  complex*16               :: ck, cl, ckl
  double precision               :: phase
  integer                        :: h1,h2,p1,p2,s1,s2, degree
  integer(bit_kind)              :: tmp_det(N_int,2), tmp_det2(N_int)
  integer                        :: exc(0:2,2),n_occ(2)
  complex*16, allocatable  :: tmp_a(:,:,:), tmp_b(:,:,:)
  integer                        :: krow, kcol, lrow, lcol
  complex*16, allocatable :: psi_bilinear_matrix_values_cap(:,:)    
  complex*16, allocatable :: psi_bilinear_matrix_transp_values_cap(:,:)    

  PROVIDE psi_det

  allocate(psi_bilinear_matrix_values_cap(N_det,N_states))
  allocate(psi_bilinear_matrix_transp_values_cap(N_det,N_states))

  call get_psi_bilinear_matrix_values_cap(psi_coef_cap,psi_bilinear_matrix_values_cap)
  call get_psi_bilinear_matrix_transp_values_cap(psi_coef_cap,psi_bilinear_matrix_transp_values_cap) 

  one_e_dm_mo_alpha_cap = dcmplx(0.d0,0d0)
  one_e_dm_mo_beta_cap  = dcmplx(0.d0,0d0)
  !$OMP PARALLEL DEFAULT(NONE)                                      &
      !$OMP PRIVATE(j,k,k_a,k_b,l,m,occ,ck, cl, ckl,phase,h1,h2,p1,p2,s1,s2, degree,exc,&
      !$OMP  tmp_a, tmp_b, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)&
      !$OMP SHARED(psi_det,psi_coef,N_int,N_states,elec_alpha_num,  &
      !$OMP  elec_beta_num,one_e_dm_mo_alpha_cap,one_e_dm_mo_beta_cap,N_det,&
      !$OMP  mo_num,psi_bilinear_matrix_rows,psi_bilinear_matrix_columns,&
      !$OMP  psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns,&
      !$OMP  psi_bilinear_matrix_order_reverse, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP  psi_bilinear_matrix_values_cap, psi_bilinear_matrix_transp_values_cap,&
      !$OMP  N_det_alpha_unique,N_det_beta_unique,irp_here)
  allocate(tmp_a(mo_num,mo_num,N_states), tmp_b(mo_num,mo_num,N_states) )
  tmp_a = dcmplx(0.d0,0d0)
  !$OMP DO SCHEDULE(guided)
  do k_a=1,N_det
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
    do m=1,N_states
      ck = psi_bilinear_matrix_values_cap(k_a,m)*psi_bilinear_matrix_values_cap(k_a,m)
      do l=1,elec_alpha_num
        j = occ(l,1)
        tmp_a(j,j,m) += ck
      enddo
    enddo

    if (k_a == N_det) cycle
    l = k_a+1
    lrow = psi_bilinear_matrix_rows(l)
    lcol = psi_bilinear_matrix_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lcol == kcol )
      tmp_det2(:) = psi_det_alpha_unique(:, lrow)
      call get_excitation_degree_spin(tmp_det(1,1),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,1),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
        do m=1,N_states
          ckl = psi_bilinear_matrix_values_cap(k_a,m)*psi_bilinear_matrix_values_cap(l,m) * dcmplx(phase,0d0)
          tmp_a(h1,p1,m) += ckl
          tmp_a(p1,h1,m) += ckl
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_rows(l)
      lcol = psi_bilinear_matrix_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT

  !$OMP CRITICAL
  one_e_dm_mo_alpha_cap(:,:,:) = one_e_dm_mo_alpha_cap(:,:,:) + tmp_a(:,:,:)
  !$OMP END CRITICAL
  deallocate(tmp_a)

  tmp_b = dcmplx(0.d0,0d0)
  !$OMP DO SCHEDULE(guided)
  do k_b=1,N_det
    krow = psi_bilinear_matrix_transp_rows(k_b)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_transp_columns(k_b)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)

    ! Diagonal part
    ! -------------

    call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
    do m=1,N_states
      ck = psi_bilinear_matrix_transp_values_cap(k_b,m)*psi_bilinear_matrix_transp_values_cap(k_b,m)
      do l=1,elec_beta_num
        j = occ(l,2)
        tmp_b(j,j,m) += ck
      enddo
    enddo

    if (k_b == N_det) cycle
    l = k_b+1
    lrow = psi_bilinear_matrix_transp_rows(l)
    lcol = psi_bilinear_matrix_transp_columns(l)
    ! Fix beta determinant, loop over alphas
    do while ( lrow == krow )
      tmp_det2(:) = psi_det_beta_unique(:, lcol)
      call get_excitation_degree_spin(tmp_det(1,2),tmp_det2,degree,N_int)
      if (degree == 1) then
        exc = 0
        call get_single_excitation_spin(tmp_det(1,2),tmp_det2,exc,phase,N_int)
        call decode_exc_spin(exc,h1,p1,h2,p2)
        do m=1,N_states
          ckl = psi_bilinear_matrix_transp_values_cap(k_b,m)*psi_bilinear_matrix_transp_values_cap(l,m) * dcmplx(phase,0d0)
          tmp_b(h1,p1,m) += ckl
          tmp_b(p1,h1,m) += ckl
        enddo
      endif
      l = l+1
      if (l>N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l)
      lcol = psi_bilinear_matrix_transp_columns(l)
    enddo

  enddo
  !$OMP END DO NOWAIT
  !$OMP CRITICAL
  one_e_dm_mo_beta_cap(:,:,:)  = one_e_dm_mo_beta_cap(:,:,:)  + tmp_b(:,:,:)
  !$OMP END CRITICAL

  deallocate(tmp_b)
  !$OMP END PARALLEL

end

