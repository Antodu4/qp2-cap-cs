subroutine mo_one_rdm_cap(coefs, N_st, sze, rdm)
  
  BEGIN_DOC
  ! Part not called by the code. Retained if necessary.
  END_DOC

  implicit none

  integer, intent(in) ::  N_st, sze
  complex*16, intent(in) :: coefs(sze,N_st)
  complex*16, intent(out) :: rdm(mo_num,mo_num,N_st)

  double precision :: phase
  complex*16 :: trace
  integer :: i,j,k,l,s,p1,h1,s1,p2,h2,s2
  integer :: degree,exc(0:2,2,2)
  integer :: n(2), occ(N_int*bit_kind_size,2)

  n(1) = elec_alpha_num
  n(2) = elec_beta_num

  do k = 1, N_st
    rdm(:,:,k) = 0d0
    do j = 1, N_det
      do i = 1, N_det
        call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
        if (degree == 0) then
          call bitstring_to_list_ab(psi_det(1,1,i), occ, n, N_int)
          do s = 1, 2
            do l = 1, n(s)
              rdm(occ(l,s),occ(l,s),k) += coefs(i,k)**2
            enddo
          enddo
        else if (degree == 1) then
          call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
          call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
          rdm(h1,p1,k) += phase * coefs(i,k) * coefs(j,k)
        endif
      enddo
    enddo
  enddo

  trace = (0d0,0d0)
  do i = 1, mo_num
    trace += rdm(i,i,1)
  enddo
  print*,'trace',trace
  print*,'rdm'
  do i = 1, mo_num
    write(*,'(100(ES12.4))') rdm(i,:,1)
  enddo

  complex*16, allocatable :: one_e_dm_mo_cap(:,:,:)

  allocate(one_e_dm_mo_cap(mo_num,mo_num,N_states))

  call get_one_e_rdm_mo_cap(coefs, one_e_dm_mo_cap)
  print*,'rdm'
  do i = 1, mo_num
    write(*,'(100(ES12.4))') one_e_dm_mo_cap(i,:,1)
  enddo
  print*,'qp rdm'
  do i = 1, mo_num
    write(*,'(100(ES12.4))') dm(i,:,1,1)
  enddo

end

BEGIN_PROVIDER [ double precision, dm, (mo_num,mo_num,N_states,N_states)]
  implicit none
  BEGIN_DOC
  ! Real one-body density matrix in the MO basis.
  !
  ! dm(p,q,istate,jstate) = sum_{i,j} C_{i,istate} * C_{j,jstate} * <i|E_pq|j>
  END_DOC

  integer :: i,j,istate,jstate
  integer :: p,q,tmp_p
  integer :: degree, exc(0:2,2,2)
  integer :: h1,p1,h2,p2,s1,s2
  double precision :: phase, max_error
  integer :: nb_error
  integer, allocatable :: list_occ_a(:), list_occ_b(:)

  allocate(list_occ_a(elec_alpha_num), list_occ_b(elec_beta_num))

  dm = 0d0

  do jstate = 1, N_states
    do istate = 1, N_states

      do i = 1, N_det
        do j = 1, N_det
          if (i == j) then

            ! diag part
            call bitstring_to_list(psi_det_alpha(N_int,i), list_occ_a, elec_alpha_num, N_int)
            call bitstring_to_list(psi_det_beta(N_int,i), list_occ_b, elec_beta_num, N_int)

            ! alpha
            do p = 1, elec_alpha_num
              tmp_p = list_occ_a(p)
              dm(tmp_p,tmp_p,istate,jstate) = dm(tmp_p,tmp_p,istate,jstate) + psi_coef(i,istate) * psi_coef(i,jstate)
            enddo

            ! beta
            do p = 1, elec_beta_num
              tmp_p = list_occ_b(p)
              dm(tmp_p,tmp_p,istate,jstate) = dm(tmp_p,tmp_p,istate,jstate) + psi_coef(i,istate) * psi_coef(i,jstate)
            enddo

          else

            ! extra diag part
            call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
            if (degree == 1) then
              exc = 0
              call get_single_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,phase,N_int)
              call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
              dm(h1,p1,istate,jstate) = dm(h1,p1,istate,jstate) + psi_coef(i,istate) * psi_coef(j,jstate) * phase
            endif
          endif
        enddo
      enddo

    enddo
  enddo

  deallocate(list_occ_a,list_occ_b)

END_PROVIDER
