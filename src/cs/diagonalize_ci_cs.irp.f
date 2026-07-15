subroutine diagonalize_ci_cs(u_in, energy, corr)

   implicit none
   BEGIN_DOC
   ! Diagonalizes the complex-scaled CI Hamiltonian H_CS and updates the wavefunction.
   !
   ! Uses either iterative Davidson or direct LAPACK diagonalization (controlled by
   ! diag_algorithm). After diagonalization, eigenvectors are c-normalized (sum v_i^2 = 1),
   ! energies and S^2 values are computed via (Psi|H_CS|Psi)/(Psi|Psi), and a first-order
   ! correction to the CS energy is evaluated from the one-body density matrix.
   !
   ! u_in  : on input, initial guess vectors; on output, converged CS eigenvectors
   ! energy : total CS energy for each state, E^{CS} + e^{-i*theta} * V_nn
   ! corr   : first-order kinetic correction, e^{-i*theta} * Tr[gamma * (e^{-i*theta}-1)*T]
   END_DOC
   complex*16, intent(inout) :: u_in(N_det,N_states)
   complex*16, intent(out) :: energy(N_states), corr(N_states)

   complex*16, allocatable :: ci_eigenvectors_cs(:,:)
   complex*16, allocatable :: ci_s2_cs(:)
   complex*16, allocatable :: ci_electronic_energy_cs(:)
   complex*16, allocatable :: h_tmp(:,:)
   
   integer, allocatable           :: index_good_state_array(:)
   logical, allocatable           :: good_state_array(:)
   integer                        :: i_other_state
   complex*16, allocatable  :: eigenvectors(:,:), eigenvalues(:), H_prime(:,:)
   integer                        :: i_state
   integer                        :: i,j,k, info
   complex*16, allocatable  :: s2_eigvalues(:)
   logical                        :: converged

   complex*16 :: e_itheta
   e_itheta = dcmplx(dcos(theta_cs), -dsin(theta_cs)) 

   PROVIDE threshold_davidson nthreads_davidson distributed_davidson

   ! DEBUG: dimension/consistency snapshot at subroutine entry, before any
   ! allocation or copy. Compares raw N_det/N_states/N_states_diag against
   ! the min(...) combinations actually used below and in the Davidson calls,
   ! to catch a possible large-N_det-specific mismatch (e.g. N_states_diag
   ! not actually being > N_states, or N_det exceeding an implicit assumption
   ! somewhere downstream) that would not show up at the smaller debug/mid
   ! scales where everything stayed comfortably inside typical bounds.
   write(*,'(A)') ' [DEBUG diagonalize_ci_cs] entry dimensions:'
   write(*,'(A,I10,A,I6,A,I6)') '   N_det=', N_det, '  N_states=', N_states, &
     '  N_states_diag=', N_states_diag
   write(*,'(A,I6,A,I6)') '   min(N_det,N_states)=', min(N_det,N_states), &
     '  min(N_det,N_states_diag)=', min(N_det,N_states_diag)
   write(*,'(A,I10,A,I10)') '   size(u_in,1)=', size(u_in,1), '  size(u_in,2)=', size(u_in,2)
   write(*,'(A,A)') '   diag_algorithm=', trim(diag_algorithm)

   allocate(ci_eigenvectors_cs(N_det,N_states_diag))
   allocate(ci_s2_cs(N_states_diag))
   allocate(ci_electronic_energy_cs(N_states_diag))

   ! Guess values for the "N_states" states of the |ci| eigenvectors
   do j=1,min(N_states,N_det)
     do i=1,N_det
       ci_eigenvectors_cs(i,j) = u_in(i,j)
     enddo
   enddo

   do j=min(N_states,N_det)+1,N_states_diag
     do i=1,N_det
       ci_eigenvectors_cs(i,j) = dcmplx(0.d0,0d0)
     enddo
   enddo

   if (diag_algorithm == "Davidson") then

     call davidson_diag_HS2_complex_cs(psi_det,ci_eigenvectors_cs, ci_s2_cs, &
       size(ci_eigenvectors_cs,1),ci_electronic_energy_cs,               &
       N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,converged)

     ! DEBUG: state of ci_eigenvectors_cs as returned by the FIRST Davidson call,
     ! before the "not converged -> double N_states_diag and retry" logic below
     ! has any chance to touch it.
     write(*,'(A)') ' [DEBUG diagonalize_ci_cs] after first davidson_diag_HS2_complex_cs call:'
     write(*,'(A,L1)') '   converged=', converged
     do i_state = 1, N_states
       write(*,'(A,I3,A,ES16.8,A,ES16.8)') '   state=', i_state,          &
         '  |ci_eigenvectors_cs(1,state)|=', cdabs(ci_eigenvectors_cs(1,i_state)), &
         '  |ci_eigenvectors_cs(N_det,state)|=', cdabs(ci_eigenvectors_cs(N_det,i_state))
     enddo

     integer :: N_states_diag_save
     N_states_diag_save = N_states_diag
     do while (.not.converged)
       complex*16, allocatable :: ci_electronic_energy_cs_tmp (:)
       complex*16, allocatable :: ci_eigenvectors_cs_tmp (:,:)
       complex*16, allocatable :: ci_s2_cs_tmp (:)

       N_states_diag  *= 2

       allocate (ci_electronic_energy_cs_tmp (N_states_diag) )
       allocate (ci_eigenvectors_cs_tmp (N_det,N_states_diag) )
       allocate (ci_s2_cs_tmp (N_states_diag) )

       ci_electronic_energy_cs_tmp(1:N_states_diag_save) = ci_electronic_energy_cs(1:N_states_diag_save)
       ci_eigenvectors_cs_tmp(1:N_det,1:N_states_diag_save) = ci_eigenvectors_cs(1:N_det,1:N_states_diag_save)
       ci_s2_cs_tmp(1:N_states_diag_save) = ci_s2_cs(1:N_states_diag_save)

       call davidson_diag_HS2_complex_cs(psi_det,ci_eigenvectors_cs_tmp, ci_s2_cs_tmp, &
         size(ci_eigenvectors_cs_tmp,1),ci_electronic_energy_cs_tmp,               &
         N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,converged)

       ! DEBUG: this whole "doubling" retry block only runs if the FIRST
       ! davidson call (Point B) returned converged=.False. -- if this print
       ! ever fires, it means convergence was NOT reached on the first try,
       ! which would be a major, previously-unobserved finding on its own.
       ! Note it also only ever re-seeds from the ORIGINAL N_states_diag_save
       ! columns of ci_eigenvectors_cs (see copy-back below), not from the
       ! larger subspace just explored -- worth revisiting if this fires.
       write(*,'(A,I3,A,L1,A,ES16.8)') ' [DEBUG diagonalize_ci_cs] retry loop: N_states_diag=', &
         N_states_diag, '  converged=', converged, '  |ci_eigenvectors_cs_tmp(1,1)|=', &
         cdabs(ci_eigenvectors_cs_tmp(1,1))

       ci_electronic_energy_cs(1:N_states_diag_save) = ci_electronic_energy_cs_tmp(1:N_states_diag_save)
       ci_eigenvectors_cs(1:N_det,1:N_states_diag_save) = ci_eigenvectors_cs_tmp(1:N_det,1:N_states_diag_save)
       ci_s2_cs(1:N_states_diag_save) = ci_s2_cs_tmp(1:N_states_diag_save)

       deallocate (ci_electronic_energy_cs_tmp)
       deallocate (ci_eigenvectors_cs_tmp)
       deallocate (ci_s2_cs_tmp)
     enddo

     if (N_states_diag > N_states_diag_save) then
       N_states_diag = N_states_diag_save
      deallocate(ci_eigenvectors_cs, ci_electronic_energy_cs, ci_s2_cs)
       allocate(ci_eigenvectors_cs(N_det,N_states_diag))
       allocate(ci_s2_cs(N_states_diag))
       allocate(ci_electronic_energy_cs(N_states_diag))
       
     endif

   else if (diag_algorithm == "Lapack") then

     print *,  'Diagonalization of H_CS using Lapack'
     allocate (eigenvectors(size(H_matrix_all_dets_cs,1),N_det))
     allocate (eigenvalues(N_det))

     if (s2_eig) then

       complex*16, parameter :: alpha = dcmplx(0.1d0,0d0)
       allocate (H_prime(N_det,N_det) )

       H_prime(1:N_det,1:N_det) = H_matrix_all_dets_cs(1:N_det,1:N_det) +  &
         alpha * dcmplx(S2_matrix_all_dets(1:N_det,1:N_det),0d0)

       do j=1,N_det
         H_prime(j,j) = H_prime(j,j) - alpha * dcmplx(expected_s2,0d0)
       enddo

       call diag_general_complex(eigenvalues,eigenvectors,H_prime,size(H_prime,1),N_det,info)
       ! Do NOT call qr_decomposition_c here: LAPACK returns Hermitian-normalised
       ! eigenvectors (||v||_2 = 1) whose c-norms (sum v_i^2) may be near-zero.
       ! qr_decomposition_c has no guard against near-zero c-norm, so it divides
       ! by ~0 and corrupts every subsequent column via cascading Gram-Schmidt.
       ! State selection below uses only dble(eigenvectors(i,k)); c-normalisation
       ! of the final eigvec is handled by the per-column loop further down.
       call nullify_small_elements_complex(N_det,N_det,eigenvectors,size(eigenvectors,1),1.d-12)

       ci_electronic_energy_cs(:) = dcmplx(0.d0,0d0)
       i_state = 0

       allocate (s2_eigvalues(N_det))
       allocate(index_good_state_array(N_det),good_state_array(N_det))

       good_state_array = .False.
       call u_0_S2_u_0_complex(s2_eigvalues,eigenvectors,N_det,psi_det,N_int,&
         N_det,size(eigenvectors,1))

       if (state_following .and. only_expected_s2) then

         integer :: state(N_states), idx,l
         double precision :: omax
         double precision, allocatable :: overlp(:)
         logical :: used
         logical, allocatable :: ok(:)

         allocate(overlp(N_det),ok(N_det))

         i_state = 0
         state = 0
         do l = 1, N_states

           ! Overlap wrt each state
           overlp = 0d0
           do k = 1, N_det
             do i = 1, N_det
               overlp(k) = overlp(k) + psi_coef(i,l) * dble(eigenvectors(i,k))
             enddo
           enddo

           ! Idx of the state with the maximum overlap not already "used"
           omax = 0d0
           idx = 0
           do k = 1, N_det

             ! Already used ?
             used = .False.
             do i = 1, N_states
               if (state(i) == k) then
                 used = .True.
               endif
             enddo

             ! Maximum overlap
             if (dabs(overlp(k)) > omax .and. .not. used) then
               if (dabs(cdabs(s2_eigvalues(k))-expected_s2) > 0.5d0) cycle
               omax = dabs(overlp(k))
               idx = k
             endif
           enddo

           state(l) = idx
           i_state +=1
         enddo

         deallocate(overlp,ok)

         do i = 1, i_state
           index_good_state_array(i) = state(i)
           good_state_array(i) = .True.
         enddo

       else if (only_expected_s2) then

         do j=1,N_det
           ! Select at least n_states states with S^2 values closed to "expected_s2"
           if(dabs(dble(s2_eigvalues(j))-expected_s2).le.0.5d0)then
             i_state +=1
             index_good_state_array(i_state) = j
             good_state_array(j) = .True.
           endif

           if(i_state.eq.N_states) then
             exit
           endif
         enddo

       else

         do j=1,N_det
           index_good_state_array(j) = j
           good_state_array(j) = .True.
         enddo

       endif

       if(i_state .ne.0)then

         ! Fill the first "i_state" states that have a correct S^2 value
         do j = 1, i_state
           do i=1,N_det
             ci_eigenvectors_cs(i,j) = eigenvectors(i,index_good_state_array(j))
           enddo
           ci_electronic_energy_cs(j) = eigenvalues(index_good_state_array(j))
           ci_s2_cs(j) = s2_eigvalues(index_good_state_array(j))
         enddo

         i_other_state = 0
         do j = 1, N_det
           if(good_state_array(j))cycle
           i_other_state +=1
           if(i_state+i_other_state.gt.n_states_diag)then
             exit
           endif
           do i=1,N_det
             ci_eigenvectors_cs(i,i_state+i_other_state) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_cs(i_state+i_other_state) = eigenvalues(j)
           ci_s2_cs(i_state+i_other_state) = s2_eigvalues(i_state+i_other_state)
         enddo

       else
         print*,''
         print*,'!!!!!!!!   WARNING  !!!!!!!!!!'
         print*,'  Within the ',N_det,'determinants selected'
         print*,'  and the ',N_states_diag,'states requested'
         print*,'  We did not find only states with S^2 values close to ',expected_s2
         print*,'  We will then set the first N_states eigenvectors of the H_CS matrix'
         print*,'  as the ci_eigenvectors_cs'
         print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the ci space'
         print*,''
         do j=1,min(N_states_diag,N_det)
           do i=1,N_det
             ci_eigenvectors_cs(i,j) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_cs(j) = eigenvalues(j)
           ci_s2_cs(j) = s2_eigvalues(j)
         enddo
       endif

       deallocate(index_good_state_array,good_state_array)
       deallocate(s2_eigvalues)

     else

       allocate(h_tmp(N_det,N_det))
       
       h_tmp = H_matrix_all_dets_cs
       call diag_general_complex(eigenvalues,eigenvectors,                    &
           h_tmp,size(h_tmp,1),N_det,info)
       deallocate(h_tmp)

       allocate(s2_eigvalues(N_det))
       call u_0_S2_u_0_complex(s2_eigvalues,eigenvectors,N_det,psi_det,N_int, &
           min(N_det,N_states_diag),size(eigenvectors,1))
       
       ! Select the "N_states_diag" states of lowest energy
       do j=1,min(N_det,N_states_diag)
         do i=1,N_det
           ci_eigenvectors_cs(i,j) = eigenvectors(i,j)
         enddo
         ci_electronic_energy_cs(j) = eigenvalues(j)
         write(*,'(I8,2(F22.10))') j, ci_electronic_energy_cs(j) + e_itheta*dcmplx(nuclear_repulsion,0d0)
      enddo
     
     endif

     do k=1,N_states_diag
       ci_electronic_energy_cs(k) = dcmplx(0.d0,0d0)
       do j=1,N_det
         do i=1,N_det
           ci_electronic_energy_cs(k) +=                                &
               ci_eigenvectors_cs(i,k) * ci_eigenvectors_cs(j,k) *         &
               H_matrix_all_dets_cs(i,j)
         enddo
       enddo
     enddo

     deallocate(eigenvectors,eigenvalues)
   endif

   do k = 1, N_states
     do i = 1, N_det
       u_in(i,k) = ci_eigenvectors_cs(i,k)
     enddo
   enddo

   ! c-normalization and calculation of the energy
   complex*16, allocatable :: eigvec(:,:)
   complex*16 :: res
   integer :: k_col
   complex*16 :: cnorm_col
   allocate(eigvec(N_det,N_states))

   eigvec = ci_eigenvectors_cs(:,1:N_states)

   ! DEBUG: eigvec right after the copy from ci_eigenvectors_cs, before any
   ! normalization. Compares directly against Point A/B: if this is already
   ! zero while Point A/B showed non-zero, the corruption happens in the
   ! "doubling" retry block or in the u_in -> ci_eigenvectors_cs -> eigvec
   ! copy chain above.
   write(*,'(A)') ' [DEBUG diagonalize_ci_cs] eigvec right after copy, pre-normalization:'
   do k_col = 1, N_states
     write(*,'(A,I3,A,ES16.8,A,F18.10,A,F18.10)') '   state=', k_col,     &
       '  |eigvec(1,state)|=', cdabs(eigvec(1,k_col)),                     &
       '  Re(E)=', dble(ci_electronic_energy_cs(k_col)),                   &
       '  Im(E)=', dimag(ci_electronic_energy_cs(k_col))
   enddo

   ! Per-column c-normalization with guard against near-zero c-norm.
   ! We do NOT use qr_decomposition_c (Gram-Schmidt) here because it mixes all
   ! states via subtraction of projections: if one state has a near-zero c-norm
   ! (e.g. Davidson tracked a wrong eigenvalue), the Gram-Schmidt blows up and
   ! corrupts every other state through cascading garbage.
   ! The energy/S2 computations below only read diagonal elements of U^T H U,
   ! so inter-state c-orthogonality is not required.
   do k_col = 1, N_states
     call inner_prod_c(eigvec(1,k_col), eigvec(1,k_col), N_det, cnorm_col)
     ! DEBUG: unconditional print of the c-norm actually seen for every state,
     ! even when the primary (c-norm) branch succeeds normally. This is the
     ! quantity being tested against the 1d-12 threshold below.
     write(*,'(A,I3,A,ES16.8,A,ES16.8)') ' [DEBUG diagonalize_ci_cs] state=', k_col, &
       '  c-norm Re=', dble(cnorm_col), '  Im=', dimag(cnorm_col)
     if (cdabs(cnorm_col) > 1d-12) then
       eigvec(:,k_col) = eigvec(:,k_col) / cdsqrt(cnorm_col)
     else
       ! Near-zero c-norm: fall back to Hermitian normalization to keep the
       ! vector finite and flag the state as unconverged via its residue.
       call inner_product_complex(eigvec(1,k_col), eigvec(1,k_col), N_det, cnorm_col)
       write(*,'(A,I3,A,ES16.8)') ' [DEBUG diagonalize_ci_cs] state=', k_col, &
         '  fallback Hermitian norm=', dble(cnorm_col)
       if (dble(cnorm_col) > 1d-20) then
         eigvec(:,k_col) = eigvec(:,k_col) / dsqrt(dble(cnorm_col))
       else
         write(*,'(A,I3,A)') ' [DEBUG diagonalize_ci_cs] state=', k_col, '  ZEROING (both norms below threshold)'
         eigvec(:,k_col) = (0d0, 0d0)
       endif
     endif
   enddo

   call overlap_cs(eigvec)
   if (overlap_analysis) then
     call overlap_cs_analysis(eigvec)
   endif

   complex*16, allocatable :: W(:,:), h(:,:), S_d(:,:), residue(:,:), H_jj(:)
   double precision, external     :: diag_H_mat_elem
   complex*16, external           :: diag_H_mat_elem_cs
   allocate(W(N_det,N_states),h(N_states,N_states),S_d(N_det,N_states),residue(N_det,N_states), H_jj(N_det))

  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(N_det, H_jj, psi_det, N_int, e_itheta)                    &
      !$OMP  PRIVATE(i)
  !$OMP DO SCHEDULE(static)
  do i=1, N_det
    H_jj(i) = e_itheta * (dcmplx(diag_H_mat_elem(psi_det(1,1,i),N_int), 0d0) + diag_H_mat_elem_cs(psi_det(1,1,i),N_int))
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


   if (diag_algorithm == "Davidson") then

     ! E = U^T H_CS U
     call H_S2_u_0_nstates_openmp_complex_cs(W,S_d,eigvec,N_states,N_det)

     call zgemm('T','N', N_states, N_states, N_det,                       &
       (1.d0,0d0), eigvec, size(eigvec,1), W, size(W,1),                          &
       (0.d0,0d0), h, size(h,1))

     do i = 1, N_states
       ci_electronic_energy_cs(i) = h(i,i)
       ! ci_electronic_energy_cs holds E^{CS} = e^{-i*theta} * E^{CAP/CS}.
       ! H_matrix_all_dets_cs is already expressed in the H^{CS} space.
     enddo

     ! Check the residue of U^T H_CS U
     do i = 1, N_states
       residue(:,i) = (h(i,i) * eigvec(:,i) - W(:,i)) / (H_jj(:) - h(i,i))
     enddo

     write(*,'(A)') ''
     write(*,'(A)') ' Residue of (Psi|H_CS|Psi)/(Psi|Psi):'
     write(*,'(A)') ' ======= ========= ========='
     write(*,'(A)') '  State   Re(Res)   Im(Res) '
     write(*,'(A)') ' ======= ========= ========='
     do i = 1, N_states
       call inner_prod_c(residue(1,i),residue(1,i), N_det, res)
       write(*,'(1X,I4,4X,ES8.1,2X,ES8.1)') i, dble(res), dimag(res)
     enddo
     write(*,'(A)') ' ======= ========= ========='

     ! S^2
     call zgemm('T','N', N_states, N_states, N_det,                       &
       (1.d0,0d0), eigvec, size(eigvec,1), S_d, size(S_d,1),                          &
       (0.d0,0d0), h, size(h,1))

     do i = 1, N_states
       ci_s2_cs(i) = h(i,i)
     enddo

   else
     call hpsi_complex(W,eigvec,N_states,N_det,H_matrix_all_dets_cs) 
     call zgemm('T','N', N_states, N_states, N_det,                       &
       (1.d0,0d0), eigvec, size(eigvec,1), W, size(W,1),                          &
       (0.d0,0d0), h, size(h,1))

     do i = 1, N_states
       ci_electronic_energy_cs(i) = h(i,i)
       ! ci_electronic_energy_cs holds E^{CS} = e^{-i*theta} * E^{CAP/CS}.
       ! H_matrix_all_dets_cs is already expressed in the H^{CS} space.
     enddo
   endif

   ! First-order correction via the 1-RDM and the kinetic energy operator T.
   ! In the H^{CS} space:
   !   corr(i) = e^{-i*theta} * Tr[ gamma_i * (e^{-i*theta}-1) * T ]
   !   first_order_e(i) = E^{CS}(i) - corr(i) ≈ e^{-i*theta} * E^H
   complex*16, allocatable :: rdm_cs(:,:,:), gw(:,:), tmp_t(:,:)
   allocate(rdm_cs(mo_num,mo_num,N_states),gw(mo_num,mo_num),tmp_t(mo_num,mo_num))
   call get_one_e_rdm_mo_cs(eigvec, rdm_cs)

   complex*16 :: first_order_e(N_states), trace
   complex*16, external :: trace_complex
   complex*16 :: f_Tprime

   ! f_T' = e^{-i*theta} - 1 = (cos(theta) - 1) - i*sin(theta)
   f_Tprime = dcmplx(dcos(theta_cs) - 1.d0, -dsin(theta_cs))
   tmp_t = f_Tprime * dcmplx(mo_kinetic_integrals, 0.d0)

   write(*,*) ''
   do i = 1, N_states
     gw = matmul(rdm_cs(1:mo_num,1:mo_num,i),tmp_t)
     trace = trace_complex(gw,mo_num)
     ! corr(i) = e^{-i*theta} * Tr[ gamma_i * f_T' * T ]
     corr(i) = e_itheta * trace
     ! first_order_e(i) ≈ e^{-i*theta} * E^H :
     ! E^{CS} - corr = e^{-i*theta}*(E^{CAP/CS} - trace) = e^{-i*theta} * E^H
     first_order_e(i) = ci_electronic_energy_cs(i) - corr(i)
     write(*,'(A,I8,2F16.10)') 'First order correction:', i, dble(corr(i)), dimag(corr(i))
   enddo
   write(*,*) ''

   ! energy() = E^{CS}_tot = E^{CS}_el + e^{-i*theta} * V_nn
   ! (ci_electronic_energy_cs already holds E^{CS}_el)
   do i = 1, N_states
     energy(i) = ci_electronic_energy_cs(i) + e_itheta * dcmplx(nuclear_repulsion,0d0)
   enddo

   write(6,*) ''
   write(6,'(A35,ES12.3)') ' Energy of the states for theta = ', theta_cs
   write(6,'(A65)') '=======  ================ ================  ========== =========='
   write(6,'(A65)') ' State      Energy (Re)      Energy (Im)     S^2 (Re)   S^2 (Im) '
   write(6,'(A65)') '=======  ================ ================  ========== =========='
   do k = 1, N_states
     ! E^{CS}_tot = E^{CS}_el + e^{-i*theta} * V_nn
     write(6,'(1X,I4,3X,F16.10,1X,F16.10,2X,ES10.2,1X,ES10.2)') k, &
       ci_electronic_energy_cs(k) + e_itheta * dcmplx(nuclear_repulsion,0d0), ci_s2_cs(k)
   enddo
   write(6,'(A65)') '=======  ================ ================  ========== =========='
   write(6,*) ''

   write(6,'(A42)') '=======  ================ ================'
   write(6,'(A42)') ' State        U (Re)           U (Im)     '
   write(6,'(A42)') '=======  ================ ================'
   do k = 1, N_states
     ! first_order_e = E^{CS} - corr ≈ e^{-i*theta} * E^H
     write(6,'(1X,I4,3X,F16.10,1X,F16.10)') k, &
       first_order_e(k) + e_itheta * dcmplx(nuclear_repulsion,0d0)
   enddo
   write(6,'(A42)') '=======  ================ ================'
   write(6,*) ''

   call ref_idx_complex_cs(psi_det,eigvec,N_det,N_states,N_int)

   ! psi_energy_cs stores E^{CS} for consistency with Hii (both in H^{CS} space) in the stochastic PT2
   ! (delta_E = psi_energy_cs - Hii, both quantities in the H^{CS} space)
   psi_energy_cs = ci_electronic_energy_cs
   touch psi_energy_cs
   psi_cs_coef = eigvec
   touch psi_cs_coef
   if (cs_write_wf) then
     call write_wf_cs()
   endif

   deallocate(W,eigvec,h,S_d,rdm_cs,tmp_t,residue)

end
