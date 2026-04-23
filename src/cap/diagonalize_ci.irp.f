
subroutine diagonalize_ci_cap(u_in, energy, corr)

   implicit none
   BEGIN_DOC
   ! Diagonalizes the complex CAP Hamiltonian H + i*eta*W in the CI basis
   ! and computes the complex energies and first-order CAP corrections.
   !
   ! Uses either Davidson (iterative, preferred for large CI spaces) or
   ! Lapack (full diagonalization) according to diag_algorithm.
   ! After diagonalization, eigenvectors are c-normalized (biorthogonal
   ! convention: (u,v) = sum_i u_i*v_i without conjugation) via QR decomposition.
   !
   ! The first-order CAP correction is computed as:
   !   corr(i) = -i * eta * Tr[gamma * W_mo]
   ! where gamma is the one-body density matrix and W_mo the CAP in MO basis.
   !
   ! Inputs/outputs:
   !   u_in(N_det,N_states)   : on input, initial guess vectors; on output,
   !                            the c-normalized CAP eigenvectors
   ! Outputs:
   !   energy(N_states)       : complex total energies (electronic + nuclear repulsion)
   !   corr(N_states)         : first-order energy corrections
   END_DOC
   complex*16, intent(inout) :: u_in(N_det,N_states)
   complex*16, intent(out) :: energy(N_states), corr(N_states)

   complex*16, allocatable :: ci_eigenvectors_cap(:,:)
   complex*16, allocatable :: ci_s2_cap(:)
   complex*16, allocatable :: ci_electronic_energy_cap(:)
   complex*16, allocatable :: h_tmp(:,:)
   
   integer                        :: i_good_state
   integer, allocatable           :: index_good_state_array(:)
   logical, allocatable           :: good_state_array(:)
   complex*16, allocatable  :: s2_values_tmp(:)
   integer                        :: i_other_state
   complex*16, allocatable  :: eigenvectors(:,:), eigenvalues(:), H_prime(:,:)
   integer                        :: i_state
   complex*16               :: e_0
   integer                        :: i,j,k, info
   complex*16, allocatable  :: s2_eigvalues(:)
   integer, allocatable           :: iorder(:)
   logical                        :: converged
   logical                        :: do_csf

   PROVIDE threshold_davidson nthreads_davidson distributed_davidson

   allocate(ci_eigenvectors_cap(N_det,N_states_diag))
   allocate(ci_s2_cap(N_states_diag))
   allocate(ci_electronic_energy_cap(N_states_diag))

   ! Guess values for the "N_states" states of the |ci| eigenvectors
   do j=1,min(N_states,N_det)
     do i=1,N_det
       ci_eigenvectors_cap(i,j) = u_in(i,j)
     enddo
   enddo

   do j=min(N_states,N_det)+1,N_states_diag
     do i=1,N_det
       ci_eigenvectors_cap(i,j) = dcmplx(0.d0,0d0)
     enddo
   enddo

   if (diag_algorithm == "Davidson") then

     call davidson_diag_HS2_complex(psi_det,ci_eigenvectors_cap, ci_s2_cap, &
       size(ci_eigenvectors_cap,1),ci_electronic_energy_cap,               &
       N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,converged)

     integer :: N_states_diag_save
     N_states_diag_save = N_states_diag
     do while (.not.converged)
       complex*16, allocatable :: ci_electronic_energy_cap_tmp (:)
       complex*16, allocatable :: ci_eigenvectors_cap_tmp (:,:)
       complex*16, allocatable :: ci_s2_cap_tmp (:)

       N_states_diag  *= 2

       allocate (ci_electronic_energy_cap_tmp (N_states_diag) )
       allocate (ci_eigenvectors_cap_tmp (N_det,N_states_diag) )
       allocate (ci_s2_cap_tmp (N_states_diag) )

       ci_electronic_energy_cap_tmp(1:N_states_diag_save) = ci_electronic_energy_cap(1:N_states_diag_save)
       ci_eigenvectors_cap_tmp(1:N_det,1:N_states_diag_save) = ci_eigenvectors_cap(1:N_det,1:N_states_diag_save)
       ci_s2_cap_tmp(1:N_states_diag_save) = ci_s2_cap(1:N_states_diag_save)

       call davidson_diag_HS2_complex(psi_det,ci_eigenvectors_cap_tmp, ci_s2_cap_tmp, &
         size(ci_eigenvectors_cap_tmp,1),ci_electronic_energy_cap_tmp,               &
         N_det,min(N_det,N_states),min(N_det,N_states_diag),N_int,converged)

       ci_electronic_energy_cap(1:N_states_diag_save) = ci_electronic_energy_cap_tmp(1:N_states_diag_save)
       ci_eigenvectors_cap(1:N_det,1:N_states_diag_save) = ci_eigenvectors_cap_tmp(1:N_det,1:N_states_diag_save)
       ci_s2_cap(1:N_states_diag_save) = ci_s2_cap_tmp(1:N_states_diag_save)

       deallocate (ci_electronic_energy_cap_tmp)
       deallocate (ci_eigenvectors_cap_tmp)
       deallocate (ci_s2_cap_tmp)
     enddo

     if (N_states_diag > N_states_diag_save) then
       N_states_diag = N_states_diag_save
       deallocate(ci_eigenvectors_cap, ci_electronic_energy_cap, ci_s2_cap)
       allocate(ci_eigenvectors_cap(N_det,N_states_diag))
       allocate(ci_s2_cap(N_states_diag))
       allocate(ci_electronic_energy_cap(N_states_diag))
     endif

   else if (diag_algorithm == "Lapack") then

     print *,  'Diagonalization of H using Lapack'
     allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
     allocate (eigenvalues(N_det))

     if (s2_eig) then

       complex*16, parameter :: alpha = dcmplx(0.1d0,0d0)
       allocate (H_prime(N_det,N_det) )

       H_prime(1:N_det,1:N_det) = H_matrix_all_dets_complex(1:N_det,1:N_det) +  &
         alpha * dcmplx(S2_matrix_all_dets(1:N_det,1:N_det),0d0)

       do j=1,N_det
         H_prime(j,j) = H_prime(j,j) - alpha * dcmplx(expected_s2,0d0)
       enddo

       call diag_general_complex(eigenvalues,eigenvectors,H_prime,size(H_prime,1),N_det,info)
       call qr_decomposition_c(eigenvectors,size(eigenvectors,1),N_det,N_det)
       call nullify_small_elements_complex(N_det,N_det,eigenvectors,size(eigenvectors,1),1.d-12)

       ci_electronic_energy_cap(:) = dcmplx(0.d0,0d0)
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
             ci_eigenvectors_cap(i,j) = eigenvectors(i,index_good_state_array(j))
           enddo
           ci_electronic_energy_cap(j) = eigenvalues(index_good_state_array(j))
           ci_s2_cap(j) = s2_eigvalues(index_good_state_array(j))
         enddo

         i_other_state = 0
         do j = 1, N_det
           if(good_state_array(j))cycle
           i_other_state +=1
           if(i_state+i_other_state.gt.n_states_diag)then
             exit
           endif
           do i=1,N_det
             ci_eigenvectors_cap(i,i_state+i_other_state) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_cap(i_state+i_other_state) = eigenvalues(j)
           ci_s2_cap(i_state+i_other_state) = s2_eigvalues(i_state+i_other_state)
         enddo

       else
         print*,''
         print*,'!!!!!!!!   WARNING  !!!!!!!!!'
         print*,'  Within the ',N_det,'determinants selected'
         print*,'  and the ',N_states_diag,'states requested'
         print*,'  We did not find only states with S^2 values close to ',expected_s2
         print*,'  We will then set the first N_states eigenvectors of the H matrix'
         print*,'  as the ci_eigenvectors_cap'
         print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the ci space'
         print*,''
         do j=1,min(N_states_diag,N_det)
           do i=1,N_det
             ci_eigenvectors_cap(i,j) = eigenvectors(i,j)
           enddo
           ci_electronic_energy_cap(j) = eigenvalues(j)
           ci_s2_cap(j) = s2_eigvalues(j)
         enddo
       endif

       deallocate(index_good_state_array,good_state_array)
       deallocate(s2_eigvalues)

     else

       allocate(h_tmp(N_det,N_det))

       h_tmp = H_matrix_all_dets_complex
       call diag_general_complex(eigenvalues,eigenvectors,                    &
           h_tmp,size(h_tmp,1),N_det,info)

       deallocate(h_tmp)

       ci_electronic_energy_cap(:) = dcmplx(0.d0,0d0)
       call u_0_S2_u_0_complex(ci_s2_cap,eigenvectors,N_det,psi_det,N_int,       &
           min(N_det,N_states_diag),size(eigenvectors,1))

       ! Select the "N_states_diag" states of lowest energy
       do j=1,min(N_det,N_states_diag)
         do i=1,N_det
           ci_eigenvectors_cap(i,j) = eigenvectors(i,j)
         enddo
         ci_electronic_energy_cap(j) = eigenvalues(j)
         write(*,'(I8,2(F22.10))') j, ci_electronic_energy_cap(j) + dcmplx(nuclear_repulsion,0d0)
       enddo

     endif

     do k=1,N_states_diag
       ci_electronic_energy_cap(k) = dcmplx(0.d0,0d0)
       do j=1,N_det
         do i=1,N_det
           ci_electronic_energy_cap(k) +=                                &
               ci_eigenvectors_cap(i,k) * ci_eigenvectors_cap(j,k) *         &
               H_matrix_all_dets_complex(i,j)
         enddo
       enddo
     enddo

     deallocate(eigenvectors,eigenvalues)
   endif

   do k = 1, N_states
     do i = 1, N_det
       u_in(i,k) = ci_eigenvectors_cap(i,k)
     enddo
   enddo

   ! c-normalization and calculation of the energy
   complex*16, allocatable :: eigvec(:,:), eige(:)
   complex*16 :: res
   allocate(eigvec(N_det,N_states),eige(N_states))

   eigvec = ci_eigenvectors_cap(:,1:N_states)

   call qr_decomposition_c(eigvec,size(eigvec,1),N_det,N_states)

   call overlap_cap(eigvec)
   if (overlap_analysis) then
     call overlap_cap_analysis(eigvec)
   endif

   complex*16, allocatable :: W(:,:), h(:,:), S_d(:,:), prefactor(:), residue(:,:), H_jj(:)
   double precision, external     :: diag_H_mat_elem
   allocate(W(N_det,N_states),h(N_states,N_states),S_d(N_det,N_states),prefactor(N_states), residue(N_det,N_states), H_jj(N_det))

  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP  SHARED(N_det,H_jj, psi_det,N_int)                    &
      !$OMP  PRIVATE(i)
  !$OMP DO SCHEDULE(static)
  do i=1, N_det
    H_jj(i)  = dcmplx(diag_H_mat_elem(psi_det(1,1,i),N_int), 0d0)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL


   if (diag_algorithm == "Davidson") then

     ! E = U^T H U 
     call H_S2_u_0_nstates_openmp_complex(W,S_d,eigvec,N_states,N_det)

     call zgemm('T','N', N_states, N_states, N_det,                       &
       (1.d0,0d0), eigvec, size(eigvec,1), W, size(W,1),                          &
       (0.d0,0d0), h, size(h,1))

     do i = 1, N_states
       ci_electronic_energy_cap(i) = h(i,i)
     enddo     

     ! Check the residue of U^T H U since we are using U^* H U in Davidson...     
     do i = 1, N_states
       residue(:,i) = (h(i,i) * eigvec(:,i) - W(:,i)) / (H_jj(:) - h(i,i))
     enddo

     write(*,'(A)') ''
     write(*,'(A)') ' Residue of (Psi|H|Psi)/(Psi|Psi):'
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
       ci_s2_cap(i) = h(i,i)
     enddo

   else
     call hpsi_complex(W,eigvec,N_states,N_det,H_matrix_all_dets_complex) 
     call zgemm('T','N', N_states, N_states, N_det,                       &
       (1.d0,0d0), eigvec, size(eigvec,1), W, size(W,1),                          &
       (0.d0,0d0), h, size(h,1))

     do i = 1, N_states
       ci_electronic_energy_cap(i) = h(i,i)
     enddo
   endif

   complex*16, allocatable :: rdm_cap(:,:,:), gw(:,:), tmp_w(:,:)
   allocate(rdm_cap(mo_num,mo_num,N_states),gw(mo_num,mo_num),tmp_w(mo_num,mo_num))
   call get_one_e_rdm_mo_cap(eigvec, rdm_cap)

   complex*16 :: first_order_e(N_states), trace
   complex*16, external :: trace_complex
   double precision :: corr_re, corr_im

   tmp_w = dcmplx(mo_one_e_integrals_cap,0d0)

   write(*,*) ''
   do i = 1, N_states
     gw = matmul(rdm_cap(1:mo_num,1:mo_num,i),tmp_w)
     trace = eta_cap * trace_complex(gw,mo_num)
     corr_re = dimag(trace)
     corr_im = - dble(trace)
     corr(i) = - dcmplx(corr_re,corr_im)
     first_order_e(i) = ci_electronic_energy_cap(i) - dcmplx(corr_re,corr_im)
     write(*,'(A,I8,F16.10,F16.10)') 'First order correction:', i, dcmplx(corr_re,corr_im)
   enddo
   write(*,*) ''

   do i = 1, N_states
     energy(i) = ci_electronic_energy_cap(i) + dcmplx(nuclear_repulsion,0d0)
   enddo
   
   write(6,*) ''
   write(6,'(A32,ES12.3)') ' Energy of the states for eta = ', eta_cap
   write(6,'(A65)') '=======  ================ ================  ========== =========='
   write(6,'(A65)') ' State      Energy (Re)      Energy (Im)     S^2 (Re)   S^2 (Im) '
   write(6,'(A65)') '=======  ================ ================  ========== =========='
   do k = 1, N_states
     write(6,'(1X,I4,3X,F16.10,1X,F16.10,2X,ES10.2,1X,ES10.2)') k, ci_electronic_energy_cap(k) + dcmplx(nuclear_repulsion,0d0), ci_s2_cap(k)
   enddo
   write(6,'(A65)') '=======  ================ ================  ========== =========='
   write(6,*) ''

   write(6,'(A42)') '=======  ================ ================'
   write(6,'(A42)') ' State        U (Re)           U (Im)     '
   write(6,'(A42)') '=======  ================ ================'
   do k = 1, N_states
     write(6,'(1X,I4,3X,F16.10,1X,F16.10)') k, first_order_e(k) + dcmplx(nuclear_repulsion,0d0)
   enddo
   write(6,'(A42)') '=======  ================ ================'
   write(6,*) ''

   ! We also need to know the dominant determinants for the states since they are
   ! swapping because of the CAP...
   call ref_idx_complex(psi_det,eigvec,N_det,N_states,N_int) 

   psi_energy_cap = ci_electronic_energy_cap
   touch psi_energy_cap
   psi_cap_coef = eigvec
   touch psi_cap_coef
   if (cap_write_wf) then
     call write_wf()
   endif

   deallocate(W,eigvec,eige,h,S_d,rdm_cap,tmp_w,residue)

end
