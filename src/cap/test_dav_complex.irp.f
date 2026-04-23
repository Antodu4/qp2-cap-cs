program test_dav_complex
  implicit none
  BEGIN_DOC
  ! To test the complex davidson routine
  END_DOC
  read_wf = .True.
  touch read_wf
  PROVIDE threshold_davidson nthreads_davidson
  call dav_complex
end

subroutine dav_complex
 implicit none
 BEGIN_DOC
 ! Runs the CAP-CI diagonalization using the complex Davidson algorithm.
 !
 ! First computes exact eigenvalues via full LAPACK diagonalization (for reference),
 ! then runs the iterative complex Davidson solver followed by diagonalize_ci_cap.
 END_DOC
 complex*16, allocatable :: u_in(:,:), H_jj(:), energies(:), s2_out(:), corr(:)
 integer :: sze,N_st,N_st_diag_in
 logical :: converged
 integer :: i,j,info
 double precision, allocatable :: u_in2(:,:), s2_out2(:), energies2(:)

 N_st = N_states
 N_st_diag_in = N_states_diag
 sze = N_det
 PROVIDE nuclear_repulsion

 !!! MARK THAT u_in mut dimensioned with "N_st_diag_in" as a second dimension
 allocate(u_in2(sze,N_st_diag_in),energies2(N_st_diag_in),s2_out2(N_st_diag_in))
 allocate(u_in(sze,N_st_diag_in),H_jj(sze),energies(N_st_diag_in),s2_out(N_st_diag_in),corr(N_st_diag_in))

 u_in = 0.d0
 do i = 1, N_st
  u_in(:,i) = dcmplx(psi_coef(:,i),0d0)
 enddo

 ! Lapack
 complex*16, allocatable :: eigvalues(:), eigvectors(:,:), h(:,:),h_cp(:,:)
 double precision, allocatable :: h_im(:,:), tmp(:,:)
 allocate(eigvalues(sze),eigvectors(sze,sze))

 call diag_general_complex(eigvalues,eigvectors,H_matrix_all_dets_complex,sze,sze,info)
 eigvalues += dcmplx(nuclear_repulsion,0d0)
 print*,'Exact energies:'
 do i = 1, min(10,N_det)
    write(*,'(I6,100(F18.10))') i, eigvalues(i)
 enddo

 u_in = (0d0,0d0)
 do i = 1, N_st
  u_in(1:sze,1:N_st) = dcmplx(psi_coef(1:sze,1:N_st),0d0)
 enddo
 call davidson_diag_hs2_complex(psi_det,u_in,s2_out,sze,energies,sze,N_st,N_st_diag_in,N_int,converged)
 call  diagonalize_ci_cap(u_in, energies,corr)

end
