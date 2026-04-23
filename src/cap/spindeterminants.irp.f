subroutine get_psi_bilinear_matrix_values_cap(psi_coef_cap,psi_bilinear_matrix_values_cap)

  use bitmasks
  implicit none
  BEGIN_DOC
  ! Sparse coefficient matrix if the wave function is expressed in a bilinear form :
  !  $D_\alpha^\dagger.C.D_\beta$
  !
  ! Rows are $\alpha$ determinants and columns are $\beta$.
  !
  ! Order refers to psi_det
  END_DOC
  complex*16, intent(in)  :: psi_coef_cap(N_det,N_states)
  complex*16, intent(out) :: psi_bilinear_matrix_values_cap(N_det,N_states)   
  integer :: i,k

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
  do k = 1, N_states
    !$OMP DO
    do i = 1, N_det
      psi_bilinear_matrix_values_cap(i,k) = psi_coef_cap(psi_bilinear_matrix_order(i),k)
    enddo
    !$OMP ENDDO NOWAIT
  enddo
  !$OMP END PARALLEL

end

subroutine get_psi_bilinear_matrix_transp_values_cap(psi_coef_cap,psi_bilinear_matrix_transp_values_cap)

  use bitmasks
  implicit none
  BEGIN_DOC
  ! Transpose of :c:data:`psi_bilinear_matrix`
  !
  ! $D_\beta^\dagger.C^\dagger.D_\alpha$
  !
  ! Rows are $\alpha$ determinants and columns are $\beta$, but the matrix is stored in row major
  ! format.
  END_DOC
  complex*16, intent(in)  :: psi_coef_cap(N_det,N_states)
  complex*16, intent(out) :: psi_bilinear_matrix_transp_values_cap(N_det,N_states)
  complex*16, allocatable :: psi_bilinear_matrix_values_cap(:,:)
  integer                        :: i,k

  allocate(psi_bilinear_matrix_values_cap(N_det,N_states))
  call get_psi_bilinear_matrix_values_cap(psi_coef_cap,psi_bilinear_matrix_values_cap)

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k)
  do k=1,N_states
    !$OMP DO
    do i=1,N_det
      psi_bilinear_matrix_transp_values_cap(i,k) = psi_bilinear_matrix_values_cap(psi_bilinear_matrix_transp_order(i),k)
    enddo
    !$OMP ENDDO NOWAIT
  enddo
  !$OMP END PARALLEL

  deallocate(psi_bilinear_matrix_values_cap)
end
