BEGIN_PROVIDER [ complex*16, H_matrix_all_dets_cs,(N_det,N_det) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Complex-scaled Hamiltonian matrix in the determinant basis.
  ! H_cs(i,j) = e^{-i*theta} * ( <i|H|j> + (e^{-i*theta}-1) * <i|T|j> )
  ! where T is the kinetic energy operator and theta is the CS rotation angle.
  END_DOC
  integer :: i,j,k,h1,p1,h2,p2,s1,s2
  double precision :: hij, phase
  complex*16       :: tij, e_itheta
  complex*16, external :: diag_H_mat_elem_cs
  integer :: degree, exc(0:2,2,2)

  ! Global phase factor
  e_itheta = dcmplx(cos(theta_cs), -sin(theta_cs))

  !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,tij,hij,degree,k,phase,h1,p1,h2,p2,s1,s2,exc) &
  !$OMP SHARED (N_det, psi_det, N_int, H_matrix_all_dets_cs, e_itheta)
  do i = 1, N_det
    do j = i, N_det
      call i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
      call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)

      if (degree == 0) then
        tij = diag_H_mat_elem_cs(psi_det(1,1,i),N_int)
      else if (degree == 1) then
        call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
        call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
        call i_T_j_single_spin_cs(psi_det(1,1,i),psi_det(1,1,j),N_int,s1,tij)
      else
        tij = (0.d0, 0.d0)
      endif

      H_matrix_all_dets_cs(i,j) = e_itheta * (dcmplx(hij, 0.0d0) + tij)
      H_matrix_all_dets_cs(j,i) = e_itheta * (dcmplx(hij, 0.0d0) + tij)
    enddo
  enddo
  !$OMP END PARALLEL DO
END_PROVIDER
