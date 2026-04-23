subroutine i_W_j_single_spin_cap(key_i,key_j,Nint,spin,wij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|W|j \rangle$ where $i$ and $j$ are determinants differing by
  ! a single excitation.
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: wij

  integer                        :: exc(0:2,2)
  double precision               :: phase

  call get_single_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  if (spin == 1) then
     wij = mo_ints_cap_alpha(exc(1,1),exc(1,2)) * phase
  else
     wij = mo_ints_cap_beta(exc(1,1),exc(1,2)) * phase
  endif
  wij = - eta_cap * wij
end

double precision function diag_H_mat_elem_cap(det_in,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|W|i \rangle$.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)

  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin,iorb
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb

  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)

  integer                        :: tmp(2)
  tmp(1) = elec_alpha_num
  tmp(2) = elec_beta_num
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ_particle, tmp, Nint)

  diag_H_mat_elem_cap = 0d0

  ! alpha part
  do i=1,sum(popcnt(det_in(:,1)))
    iorb = occ_particle(i,1)
    diag_H_mat_elem_cap += mo_ints_cap_alpha(iorb,iorb)
  enddo

  ! beta part
  do i=1,sum(popcnt(det_in(:,2)))
    iorb = occ_particle(i,2)
    diag_H_mat_elem_cap += mo_ints_cap_beta(iorb,iorb)
  enddo
  
  diag_H_mat_elem_cap = - eta_cap * diag_H_mat_elem_cap
end

BEGIN_PROVIDER [ complex*16, H_matrix_all_dets_complex,(N_det,N_det) ]
  use bitmasks
 implicit none
 BEGIN_DOC
 ! |H| matrix on the basis of the Slater determinants defined by psi_det
 END_DOC
 integer :: i,j,k,h1,p1,h2,p2,s1,s2
 double precision :: hij, wij, phase
 double precision, external :: diag_H_mat_elem_cap
 integer :: degree, exc(0:2,2,2)
 call  i_H_j(psi_det(1,1,1),psi_det(1,1,1),N_int,hij)
 !$OMP PARALLEL DO SCHEDULE(GUIDED) DEFAULT(NONE) PRIVATE(i,j,wij,hij,degree,k,&
 !$OMP phase, h1,p1,h2,p2,s1,s2,exc) &
 !$OMP SHARED (N_det, psi_det, N_int,H_matrix_all_dets_complex)
 do i =1,N_det
   do j = i, N_det
    call  i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)

    if (degree == 0) then
      wij = diag_H_mat_elem_cap(psi_det(1,1,i),N_int)
    else if (degree == 1) then
      call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
      call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
      call i_W_j_single_spin_cap(psi_det(1,1,i),psi_det(1,1,j),N_int,s1,wij)
    else
      wij = 0d0
    endif

    H_matrix_all_dets_complex(i,j) = dcmplx(hij,wij)
    H_matrix_all_dets_complex(j,i) = dcmplx(hij,wij)
  enddo
 enddo
 !$OMP END PARALLEL DO
END_PROVIDER
  
