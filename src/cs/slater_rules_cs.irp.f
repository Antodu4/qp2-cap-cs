subroutine i_T_j_single_spin_cs(key_i,key_j,Nint,spin,tij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the CS contribution $f_T \langle i|T|j \rangle$ where $i$ and $j$
  ! are determinants differing by a single excitation.
  !
  ! The CS Hamiltonian is built in two steps:
  !   1) Diagonalize H^{CAP/CS} = H - (1 - e^{-i*theta}) * T  (this routine)
  !   2) Multiply eigenvalues by e^{-i*theta} after diagonalization
  !
  ! So the kinetic prefactor applied here is:
  !   f_T = -(1 - e^{-i*theta}) = e^{-i*theta} - 1
  !       = (cos(theta) - 1) - i*sin(theta)
  !
  ! T = mo_kinetic_integrals (real, one-body operator).
  ! tij is complex*16.
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  complex*16, intent(out)        :: tij

  integer                        :: exc(0:2,2)
  double precision               :: phase
  double precision               :: t_mo
  complex*16                     :: f_T

  ! f_T = e^{-i*theta} - 1 = (cos(theta) - 1) - i*sin(theta)
  ! (convention: e_itheta = e^{-i*theta} everywhere in the code; Im(E_resonance) < 0)
  f_T = dcmplx(dcos(theta_cs) - 1.d0, -dsin(theta_cs))

  call get_single_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  t_mo = mo_kinetic_integrals(exc(1,1),exc(1,2)) * phase

  tij = f_T * dcmplx(t_mo, 0.d0)
end


complex*16 function diag_H_mat_elem_cs(det_in,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes the CAP/CS diagonal kinetic contribution $f_T \langle i|T|i \rangle$.
  !
  ! This is the diagonal element of H^{CAP/CS} = H - (1 - e^{-i*theta}) * T.
  ! The global e^{-i*theta} factor is applied AFTER diagonalization.
  !
  ! f_T = -(1 - e^{-i*theta}) = e^{-i*theta} - 1 = (cos(theta) - 1) - i*sin(theta)
  !
  ! T is the one-body kinetic operator (real, diagonal in MO basis).
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det_in(Nint,2)

  integer                        :: i, iorb
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: t_diag
  complex*16                     :: f_T

  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)

  integer :: tmp(2)
  tmp(1) = elec_alpha_num
  tmp(2) = elec_beta_num
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(det_in, occ, tmp, Nint)

  ! f_T = e^{-i*theta} - 1 = (cos(theta) - 1) - i*sin(theta)
  ! (convention: e_itheta = e^{-i*theta} everywhere in the code; Im(E_resonance) < 0)
  f_T = dcmplx(dcos(theta_cs) - 1.d0, -dsin(theta_cs))

  t_diag = 0.d0

  ! alpha part
  do i = 1, sum(popcnt(det_in(:,1)))
    iorb = occ(i,1)
    t_diag = t_diag + mo_kinetic_integrals(iorb,iorb)
  enddo

  ! beta part
  do i = 1, sum(popcnt(det_in(:,2)))
    iorb = occ(i,2)
    t_diag = t_diag + mo_kinetic_integrals(iorb,iorb)
  enddo

  diag_H_mat_elem_cs = f_T * dcmplx(t_diag, 0.d0)
end
