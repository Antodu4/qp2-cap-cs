
BEGIN_PROVIDER [ complex*16, psi_energy_cs, (N_states) ]
  implicit none
  BEGIN_DOC
    ! psi_energy_cs(i) = (Psi_i | H_CS(theta) | Psi_i)
  END_DOC
  integer :: i
  do i=1,N_states
    psi_energy_cs(i) = (0.d0,0.d0)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, psi_s2_cs, (N_states) ]
  implicit none
  BEGIN_DOC
    ! psi_s2_cs(i) = (Psi_i | S^2 | Psi_i)
  END_DOC
  integer :: i
  do i=1,N_states
    psi_s2_cs(i) = (0.d0,0.d0)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, psi_cs_coef, (N_det,N_states) ]
   implicit none
   BEGIN_DOC
   ! CS wave function coefficients
   END_DOC

   psi_cs_coef = (0d0,0d0)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pt2_im_match_coef, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Coefficient to match real and imaginary part of the pt2 in the selection process
   END_DOC

   pt2_im_match_coef = 1d0

END_PROVIDER
