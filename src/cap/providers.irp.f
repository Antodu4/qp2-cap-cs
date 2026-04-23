 BEGIN_PROVIDER [ complex*16, psi_energy_cap, (N_states) ]
&BEGIN_PROVIDER [ complex*16, psi_s2_cap, (N_states) ]
  implicit none
  BEGIN_DOC
! psi_energy_cap(i) = $(\Psi_i | H(eta) | \Psi_i)$
!    
! psi_s2_cap(i) = $(\Psi_i | S^2 | \Psi_i)$
  END_DOC 
  integer :: i
  do i=1,N_states
    psi_energy_cap(i) = (0.d0,0d0)
    psi_s2_cap(i) = (0.d0,0d0)
  enddo
END_PROVIDER

BEGIN_PROVIDER [ complex*16, psi_cap_coef, (N_det,N_states) ]
   implicit none
   BEGIN_DOC
   ! Cap wave function
   END_DOC

   psi_cap_coef = (0d0,0d0)

END_PROVIDER

BEGIN_PROVIDER [ double precision, pt2_im_match_coef, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Coefficient to match real and imaginary part of the pt2 in the selection process
   END_DOC

   pt2_im_match_coef = 1d0

END_PROVIDER
