double precision function pt2_stop_value(pt2_data)
  use selection_types
  implicit none
  BEGIN_DOC
  ! Value compared with pt2_max in the stopping test of the CIPSI loops of the
  ! complex methods (CAP and CS).
  !
  ! pt2_max_complex = .False. : largest |Re(PT2)| over the states, as in the
  !                             real CIPSI (former behaviour).
  ! pt2_max_complex = .True.  : largest modulus of the complex PT2 over the
  !                             states, sqrt(Re(PT2)**2 + Im(PT2)**2).
  !
  ! With a complex Hamiltonian, the real part of the PT2 changes sign from one
  ! CIPSI iteration to the next and can be close to zero by chance while the
  ! imaginary part is still large. The modulus is small only when both parts
  ! are small, so the loop no longer stops on a space that is not converged
  ! for Im(E).
  !
  ! At the start of a loop, pt2 and pt2_im are set to -huge: both values are
  ! then huge and the loop is entered.
  END_DOC
  type(pt2_type), intent(in)     :: pt2_data
  integer                        :: k

  pt2_stop_value = 0.d0
  if (pt2_max_complex) then
    do k = 1, N_states
      pt2_stop_value = max(pt2_stop_value,                           &
          dsqrt(pt2_data % pt2(k)**2 + pt2_data % pt2_im(k)**2))
    enddo
  else
    do k = 1, N_states
      pt2_stop_value = max(pt2_stop_value, dabs(pt2_data % pt2(k)))
    enddo
  endif

end
