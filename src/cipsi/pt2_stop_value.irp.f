double precision function pt2_stop_value(pt2_data)
  use selection_types
  implicit none
  BEGIN_DOC
  ! Value compared with pt2_max in the stopping test of the CIPSI loops of the
  ! complex methods (CAP and CS): largest modulus of the complex PT2 over the
  ! states, sqrt(Re(PT2)**2 + Im(PT2)**2).
  !
  ! With a complex Hamiltonian, the real part of the PT2 changes sign from one
  ! CIPSI iteration to the next and can be close to zero by chance while the
  ! imaginary part is still large. The modulus is small only when both parts
  ! are small, so the loop does not stop on a space that is not converged
  ! for Im(E).
  !
  ! At the start of a loop, pt2 and pt2_im are set to -huge: the value is
  ! then huge and the loop is entered.
  END_DOC
  type(pt2_type), intent(in)     :: pt2_data
  integer                        :: k

  pt2_stop_value = 0.d0
  do k = 1, N_states
    pt2_stop_value = max(pt2_stop_value,                             &
        dsqrt(pt2_data % pt2(k)**2 + pt2_data % pt2_im(k)**2))
  enddo

end

logical function cs_cipsi_continue(pt2_data)
  use selection_types
  implicit none
  BEGIN_DOC
  ! PT2 and variance part of the stopping test of the CS-CIPSI (both loop
  ! orders): true while the CIPSI must go on. The loops also test N_det_max
  ! and correlation_energy_ratio_max.
  !
  ! cs_stop_abs_sum = .True. (default), with the complex PT2 (cs_pt2) : the
  !   CIPSI goes on while the largest S = sum_alpha |e_alpha| over the states
  !   (pt2_data % pt2_abs, accumulated by selection_cs) is larger than
  !   cs_abs_sum_max; pt2_max and variance_max are ignored. At the start of a
  !   loop, pt2 is set to -huge and pt2_abs is not reset: the test is then true.
  ! Otherwise : the former test, modulus of the complex PT2 larger than
  !   pt2_max (pt2_stop_value) and largest variance larger than variance_max.
  END_DOC
  type(pt2_type), intent(in)     :: pt2_data
  double precision, external     :: pt2_stop_value
  double precision               :: s_max
  integer                        :: k

  if (cs_stop_abs_sum .and. cs_pt2) then
    s_max = 0.d0
    do k = 1, N_states
      if (pt2_data % pt2(k) < -1.d30) then     ! start of a loop: no PT2 computed yet
        cs_cipsi_continue = .True.
        return
      endif
      s_max = max(s_max, pt2_data % pt2_abs(k))
    enddo
    cs_cipsi_continue = (s_max > cs_abs_sum_max)
  else
    cs_cipsi_continue = (pt2_stop_value(pt2_data) > pt2_max) .and.                     &
        (maxval(abs(pt2_data % variance(1:N_states))) > variance_max)
  endif

end
