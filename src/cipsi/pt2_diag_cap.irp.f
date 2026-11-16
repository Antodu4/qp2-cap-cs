! Debug diagnostics of the CAP PT2 (keyword pt2_diag). The accumulation
! (pt2_diag_accumulate) and the PT2DIAG lines (print_pt2_diag_cs) are shared
! with the CS; this file only adds the per-eta record of the CAP-CIPSI.


subroutine print_caprec(r)
  implicit none
  BEGIN_DOC
! Prints one record of the per-eta CAP-CIPSI (debug pt2_diag), line CAPREC :
! eta, state, N_det, CIPSI iterations, E (Re, Im), E+PT2 (Re, Im), then for
! the final space variance, sum|Im e|, S = sum|e|, then for the space of the
! stopping decision variance, |PT2|, sum|Im e|, S, and a flag.
! Same layout as the CSREC line of the per-theta CS-CIPSI (print_csrec), with
! eta written in scientific notation.
! The final space has a PT2 only when the loop stopped on N_det_max; otherwise
! E+PT2 = E and its variance, sum|Im e| and S are -1.
! Flag : ok if the variance is below 1e-3 (final space if known, else space of
! the stopping decision, whose variance is larger), suspect otherwise.
  END_DOC
  double precision, intent(in) :: r(15)
  character(len=8)             :: flag
  double precision             :: var

  var = r(9)
  if (var < 0.d0) var = r(12)
  if (var < 1.d-3) then
    flag = 'ok'
  else
    flag = 'suspect'
  endif
  write(*,'(A,ES11.4,1X,I2,1X,I10,1X,I3,4(1X,F18.10),7(1X,ES11.3),1X,A)') ' CAPREC', r(1), nint(r(2)), &
      nint(r(3)), nint(r(4)), r(5:8), r(9:15), trim(flag)

end
