! Debug diagnostics of the CS PT2 (keyword pt2_diag). They do not change any
! result: they accumulate and print quantities that tell whether the PT2 is a
! sum of aligned contributions (as for a real Hamiltonian) or a sum of
! contributions with scattered phases that cancel each other, which makes the
! stopping test on |PT2| unreliable. See the layout of pt2_type % diag in
! selection_types.f90.


subroutine pt2_diag_accumulate(d, e, re_delta, det, i_generator)
  use bitmasks
  use selection_types
  implicit none
  BEGIN_DOC
! Adds the contribution e of the external determinant det, created from the
! generator i_generator, to the diagnostics d(1:pt2_ndiag) of one state.
! re_delta = Re(E - H_det,det): its sign separates the determinants above E
! (re_delta < 0, as for a real ground state) from those below E (re_delta > 0,
! e.g. a discretized continuum below a resonance).
  END_DOC
  double precision, intent(inout) :: d(pt2_ndiag)
  complex*16, intent(in)          :: e
  double precision, intent(in)    :: re_delta
  integer(bit_kind), intent(in)   :: det(N_int,2)
  integer, intent(in)             :: i_generator

  double precision, parameter     :: pi = 3.14159265358979323846d0
  double precision                :: a
  integer                         :: degree, ic, ib, ifam

  a = cdabs(e)
  d(1) += dabs(dble(e))
  d(2) += a
  d(3) += a*a

  ! Class: excitation degree with respect to the generator, sign of Re(Delta)
  call get_excitation_degree(psi_det_generators(1,1,i_generator), det, degree, N_int)
  if (re_delta > 0.d0) then
    ifam = 2
  else
    ifam = 1
  endif
  if (degree <= 1) then
    ic = ifam
  else
    ic = 2 + ifam
  endif
  ib = pt2_diag_iclass + 4*(ic-1)
  d(ib  ) += dble(e)
  d(ib+1) += dimag(e)
  d(ib+2) += a
  d(ib+3) += 1.d0

  ! Phase histogram, weighted by |e|, in ]-pi,pi]
  ib = int( (datan2(dimag(e), dble(e)) + pi) / (2.d0*pi) * dble(pt2_diag_nbin) ) + 1
  ib = min(max(ib,1), pt2_diag_nbin)
  ib = pt2_diag_ihist + (ifam-1)*pt2_diag_nbin + ib - 1
  d(ib) += a

end


subroutine print_pt2_diag_cs(pt2_data, label)
  use selection_types
  implicit none
  BEGIN_DOC
! Prints the diagnostics of the CS PT2 (lines starting with PT2DIAG), for each
! state. label tells which PT2 it is (e.g. iter or final).
!   sums  : S = sum |e|, sum |Re e|, sum |Im e|, sum |e|**2, n_eff = S**2/sum|e|**2,
!           R = |PT2|/S (1 if no cancellation), |Re PT2|/sum|Re e|, variance
!   class : sum Re e, sum Im e, sum |e|, number of determinants, for
!           single/double excitations above (A) or below (B) E
!   histA, histB : sum of |e| in pt2_diag_nbin sectors of arg(e) in ]-pi,pi],
!           for Re(E - H_det,det) < 0 (A) and > 0 (B)
  END_DOC
  type(pt2_type), intent(in)  :: pt2_data
  character(len=*), intent(in):: label
  integer                     :: k, c, ib
  double precision            :: S, s2, neff, R, rre
  character(len=8), parameter :: cname(4) = (/ 'single_A', 'single_B', 'double_A', 'double_B' /)
  character(len=64)           :: fmt

  if (pt2_ndiag == 0) return

  do k = 1, N_states
    S  = pt2_data % diag(2,k)
    s2 = pt2_data % diag(3,k)
    neff = 0.d0
    if (s2 > 0.d0) neff = S*S/s2
    R = 0.d0
    if (S > 0.d0) R = dsqrt(pt2_data % pt2(k)**2 + pt2_data % pt2_im(k)**2) / S
    rre = 0.d0
    if (pt2_data % diag(1,k) > 0.d0) rre = dabs(pt2_data % pt2(k)) / pt2_data % diag(1,k)
    write(*,'(A,A,A,I2,A,9(1X,ES13.5))') ' PT2DIAG ', trim(label), ' st', k,       &
        ' sums S sumRe sumIm sum2 neff R RRe var :', S, pt2_data % diag(1,k),      &
        pt2_data % pt2_norm(k), s2, neff, R, rre, pt2_data % variance(k)
    do c = 1, 4
      ib = pt2_diag_iclass + 4*(c-1)
      write(*,'(A,A,A,I2,A,A,A,4(1X,ES13.5))') ' PT2DIAG ', trim(label), ' st', k, &
          ' class ', cname(c), ' ReSum ImSum AbsSum count :', pt2_data % diag(ib:ib+3,k)
    enddo
    write(fmt,'(A,I3,A)') '(A,A,A,I2,A,', pt2_diag_nbin, '(1X,ES10.3))'
    ib = pt2_diag_ihist
    write(*,fmt) ' PT2DIAG ', trim(label), ' st', k, ' histA :', &
        pt2_data % diag(ib:ib+pt2_diag_nbin-1,k)
    ib = pt2_diag_ihist + pt2_diag_nbin
    write(*,fmt) ' PT2DIAG ', trim(label), ' st', k, ' histB :', &
        pt2_data % diag(ib:ib+pt2_diag_nbin-1,k)
  enddo

end


subroutine write_wf_cs_theta(theta)
  implicit none
  BEGIN_DOC
! Writes the converged wave function of one theta of the per-theta CS-CIPSI
! (keyword cs_write_wf_theta) in work/wf_cs_theta_<theta>.txt :
!   line 1 : theta, N_det, N_int, N_states
!   then, for each determinant : the 2*N_int integers of its alpha and beta
!   bit strings, then Re and Im of psi_cs_coef for each state.
! The file allows to recompute or analyse the PT2 of this space afterwards.
  END_DOC
  double precision, intent(in) :: theta
  integer                      :: i, k
  character(len=256)           :: fname, fmt

  write(fname,'(A,F6.4,A)') trim(ezfio_filename)//'/work/wf_cs_theta_', theta, '.txt'
  write(fmt,'(A,I3,A,I3,A)') '(', 2*N_int, '(1X,I21),', 2*N_states, '(1X,ES24.16))'
  open(unit=17, file=trim(fname))
  write(17,'(ES24.16,3(1X,I12))') theta, N_det, N_int, N_states
  do i = 1, N_det
    write(17,fmt) psi_det(1:N_int,1,i), psi_det(1:N_int,2,i), &
        (dble(psi_cs_coef(i,k)), dimag(psi_cs_coef(i,k)), k=1,N_states)
  enddo
  close(17)
  write(*,'(A,A)') ' CS wave function written in ', trim(fname)

end


subroutine print_csrec(r)
  implicit none
  BEGIN_DOC
! Prints one record of the per-theta CS-CIPSI (debug pt2_diag), line CSREC :
! theta, state, N_det, CIPSI iterations, E (Re, Im), E+PT2 (Re, Im), then for
! the final space variance, sum|Im e|, S = sum|e|, then for the space of the
! stopping decision variance, |PT2|, sum|Im e|, S, and a flag.
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
  write(*,'(A,F7.4,1X,I2,1X,I10,1X,I3,4(1X,F18.10),7(1X,ES11.3),1X,A)') ' CSREC', r(1), nint(r(2)), &
      nint(r(3)), nint(r(4)), r(5:8), r(9:15), trim(flag)

end
