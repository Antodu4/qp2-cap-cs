subroutine cap_print_traj(algo, n_eta, eta0, deta, nstates, ndet_max, pt2max, nrec, rec, e_eta)
  implicit none
  BEGIN_DOC
  ! Shared CAP trajectory output, used by BOTH CAP-CIPSI drivers
  ! (run_stochastic_cap_cipsi and run_stochastic_cap_cipsi_per_eta) so that the
  ! two algorithms produce an identical layout read by a single extractor.
  ! Mirror of cs_print_traj with eta in place of theta; machine-readable lines
  ! are tagged TRAJ_CAP (instead of TRAJ) to tell them apart from CS output.
  !
  ! Prints a banner followed by 3 blocks:
  !   A: convergence vs N_det  (one row per state,eta,iter)
  !   B: final CAP energies per eta  (format identical to cap())
  !   C: long-format machine-readable TRAJ_CAP records
  !
  ! rec columns: 1=iter 2=n_det 3=eta 4=state 5=ReE 6=ImE 7=ReE+PT2 8=ImE+PT2
  !              9=RePT2 10=ImPT2 11=variance
  !
  ! All inputs are passed explicitly (no providers, no module state) so the
  ! routine needs no explicit interface at the call sites.
  END_DOC

  character(len=*), intent(in) :: algo
  integer, intent(in)          :: n_eta, nstates, ndet_max, nrec
  double precision, intent(in) :: eta0, deta, pt2max
  double precision, intent(in) :: rec(11, nrec)
  complex*16, intent(in)       :: e_eta(nstates, n_eta)

  integer                      :: i, k, kt, jrec, keyidx
  integer, allocatable         :: order(:)
  double precision             :: eta_val

  ! --- Order records by (state, eta, iter) for blocks A and C ---
  if (nrec > 0) then
    allocate(order(nrec))
    do i = 1, nrec
      order(i) = i
    enddo
    do i = 1, nrec-1
      jrec = i
      do k = i+1, nrec
        if ( rec(4,order(k)) <  rec(4,order(jrec)) .or.                                     &
            (rec(4,order(k)) == rec(4,order(jrec)) .and. rec(3,order(k)) <  rec(3,order(jrec))) .or. &
            (rec(4,order(k)) == rec(4,order(jrec)) .and. rec(3,order(k)) == rec(3,order(jrec)) .and. &
             rec(1,order(k)) <  rec(1,order(jrec))) ) then
          jrec = k
        endif
      enddo
      keyidx = order(i) ; order(i) = order(jrec) ; order(jrec) = keyidx
    enddo
  endif

  ! --- Banner ---
  write(*,'(A)') ''
  write(*,'(A)') '================================================================================'
  write(*,'(A,A)') '  CAP SCAN SUMMARY     algorithm = ', trim(algo)
  write(*,'(A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') '  N_states = ', nstates, &
    '     N_eta = ', n_eta, '     eta in [', eta0, ', ', &
    eta0 + dble(n_eta-1)*deta, '] step ', deta
  write(*,'(A,I0,A,ES9.1)') '  stopping: N_det_max = ', ndet_max, '   pt2_max = ', pt2max
  write(*,'(A)') '================================================================================'
  write(*,'(A)') ''

  ! --- Block A: convergence vs N_det ---
  write(*,'(A)') ' ### CAP CONVERGENCE vs N_det ###'
  write(*,'(A)') ' # Determinant-plane trajectory. One row per (state,eta,iter). E, E+PT2, PT2 complex (Ha).'
  write(*,'(A)') '  state     eta       iter   N_det       Re(E)          Im(E)        Re(E+PT2)       Im(E+PT2)       Re(PT2)        Im(PT2)       variance'
  do i = 1, nrec
    jrec = order(i)
    write(*,'(I7,2X,ES12.4,I6,I9,3X,F15.10,1X,F14.10,3X,F15.10,1X,F14.10,3X,ES15.7,1X,ES15.7,2X,ES10.3)') &
      nint(rec(4,jrec)), rec(3,jrec), nint(rec(1,jrec)), nint(rec(2,jrec)), &
      rec(5,jrec), rec(6,jrec), rec(7,jrec), rec(8,jrec), rec(9,jrec), rec(10,jrec), rec(11,jrec)
  enddo
  write(*,'(A)') ''

  ! --- Block B: final CAP energies per eta (format identical to cap()) ---
  write(*,'(A)') ' ### CAP energies (final, converged per eta) ###'
  write(*,'(A)') ' # One row per eta; columns = (Re,Im) of the converged E for each state.'
  do kt = 1, n_eta
    eta_val = eta0 + dble(kt-1)*deta
    write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') eta_val, &
      (dble(e_eta(k,kt)), dimag(e_eta(k,kt)), k=1,nstates)
  enddo
  write(*,'(A)') ''

  ! --- Block C: machine-readable TRAJ_CAP records ---
  write(*,'(A)') ' ### CAP trajectories (long format, machine-readable) ###'
  write(*,'(A)') ' # columns: iter n_det eta state ReE ImE ReEPT2 ImEPT2 RePT2 ImPT2 variance'
  do i = 1, nrec
    jrec = order(i)
    write(*,'(" TRAJ_CAP",1X,I5,1X,I9,1X,ES15.7,1X,I5,1X,F15.10,1X,F14.10,1X,F15.10,1X,F14.10,1X,ES15.7,1X,ES15.7,1X,ES12.4)') &
      nint(rec(1,jrec)), nint(rec(2,jrec)), rec(3,jrec), nint(rec(4,jrec)), &
      rec(5,jrec), rec(6,jrec), rec(7,jrec), rec(8,jrec), rec(9,jrec), rec(10,jrec), rec(11,jrec)
  enddo
  write(*,'(A)') ''

  if (allocated(order)) deallocate(order)

end
