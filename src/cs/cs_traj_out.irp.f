subroutine cs_print_traj(algo, n_theta, theta0, dtheta, nstates, ndet_max, pt2max, nrec, rec, e_theta)
  implicit none
  BEGIN_DOC
  ! Shared CS trajectory output, used by BOTH CS-CIPSI drivers
  ! (run_stochastic_cs_cipsi and run_stochastic_cs_cipsi_per_theta) so that the
  ! two algorithms produce a byte-identical layout read by a single extractor.
  !
  ! Prints a banner followed by 3 blocks:
  !   A: convergence vs N_det  (one row per state,theta,iter)
  !   B: final CS energies per theta  (format identical to cs(), parsed by cs_analysis.py)
  !   C: long-format machine-readable TRAJ records
  !
  ! rec columns: 1=iter 2=n_det 3=theta 4=state 5=ReE 6=ImE 7=ReE+PT2 8=ImE+PT2
  !              9=RePT2 10=ImPT2 11=variance
  !
  ! All inputs are passed explicitly (no providers, no module state) so the
  ! routine needs no explicit interface at the call sites.
  END_DOC

  character(len=*), intent(in) :: algo
  integer, intent(in)          :: n_theta, nstates, ndet_max, nrec
  double precision, intent(in) :: theta0, dtheta, pt2max
  double precision, intent(in) :: rec(11, nrec)
  complex*16, intent(in)       :: e_theta(nstates, n_theta)

  integer                      :: i, k, kt, jrec, keyidx
  integer, allocatable         :: order(:)
  double precision             :: theta_val

  ! --- Order records by (state, theta, iter) for blocks A and C ---
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
  write(*,'(A,A)') '  CS SCAN SUMMARY      algorithm = ', trim(algo)
  write(*,'(A,I0,A,I0,A,F12.6,A,F12.6,A,F12.6)') '  N_states = ', nstates, &
    '     N_theta = ', n_theta, '     theta in [', theta0, ', ', &
    theta0 + dble(n_theta-1)*dtheta, '] step ', dtheta
  write(*,'(A,I0,A,ES9.1)') '  stopping: N_det_max = ', ndet_max, '   pt2_max = ', pt2max
  write(*,'(A)') '================================================================================'
  write(*,'(A)') ''

  ! --- Block A: convergence vs N_det ---
  write(*,'(A)') ' ### CS CONVERGENCE vs N_det ###'
  write(*,'(A)') ' # Determinant-plane trajectory. One row per (state,theta,iter). E, E+PT2, PT2 complex (Ha).'
  write(*,'(A)') '  state    theta      iter   N_det       Re(E)          Im(E)        Re(E+PT2)       Im(E+PT2)       Re(PT2)        Im(PT2)       variance'
  do i = 1, nrec
    jrec = order(i)
    write(*,'(I7,2X,ES12.4,I6,I9,3X,F15.10,1X,F14.10,3X,F15.10,1X,F14.10,3X,ES15.7,1X,ES15.7,2X,ES10.3)') &
      nint(rec(4,jrec)), rec(3,jrec), nint(rec(1,jrec)), nint(rec(2,jrec)), &
      rec(5,jrec), rec(6,jrec), rec(7,jrec), rec(8,jrec), rec(9,jrec), rec(10,jrec), rec(11,jrec)
  enddo
  write(*,'(A)') ''

  ! --- Block B: final CS energies per theta (format identical to cs()) ---
  write(*,'(A)') ' ### CS energies (final, converged per theta) ###'
  write(*,'(A)') ' # One row per theta; columns = (Re,Im) of the converged E for each state.'
  do kt = 1, n_theta
    theta_val = theta0 + dble(kt-1)*dtheta
    write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') theta_val, &
      (dble(e_theta(k,kt)), dimag(e_theta(k,kt)), k=1,nstates)
  enddo
  write(*,'(A)') ''

  ! --- Block C: machine-readable TRAJ records ---
  write(*,'(A)') ' ### CS trajectories (long format, machine-readable) ###'
  write(*,'(A)') ' # columns: iter n_det theta state ReE ImE ReEPT2 ImEPT2 RePT2 ImPT2 variance'
  do i = 1, nrec
    jrec = order(i)
    write(*,'(" TRAJ",1X,I5,1X,I9,1X,ES15.7,1X,I5,1X,F15.10,1X,F14.10,1X,F15.10,1X,F14.10,1X,ES15.7,1X,ES15.7,1X,ES12.4)') &
      nint(rec(1,jrec)), nint(rec(2,jrec)), rec(3,jrec), nint(rec(4,jrec)), &
      rec(5,jrec), rec(6,jrec), rec(7,jrec), rec(8,jrec), rec(9,jrec), rec(10,jrec), rec(11,jrec)
  enddo
  write(*,'(A)') ''

  if (allocated(order)) deallocate(order)

end
