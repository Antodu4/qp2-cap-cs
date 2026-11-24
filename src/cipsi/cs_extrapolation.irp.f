! Extrapolation of the complex energy of the CS-CIPSI to the FCI limit.
!
! At each CIPSI iteration, the stationary energy E of the current space and the
! sums of the contributions e_alpha of the external determinants (selection_cs)
! are stored (cs_extrap_add). The real and imaginary parts are then
! extrapolated independently by linear fits on the last k points (k = 2, 3, ...),
! as in Damour et al., J. Phys. Chem. Lett. 15, 8296 (2024) :
!   Re E  against  Re PT2  = sum Re e_alpha          -> 0
!   Im E  against  sum |Im e_alpha|                  -> 0   (Eq. 17)
! with, in addition, the absolute-value version of Eq. 16 for the real part and
! the variance, which do not assume a sign for the contributions :
!   Re E  against  sum |Re e_alpha|                  -> 0   (Eq. 16)
!   Re E, Im E  against  variance sum |<a|H|Psi>|**2 -> 0
! The fits use extrapolate_data (src/utils/extrapolation.irp.f). The energies
! of the calculation are not changed : the extrapolated values are only printed.
!
! Output lines (one per state) :
!   CSEXT    theta it N_det ReE ImE RePT2 ImPT2 sum|Re e| sum|Im e| sum|e| variance Re(sum v^2) Im(sum v^2)
!            followed by the stochastic errors of the 8 sums (lines CSEXTERR), v = <alpha|H|Psi>
!   CSEXTRAP theta st k ReE(RePT2) ImE(sum|Im e|) ReE(sum|Re e|) ReE(var) ImE(var)
!            extrapolated energies with the k last points

subroutine cs_extrap_reset
  use selection_types
  implicit none
  BEGIN_DOC
! Forgets the points of the extrapolation (start of a theta)
  END_DOC
  cs_ext_n = 0
  if (allocated(cs_ext)) cs_ext = 0.d0
end

subroutine cs_extrap_add(pt2_data, pt2_data_err)
  use selection_types
  implicit none
  BEGIN_DOC
! Stores the point of the current CIPSI iteration : energy of the current space
! (psi_energy_cs) and the sums of the PT2 contributions, and prints it (CSEXT).
  END_DOC
  type(pt2_type), intent(in) :: pt2_data, pt2_data_err
  complex*16                 :: eith, e
  integer                    :: s, k

  if (.not. allocated(cs_ext)) then
    allocate(cs_ext(cs_ext_nq, cs_ext_nmax, N_states))
    cs_ext = 0.d0
    cs_ext_n = 0
  endif
  if (size(cs_ext,3) /= N_states) then
    deallocate(cs_ext)
    allocate(cs_ext(cs_ext_nq, cs_ext_nmax, N_states))
    cs_ext = 0.d0
    cs_ext_n = 0
  endif

  ! Keep the cs_ext_nmax last points
  if (cs_ext_n == cs_ext_nmax) then
    do k = 1, cs_ext_nmax-1
      cs_ext(:,k,:) = cs_ext(:,k+1,:)
    enddo
  else
    cs_ext_n = cs_ext_n + 1
  endif

  eith = dcmplx(dcos(theta_cs), -dsin(theta_cs))
  do s = 1, N_states
    e = psi_energy_cs(s) + eith * dcmplx(nuclear_repulsion, 0.d0)
    cs_ext( 1,cs_ext_n,s) = dble(N_det)
    cs_ext( 2,cs_ext_n,s) = dble(e)
    cs_ext( 3,cs_ext_n,s) = dimag(e)
    cs_ext( 4,cs_ext_n,s) = pt2_data % pt2(s)
    cs_ext( 5,cs_ext_n,s) = pt2_data % pt2_im(s)
    cs_ext( 6,cs_ext_n,s) = pt2_data % pt2_re_abs(s)
    cs_ext( 7,cs_ext_n,s) = pt2_data % pt2_norm(s)
    cs_ext( 8,cs_ext_n,s) = pt2_data % pt2_abs(s)
    cs_ext( 9,cs_ext_n,s) = pt2_data % variance(s)
    cs_ext(10,cs_ext_n,s) = pt2_data % variance_c_re(s)
    cs_ext(11,cs_ext_n,s) = pt2_data % variance_c_im(s)
    write(6,'(A,F8.4,A,I3,A,I3,I12,2(1X,F18.10),8(1X,ES16.8))') ' CSEXT    theta', theta_cs, &
        ' st', s, ' it', cs_ext_n, N_det, (cs_ext(k,cs_ext_n,s), k=2,11)
    write(6,'(A,F8.4,A,I3,A,I3,12X,38X,8(1X,ES16.8))') ' CSEXTERR theta', theta_cs, &
        ' st', s, ' it', cs_ext_n, pt2_data_err % pt2(s), pt2_data_err % pt2_im(s),        &
        pt2_data_err % pt2_re_abs(s), pt2_data_err % pt2_norm(s), pt2_data_err % pt2_abs(s), &
        pt2_data_err % variance(s), pt2_data_err % variance_c_re(s), pt2_data_err % variance_c_im(s)
  enddo

end

subroutine cs_extrap_print
  use selection_types
  implicit none
  BEGIN_DOC
! Prints the extrapolated complex energy of each state, with the k last points
! of the current theta (k = 2 .. min(n,8)) : table and CSEXTRAP lines.
  END_DOC
  integer                       :: s, k, n, n0
  double precision, allocatable :: o_re_pt2(:), o_im_aim(:), o_re_are(:), o_re_var(:), o_im_var(:)

  n = min(cs_ext_n, 8)
  if (n < 2) return
  n0 = cs_ext_n - n
  allocate(o_re_pt2(n), o_im_aim(n), o_re_are(n), o_re_var(n), o_im_var(n))

  write(6,'(A)') ''
  write(6,'(A,F8.4)') 'CS extrapolated energies, theta = ', theta_cs
  write(6,'(A)') '----------------------------------------'
  write(6,'(A)') ' Linear fits on the k last CIPSI iterations of this theta :'
  write(6,'(A)') '   (1) Re E vs Re PT2   (2) Im E vs sum|Im e| (Eq. 17)   (3) Re E vs sum|Re e| (Eq. 16)'
  write(6,'(A)') '   (4) Re E vs variance (5) Im E vs variance'
  do s = 1, N_states
    call extrapolate_data(n, cs_ext(2,n0+1:cs_ext_n,s), cs_ext( 4,n0+1:cs_ext_n,s), o_re_pt2)
    call extrapolate_data(n, cs_ext(3,n0+1:cs_ext_n,s), cs_ext( 7,n0+1:cs_ext_n,s), o_im_aim)
    call extrapolate_data(n, cs_ext(2,n0+1:cs_ext_n,s), cs_ext( 6,n0+1:cs_ext_n,s), o_re_are)
    call extrapolate_data(n, cs_ext(2,n0+1:cs_ext_n,s), cs_ext( 9,n0+1:cs_ext_n,s), o_re_var)
    call extrapolate_data(n, cs_ext(3,n0+1:cs_ext_n,s), cs_ext( 9,n0+1:cs_ext_n,s), o_im_var)
    write(6,'(A)') ''
    write(6,'(A,I3)') ' State ', s
    write(6,'(A)') ' ===  ================  ================  ================  ================  ================'
    write(6,'(A)') '  k   (1) Re E          (2) Im E          (3) Re E          (4) Re E          (5) Im E'
    write(6,'(A)') ' ===  ================  ================  ================  ================  ================'
    do k = 2, n
      write(6,'(I4,5(2X,F16.10))') k, o_re_pt2(k), o_im_aim(k), o_re_are(k), o_re_var(k), o_im_var(k)
    enddo
    write(6,'(A)') ' ===  ================  ================  ================  ================  ================'
    do k = 2, n
      write(6,'(A,F8.4,A,I3,A,I3,5(1X,F18.10))') ' CSEXTRAP theta', theta_cs, ' st', s, ' k', k, &
          o_re_pt2(k), o_im_aim(k), o_re_are(k), o_re_var(k), o_im_var(k)
    enddo
  enddo
  write(6,'(A)') ''
  deallocate(o_re_pt2, o_im_aim, o_re_are, o_re_var, o_im_var)

end
