subroutine cs_dav_trajectory(action, e_tot, s2, accepted)
  implicit none
  BEGIN_DOC
  ! Validation of the converged roots of the complex Davidson (cs_dav_validate).
  !
  ! action = 'check'  : accepted is .false. if, for one of the N_states states,
  !                     Im(E) > cs_sf_im_tol, or |S2 - expected_s2| >= cs_sf_s2_tol,
  !                     or |E - E_pred| > cs_dav_extrap_tol * |E(theta-d) - E(theta-2d)|,
  !                     where E_pred extrapolates the energies accepted at the
  !                     previous thetas: linear 2 E(theta-d) - E(theta-2d), or
  !                     quadratic 3 E(theta-d) - 3 E(theta-2d) + E(theta-3d) when
  !                     three points are available. The extrapolation is only used
  !                     with a constant step d and the same N_det at all the points
  !                     (scan of theta in a fixed space); otherwise only Im(E) and
  !                     S2 are checked.
  ! action = 'record' : adds e_tot at theta_cs to the trajectory. A second record
  !                     at the same theta replaces the first one; a theta smaller
  !                     than the last one starts a new trajectory.
  !
  ! e_tot : total energies of the states, E^{CS}_el + e^{-i*theta} * V_nn
  END_DOC
  character*(*), intent(in) :: action
  complex*16, intent(in)    :: e_tot(N_states), s2(N_states)
  logical, intent(out)      :: accepted

  integer, parameter :: n_keep = 4
  integer, save      :: n_traj = 0
  double precision, save :: traj_theta(n_keep)
  integer, save      :: traj_ndet(n_keep)
  complex*16, allocatable, save :: traj_e(:,:)

  integer          :: l, k, n_prev, i_prev(3)
  double precision :: d, tol, de
  complex*16       :: e_pred
  character*(32)   :: why, kind_pred

  if (.not. allocated(traj_e)) then
    allocate(traj_e(N_states,n_keep))
  endif
  if (size(traj_e,1) /= N_states) then
    deallocate(traj_e)
    allocate(traj_e(N_states,n_keep))
    n_traj = 0
  endif

  if (action == 'record') then
    accepted = .True.
    if (n_traj > 0) then
      if (dabs(theta_cs - traj_theta(n_traj)) < 1.d-10) then
        n_traj -= 1
      else if (theta_cs < traj_theta(n_traj)) then
        n_traj = 0
      endif
    endif
    if (n_traj == n_keep) then
      do k = 1, n_keep-1
        traj_theta(k)  = traj_theta(k+1)
        traj_ndet(k)   = traj_ndet(k+1)
        traj_e(:,k)    = traj_e(:,k+1)
      enddo
      n_traj -= 1
    endif
    n_traj += 1
    traj_theta(n_traj) = theta_cs
    traj_ndet(n_traj)  = N_det
    traj_e(1:N_states,n_traj) = e_tot(1:N_states)
    return
  endif

  ! Previous thetas of the trajectory (a point at the current theta is skipped)
  n_prev = 0
  do k = n_traj, 1, -1
    if (theta_cs - traj_theta(k) < 1.d-10) cycle
    if (n_prev == 3) exit
    n_prev += 1
    i_prev(n_prev) = k
  enddo

  ! Number of usable points: constant step d, same N_det
  if (n_prev >= 1) then
    d = theta_cs - traj_theta(i_prev(1))
    if (traj_ndet(i_prev(1)) /= N_det) n_prev = 0
  endif
  do k = 2, n_prev
    if ((dabs(traj_theta(i_prev(k-1)) - traj_theta(i_prev(k)) - d) > 1.d-6*d) .or. &
        (traj_ndet(i_prev(k)) /= N_det)) then
      n_prev = k-1
      exit
    endif
  enddo

  accepted = .True.
  do l = 1, N_states
    why = ''
    if (dimag(e_tot(l)) > cs_sf_im_tol) why = trim(why)//' Im'
    if (only_expected_s2) then
      if (dabs(dble(s2(l)) - expected_s2) >= cs_sf_s2_tol) why = trim(why)//' S2'
    endif

    kind_pred = 'none'
    e_pred = (0.d0,0.d0)
    tol = 0.d0
    de  = 0.d0
    if ((cs_dav_extrap_tol > 0.d0) .and. (n_prev >= 2)) then
      if (n_prev >= 3) then
        e_pred = 3.d0*traj_e(l,i_prev(1)) - 3.d0*traj_e(l,i_prev(2)) + traj_e(l,i_prev(3))
        kind_pred = 'quadratic'
      else
        e_pred = 2.d0*traj_e(l,i_prev(1)) - traj_e(l,i_prev(2))
        kind_pred = 'linear'
      endif
      tol = cs_dav_extrap_tol * cdabs(traj_e(l,i_prev(1)) - traj_e(l,i_prev(2)))
      de  = cdabs(e_tot(l) - e_pred)
      if (de > tol) why = trim(why)//' extrapolation'
    endif

    if (len_trim(why) > 0) accepted = .False.
    if (cs_dav_debug .or. (len_trim(why) > 0)) then
      write(6,'(A,F8.4,A,I2,A,2F16.10,A,F8.4,A,A)') ' DV validate theta', theta_cs, ' st', l, &
          '  E =', e_tot(l), '  S2 =', dble(s2(l)), '  prediction: ', trim(kind_pred)
      if (kind_pred /= 'none') then
        write(6,'(A,2F16.10,A,ES10.2,A,ES10.2)') ' DV validate   E_pred =', e_pred,           &
            '  |E-E_pred| =', de, '  tol =', tol
      endif
      if (len_trim(why) > 0) then
        write(6,'(A,A)') ' DV validate   rejected:', trim(why)
      else
        write(6,'(A)') ' DV validate   accepted'
      endif
    endif
  enddo

end
