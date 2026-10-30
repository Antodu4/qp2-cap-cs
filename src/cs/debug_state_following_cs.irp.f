! Debug of the state following of the complex Davidson (cs_debug_state_following).
! These routines only read the Davidson data and print lines starting with SF:
! they never change the trajectory of the diagonalization.

double precision function dbg_sf_ovh(u,v,sze)
  implicit none
  BEGIN_DOC
  ! Normalized Hermitian overlap |<u|v>| / (||u|| ||v||), between 0 and 1.
  END_DOC
  integer, intent(in)    :: sze
  complex*16, intent(in) :: u(sze), v(sze)

  integer                :: i
  complex*16             :: uv
  double precision       :: uu, vv

  uv = (0d0,0d0)
  uu = 0d0
  vv = 0d0
  do i = 1, sze
    uv = uv + dconjg(u(i)) * v(i)
    uu = uu + dble(u(i))**2 + dimag(u(i))**2
    vv = vv + dble(v(i))**2 + dimag(v(i))**2
  enddo

  if (uu*vv > 0d0) then
    dbg_sf_ovh = cdabs(uv) / dsqrt(uu*vv)
  else
    dbg_sf_ovh = 0d0
  endif

end

integer function dbg_sf_argmax(v,mask,n)
  implicit none
  BEGIN_DOC
  ! Index of the largest v(k) among the k such that mask(k), 0 if there is none.
  END_DOC
  integer, intent(in)          :: n
  double precision, intent(in) :: v(n)
  logical, intent(in)          :: mask(n)

  integer                      :: k
  double precision             :: vmax

  dbg_sf_argmax = 0
  vmax = -1d0
  do k = 1, n
    if (mask(k) .and. (v(k) > vmax)) then
      vmax = v(k)
      dbg_sf_argmax = k
    endif
  enddo

end

subroutine dbg_sf_legend
  implicit none
  BEGIN_DOC
  ! Prints the meaning of the columns of the SF debug lines.
  END_DOC

  write(6,'(A)') ' SFH Debug of the state following (cs_debug_state_following)'
  write(6,'(A)') ' SFH  SF   cyc/it : Davidson cycle (restart) and iteration, as in the Davidson table'
  write(6,'(A)') ' SFH       choice : candidate selected by the current code (max |score|)'
  write(6,'(A)') ' SFH       C1     : max ov_0 (frozen input vector as reference)'
  write(6,'(A)') ' SFH       C2     : max ov_H (Hermitian overlap, no c-normalization)'
  write(6,'(A)') ' SFH       C3     : max |score| with |Re S2 - expected| < 0.05'
  write(6,'(A)') ' SFH       C123   : max ov_0 with |Re S2 - expected| < 0.05'
  write(6,'(A)') ' SFH  rows: rank by |score|, candidate k, E, Re S2, ok (S2 filter), tight (< 0.05),'
  write(6,'(A)') ' SFH        score (= selection criterion), |z| = |psi^T psi| before c-normalization,'
  write(6,'(A)') ' SFH        hn = <psi|psi> (should be 1), ov_H = |<psi|u_in>|, ov_0 = |<psi|u_0>|,'
  write(6,'(A)') ' SFH        res (same definition as the Davidson residual), tags: * choice, 1 2 3 A = C1 C2 C3 C123'
  write(6,'(A)') ' SFH  SFR  at each restart: |<u_new|u_old>| and |<u_new|u_0>| for the reference vector'
  write(6,'(A)') ' SFH  SFF  end of diagonalize_ci_cs: |<psi_out|psi_in>|'

end

subroutine dbg_sf_report(l, itertot, iter, sze, n, U, ldu, W, ldw, y, ldy, lambda, s2, H_jj, &
    score, z, hn, oh, o0, ok, chosen, s2_target, e_shift)
  implicit none
  BEGIN_DOC
  ! Prints, for the state l, the candidates of the state following at one iteration
  ! of the complex Davidson, and the choice that alternative rules would make.
  !
  ! n      : number of Ritz vectors (shift2)
  ! U, W   : basis vectors and H_CS applied to them (first n columns used)
  ! y      : Ritz vectors in the basis U
  ! score  : overlap used by the current selection (overlp(:,l))
  ! z, hn  : |psi^T psi| and <psi|psi> of each Ritz vector before c-normalization
  ! oh, o0 : Hermitian overlaps with u_in(:,l) and with the input vector u_0(:,l)
  ! ok     : state_ok before the selection (S2 filter)
  ! chosen : candidate selected by the current code
  ! e_shift: e^{-i theta} V_nn, added to the energies as in the Davidson table
  END_DOC
  integer, intent(in)          :: l, itertot, iter, sze, n, ldu, ldw, ldy, chosen
  complex*16, intent(in)       :: U(ldu,n), W(ldw,n), y(ldy,n), lambda(n), s2(n), H_jj(sze)
  complex*16, intent(in)       :: score(n), e_shift
  double precision, intent(in) :: z(n), hn(n), oh(n), o0(n), s2_target
  logical, intent(in)          :: ok(n)

  integer, parameter           :: n_top = 5
  double precision, parameter  :: tight_s2 = 0.05d0

  logical, allocatable         :: tight(:), shown(:)
  double precision, allocatable :: sc(:)
  complex*16, allocatable      :: psi(:), wpsi(:)
  integer                      :: list(n_top+5), n_list
  integer                      :: i, k, m, rank, k_c1, k_c2, k_c3, k_all, k_best
  integer, external            :: dbg_sf_argmax
  double precision             :: res, best
  complex*16                   :: delta, e
  character(len=8)             :: tag
  character(len=24)            :: flag

  allocate(tight(n), shown(n), sc(n), psi(sze), wpsi(sze))

  do k = 1, n
    sc(k) = cdabs(score(k))
    tight(k) = ok(k) .and. (dabs(dble(s2(k)) - s2_target) < tight_s2)
  enddo

  k_c1  = dbg_sf_argmax(o0, ok, n)
  k_c2  = dbg_sf_argmax(oh, ok, n)
  k_c3  = dbg_sf_argmax(sc, tight, n)
  k_all = dbg_sf_argmax(o0, tight, n)

  flag = ''
  if (chosen == 0) then
    flag = '  <== NO CANDIDATE'
  else if ((k_c1 /= chosen) .or. (k_c2 /= chosen) .or. (k_c3 /= chosen) .or. (k_all /= chosen)) then
    flag = '  <== DIFF'
  endif

  write(6,'(A,I3,A,I3,A,I2,A,I5,A,I5,4(A,I5),A)') ' SF cyc', itertot, ' it', iter-1, &
      ' st', l, ' nsub', n, '  choice', chosen, '  C1', k_c1, '  C2', k_c2,     &
      '  C3', k_c3, '  C123', k_all, trim(flag)

  ! Rows: the n_top best scores, plus the candidates chosen by each rule
  shown = .False.
  do m = 1, min(n_top, n)
    best = -1d0
    k_best = 0
    do k = 1, n
      if ((.not. shown(k)) .and. (sc(k) > best)) then
        best = sc(k)
        k_best = k
      endif
    enddo
    shown(k_best) = .True.
  enddo
  if (chosen > 0) shown(chosen) = .True.
  if (k_c1   > 0) shown(k_c1)   = .True.
  if (k_c2   > 0) shown(k_c2)   = .True.
  if (k_c3   > 0) shown(k_c3)   = .True.
  if (k_all  > 0) shown(k_all)  = .True.

  ! Sorted by decreasing score
  n_list = 0
  do k = 1, n
    if (shown(k)) then
      n_list = n_list + 1
      list(n_list) = k
    endif
  enddo
  do i = 2, n_list
    m = list(i)
    k = i - 1
    do while (k >= 1)
      if (sc(list(k)) >= sc(m)) exit
      list(k+1) = list(k)
      k = k - 1
    enddo
    list(k+1) = m
  enddo

  write(6,'(A)') ' SF  rank      k          Re(E)            Im(E)   Re(S2) ok tight     score       |z|        hn      ov_H      ov_0        res  tags'

  do i = 1, n_list
    k = list(i)
    rank = 1 + count(sc(1:n) > sc(k))

    ! Residual of the candidate, with the same definition as in the Davidson
    call zgemv('N', sze, n, (1d0,0d0), U, ldu, y(1,k), 1, (0d0,0d0), psi, 1)
    call zgemv('N', sze, n, (1d0,0d0), W, ldw, y(1,k), 1, (0d0,0d0), wpsi, 1)
    res = 0d0
    do m = 1, sze
      if (dabs(dble(H_jj(m) - lambda(k))) >= 1d-2) then
        delta = (lambda(k) * psi(m) - wpsi(m)) / (H_jj(m) - lambda(k))
      else
        delta = (lambda(k) * psi(m) - wpsi(m))                                   &
            / dcmplx(dsign(1d-2, dble(H_jj(m) - lambda(k))), dimag(H_jj(m) - lambda(k)))
      endif
      res = res + dble(delta)**2 + dimag(delta)**2
    enddo

    tag = ''
    if (k == chosen) tag = trim(tag)//'*'
    if (k == k_c1)   tag = trim(tag)//'1'
    if (k == k_c2)   tag = trim(tag)//'2'
    if (k == k_c3)   tag = trim(tag)//'3'
    if (k == k_all)  tag = trim(tag)//'A'

    e = lambda(k) + e_shift
    write(6,'(A,I5,1X,I6,2(1X,F16.10),1X,F7.3,2(2X,L1),1X,ES10.3,4(1X,F9.4),1X,ES10.2,2X,A)') &
        ' SF ', rank, k, dble(e), dimag(e), dble(s2(k)), ok(k), tight(k), sc(k), z(k), hn(k), &
        oh(k), o0(k), res, trim(tag)
  enddo

  deallocate(tight, shown, sc, psi, wpsi)

end
