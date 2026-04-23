subroutine overlap_cap(psi_cap)

  implicit none
  BEGIN_DOC
  ! Computes and prints the overlap matrix between the real part of the
  ! CAP wave function coefficients and the unperturbed (real) wave function.
  !
  ! For each pair of states (i,j), computes:
  !   overlap(i,j) = |sum_k psi_coef(k,j) * Re(psi_cap(k,i))|
  !
  ! Input:
  !   psi_cap(N_det,N_states) : complex CI coefficients of the CAP wavefunction
  !
  ! Output: prints the overlap matrix to stdout.
  END_DOC

  complex*16, intent(in) :: psi_cap(N_det,N_states)
  double precision, allocatable :: overlap(:,:)
  double precision :: tmp
  integer :: states(N_states)
  complex*16 :: res

  integer :: i,j,k

  allocate(overlap(N_states,N_states))

  overlap = 0d0

  do j = 1, N_states
    do i = 1, N_states
      tmp = 0d0
      do k = 1, N_det
        tmp += psi_coef(k,j) * dble(psi_cap(k,i))
      enddo
      overlap(i,j) = dabs(tmp)
    enddo
  enddo 

  do i = 1, N_states
    states(i) = i
  enddo

  write(*,*) ''
  write(*,*) 'Overlap of the real part with the unperturbed wave function:'
  write(*,*) ''
  write(*,'(6X,100(I5,1X))') states(1:N_states)
  do i = 1, N_states
    write(*,'(I5,1X,100(F5.2,1X))') i, overlap(i,1:N_states)
  enddo
  write(*,*) ''

  deallocate(overlap)

end

subroutine overlap_cap_analysis(psi_cap)

  implicit none
  BEGIN_DOC
  ! Hermitian-normalized weight matrix between CAP eigenstates and unperturbed states.
  !
  ! For each pair (i, j):
  !   w(i,j) = |sum_k psi_coef(k,j) * psi_cap(k,i)|^2 / sum_k |psi_cap(k,i)|^2
  !
  ! w(i,j) lies in [0,1] by Cauchy-Schwarz and measures the Hermitian weight of
  ! unperturbed state j in CAP state i.  sum_j w(i,j) <= 1 (= 1 if the N_states
  ! unperturbed states span the CAP eigenstate).
  !
  ! Also prints, for each CAP state, the dominant unperturbed state and the
  ! total weight captured by the N_states unperturbed states.
  !
  ! Input:
  !   psi_cap(N_det, N_states) : complex CI coefficients of the CAP wavefunction
  END_DOC

  complex*16, intent(in) :: psi_cap(N_det, N_states)
  double precision, allocatable :: w(:,:), hnorm(:)
  complex*16 :: proj
  integer :: i, j, k, jmax
  double precision :: wmax

  allocate(w(N_states, N_states), hnorm(N_states))

  ! Hermitian norm of each CAP eigenvector: sum_k |psi_cap(k,i)|^2
  hnorm = 0d0
  do i = 1, N_states
    do k = 1, N_det
      hnorm(i) = hnorm(i) + cdabs(psi_cap(k,i))**2
    enddo
  enddo

  ! w(i,j) = |<psi_j | psi_cap_i>|^2 / hnorm(i)
  w = 0d0
  do i = 1, N_states
    do j = 1, N_states
      proj = (0d0, 0d0)
      do k = 1, N_det
        proj = proj + psi_coef(k,j) * psi_cap(k,i)
      enddo
      w(i,j) = cdabs(proj)**2 / hnorm(i)
    enddo
  enddo

  write(*,*) ''
  write(*,*) ' Hermitian weight matrix  w(i,j) = |<j|CAP_i>|^2 / <CAP_i|CAP_i>'
  write(*,*) ' Rows = CAP states, Columns = unperturbed states'
  write(*,*) ''
  write(*,'(9X,100(I7,1X))') (j, j=1,N_states)
  do i = 1, N_states
    write(*,'(A,I4,A,100(F7.4,1X))') ' CAP ', i, ' |', w(i,1:N_states)
  enddo
  write(*,*) ''

  write(*,*) ' Dominant unperturbed state for each CAP state:'
  write(*,*) ''
  write(*,'(A)') '  CAP state | dominant | weight | sum_j w(i,j)'
  write(*,'(A)') ' -----------+----------+--------+-------------'
  do i = 1, N_states
    wmax = 0d0
    jmax = 0
    do j = 1, N_states
      if (w(i,j) > wmax) then
        wmax = w(i,j)
        jmax = j
      endif
    enddo
    write(*,'(I6,7X,I6,4X,F7.4,3X,F7.4)') i, jmax, wmax, sum(w(i,:))
  enddo
  write(*,*) ''

  deallocate(w, hnorm)

end
