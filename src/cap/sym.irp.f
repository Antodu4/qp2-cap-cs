program sym
  implicit none

  read_wf = .True.
  touch read_wf

  call run()

end

subroutine run()
  implicit none
  BEGIN_DOC
  ! Computes and prints a lower-triangular overlap matrix between CI states
  ! based on the absolute value of pairwise products of CI coefficients:
  !   abs_overlap(i,j) = sum_k |psi_coef(k,i) * psi_coef(k,j)|
  !
  ! Used to identify near-degenerate or symmetry-related states in the
  ! CIPSI wave function.
  END_DOC

  double precision :: tmp
  double precision, allocatable :: abs_overlap(:,:)

  integer :: i,j,k
  logical :: used(N_states)

  allocate(abs_overlap(N_states,N_states))

  abs_overlap = 0d0

  do j = 1, N_states
     do i = j, N_states
        tmp = 0d0
        do k = 1, N_det
          tmp += dabs(psi_coef(k,i) * psi_coef(k,j))
        enddo
        abs_overlap(i,j) = tmp
     enddo
  enddo
  do i = 1, N_states
    write(*,'(I8,100(ES12.2))') i, abs_overlap(i,:)
  enddo

end
