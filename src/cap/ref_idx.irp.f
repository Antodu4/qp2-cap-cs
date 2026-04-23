subroutine ref_idx_complex(l_det,l_coef,Ndet,N_st,Nint)

  implicit none
  BEGIN_DOC
  ! Prints the dominant determinants (up to 6) of the CIPSI-CAP wave function
  ! for each state, sorted by decreasing absolute value of the CI coefficient.
  !
  ! For each dominant determinant, the excitation degree with respect to the
  ! HF bitmask is computed, and for single and double excitations the orbital
  ! indices (particle, hole, spin) are printed.
  !
  ! Inputs:
  !   l_det(Nint,2,Ndet) : bitmask representation of the Slater determinants
  !   l_coef(Ndet,N_st)  : complex CI coefficients (c-inner product convention)
  !   Ndet               : number of determinants
  !   N_st               : number of states
  !   Nint               : number of integers per spin-orbital bitstring
  END_DOC

  integer, intent(in) :: Ndet, N_st, Nint
  integer(bit_kind), intent(in) :: l_det(Nint,2,Ndet)
  complex*16, intent(in) :: l_coef(Ndet,N_st)

  integer, allocatable :: iorder(:),orb(:,:)
  double precision, allocatable :: tmp(:)
  integer :: i,s,N, degree, exc(0:2,2,2),h1,p1,h2,p2,s1,s2
  double precision :: phase

  N = min(6,N_det)

  allocate(iorder(Ndet),tmp(Ndet),orb(6,N))

  write(*,'(A)') ''
  write(*,'(A)') ' ==========================================================='
  write(*,'(A)') '  List of dominant determinants in CIPSI-CAP wave function:'
  write(*,'(A)') ' ==========================================================='
  write(*,'(A)') ''

  do s = 1, N_st
    do i = 1, Ndet
      iorder(i) = i
      tmp(i) = - cdabs(l_coef(i,s))
    enddo
    call dsort(tmp,iorder,Ndet)
    write(*,'(A7,I4)') ' State:',s
    write(*,'(10(F10.4))') -tmp(1:N)
    write(*,'(10(I10))') iorder(1:N)
    orb = 0
    do i = 1, N
      call get_excitation_degree(l_det(1,1,iorder(i)),hf_bitmask,degree,Nint)
      orb(:,i) = 0
      if (degree == 0 .or. degree > 2) cycle
      call get_excitation(l_det(1,1,iorder(i)),hf_bitmask,exc,degree,phase,Nint)
      call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
      orb(1,i) = p1
      orb(2,i) = h1
      orb(3,i) = s1
      orb(4,i) = p2
      orb(5,i) = h2
      orb(6,i) = s2
    enddo
    do i = 1, N
      write(*,'(10(A1,I3,A1,I3,A1,I1))', advance='no') ' ', orb(1,i), ' ', orb(2,i), ' ', orb(3,i)
    enddo
    write(*,*) ''
    do i = 1, N
      write(*,'(10(A1,I3,A1,I3,A1,I1))', advance='no') ' ', orb(4,i), ' ', orb(5,i), ' ', orb(6,i)
    enddo
    write(*,*) ''

  enddo
  write(*,'(A)') ' ==========================================================='
  
  deallocate(iorder,tmp,orb)

end

