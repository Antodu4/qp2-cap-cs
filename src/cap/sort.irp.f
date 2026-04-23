recursive subroutine recursive_int_sort(l_det,key,sze,Nint,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Recursively sorts an array of Slater determinants (represented as packed
  ! integer bitstrings) lexicographically by the 2*Nint integers that compose them.
  ! The sort is done on integer index idx first, then recursively on subsequent
  ! indices for groups of equal values.
  !
  ! Inputs/outputs:
  !   l_det(2*Nint,sze) : determinant bitstrings, sorted in-place
  !   key(sze)          : integer key array reordered consistently with l_det
  !   sze               : number of determinants
  !   Nint              : number of integers per spin-orbital bitstring
  !   idx               : current integer index to sort on (1-based)
  END_DOC

  integer, intent(in) :: sze, Nint, idx
  integer(bit_kind), intent(inout) :: l_det(2*Nint, sze)
  integer, intent(inout) :: key(sze)

  integer :: i,j,k,l,nb_u
  integer, allocatable :: iorder(:),nu(:),pu(:)

  if (sze == 0) return

  if (idx < 2*Nint) then

    ! Sort
    call multiple_int_sort(l_det,key,sze,Nint,idx)

    allocate(pu(sze),nu(sze))
    ! Unique, nb and position
    call search_unique_int(l_det,sze,Nint,idx,nb_u,nu,pu)

    do i = 1, nb_u
      call recursive_int_sort(l_det(1,pu(i)),key(pu(i)),nu(i),Nint,idx+1)
    enddo
    deallocate(pu,nu)

  else

    ! Sort
    call multiple_int_sort(l_det,key,sze,Nint,idx)

  endif

end

subroutine multiple_int_sort(l_det,key,sze,Nint,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Sorts an array of determinant bitstrings on integer column idx using i8sort.
  ! Both l_det and the associated key array are reordered consistently.
  !
  ! Inputs/outputs:
  !   l_det(2*Nint,sze) : determinant bitstrings, reordered in-place
  !   key(sze)          : key array reordered with l_det
  !   sze               : number of determinants
  !   Nint              : number of integers per spin-orbital bitstring
  !   idx               : column index (1..2*Nint) used as sort key
  END_DOC

  integer, intent(in) :: sze,Nint,idx
  integer(bit_kind), intent(inout) :: l_det(2*Nint,sze)
  integer, intent(inout) :: key(sze)

  integer :: i,j,k,l,val
  integer(bit_kind), allocatable :: tmp(:),tmp_int(:,:)
  integer, allocatable :: iorder(:), tmp_k(:)

  ! Sort
  allocate(tmp(sze),tmp_int(2*Nint,sze),tmp_k(sze),iorder(sze))

  do i = 1, sze
    tmp(i) = l_det(idx,i)
    tmp_int(:,i) = l_det(:,i)
    tmp_k(i) = key(i)
    iorder(i) = i
  enddo

  call i8sort(tmp,iorder,sze)

  do i = 1, sze
    l_det(:,i) = tmp_int(:,iorder(i))
    key(i) = tmp_k(iorder(i))
  enddo

  deallocate(tmp,tmp_int,tmp_k,iorder)
end

subroutine search_unique_int(l_det,sze,Nint,idx,nb_u,nu,pu)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Finds runs of identical values in column idx of l_det (assumed sorted on that column)
  ! and returns the number of unique values (nb_u), the count of each run (nu), and
  ! the starting position of each run (pu).
  !
  ! Inputs:
  !   l_det(2*Nint,sze) : sorted determinant bitstrings
  !   sze               : number of determinants
  !   Nint              : number of integers per spin-orbital bitstring
  !   idx               : column index used for comparison
  !
  ! Outputs:
  !   nb_u              : number of unique values in column idx
  !   nu(sze)           : number of elements in each unique group
  !   pu(sze)           : starting position (1-based) of each unique group
  END_DOC

  integer, intent(in) :: sze,Nint,idx
  integer(bit_kind), intent(in) :: l_det(2*Nint,sze)
  integer(bit_kind) :: val

  integer, intent(out) :: nb_u, nu(sze), pu(sze)

  integer :: i,j,k,l

  ! Unique, nb and position
  k = 1
  pu = 0 ! starting position
  nu = 0 ! nb
  pu(1) = 1
  nu(1) = 1
  val = l_det(idx,1)
  do i = 2, sze
    if (val /= l_det(idx,i)) then
      k = k + 1
      pu(k) = i
      nu(k) = nu(k) + 1
      val = l_det(idx,i)
    else
      nu(k) = nu(k) + 1
    endif
  enddo

  nb_u = k

end

