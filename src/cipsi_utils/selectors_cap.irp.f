BEGIN_PROVIDER [ complex*16, psi_coef_cap_sorted, (psi_det_size,N_states) ]
   implicit none
   BEGIN_DOC
   ! Wave function sorted by determinants contribution to the norm (state-averaged)
   END_DOC

   psi_coef_cap_sorted = (0d0,0d0)

END_PROVIDER

BEGIN_PROVIDER [ double precision, psi_coef_tmpsave, (psi_det_size,N_states) ]
   implicit none
   BEGIN_DOC
   ! tmp storage of psi_coef that will be overwritten
   END_DOC

   psi_coef_tmpsave = 0d0

END_PROVIDER

 BEGIN_PROVIDER [ complex*16, psi_selectors_coef_cap, (psi_selectors_size,N_states) ]
  implicit none
  BEGIN_DOC
  ! Determinants on which we apply <i|H|psi> for perturbation.
  END_DOC
  integer                        :: i,k

  do k=1,N_states
    do i=1,N_det_selectors
      psi_selectors_coef_cap(i,k) = psi_coef_cap_sorted(i,k)
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ complex*16, psi_selectors_coef_cap_transp, (N_states,psi_selectors_size) ]
  implicit none
  BEGIN_DOC
  ! Transposed psi_selectors
  END_DOC
  integer                        :: i,k

  do i=1,N_det_selectors
    do k=1,N_states
      psi_selectors_coef_cap_transp(k,i) = psi_selectors_coef_cap(i,k)
    enddo
  enddo
END_PROVIDER

subroutine cap_build_psi_det_sorted(psi_cap)
    implicit none

    complex*16, intent(in) :: psi_cap(N_det,N_states)

    integer :: i,k
    double precision:: n_re, n_im

    n_re = 0d0
    n_im = 0d0
    do k = 1, N_states
        do i = 1, N_det
            n_re += dble(psi_cap(i,k))**2
            n_im += dimag(psi_cap(i,k))**2
        enddo
    enddo

    psi_coef_tmpsave(:,:) = psi_coef(:,:)
    touch psi_coef_tmpsave

    !print*,'n_re',n_re
    n_re = 1d0/dsqrt(n_re)
    if (n_im > 1d-4) then
      n_im = 1d0/dsqrt(n_im)
    else
      n_im = 1d0
    endif

    do k = 1, N_states
        do i = 1, N_det
            psi_coef(i,k) = dabs(n_re * dble(psi_cap(i,k))) + dabs(n_im * dimag(psi_cap(i,k)))
        enddo
    enddo

    touch psi_coef

end

subroutine build_psi_coef_cap_sorted(psi_cap)
    implicit none

    complex*16, intent(in) :: psi_cap(N_det,N_states)
    integer :: i,k

    PROVIDE psi_det_sorted_ordering

    !if (psi_det_size /= N_det) then
    !    print*, psi_det_size, N_det
    !    print*,'Big problem here...'
    !    call abort
    !endif

    call cap_build_psi_det_sorted(psi_cap)

    do k = 1, N_states
        do i = 1, N_det
            psi_coef_cap_sorted(i,k) = psi_cap(psi_det_sorted_ordering(i),k)
        enddo
    enddo

    touch psi_coef_cap_sorted

end

subroutine save_wavefunction_cap(coefs)
  implicit none

  complex*16, intent(out) :: coefs(N_det,N_states)
  integer :: i,k

  do k = 1, N_states
    do i = 1, N_det
      coefs(i,k) = psi_cap_coef(i,k)
    enddo
  enddo

end

subroutine read_wavefunction_cap(coefs,n)
  implicit none

  integer, intent(in) :: n
  complex*16, intent(in) :: coefs(n,N_states)
  integer :: i,k

  do k = 1, N_states
    do i = 1, n
      psi_cap_coef(i,k) = coefs(i,k)
    enddo
  enddo
  do k = 1, N_states
    do i = n+1, N_det
      psi_cap_coef(i,k) = (0d0,0d0)
    enddo
  enddo

  touch psi_cap_coef

end
