subroutine cs()

  implicit none

  BEGIN_DOC
  ! Main driver for Complex Scaling (CS) calculations.
  !
  ! Diagonalizes the CS Hamiltonian H_CS = e^{-i*theta} * (H + (e^{-i*theta} - 1)*T)
  ! for one or more values of the complex scaling angle theta_cs, collecting
  ! complex energies and energy corrections for each state.
  !
  ! If n_steps_cs > 1, theta is scanned from theta_cs upward in steps of
  ! theta_step_size, and a table of (Re(E), Im(E)) values is printed for
  ! each angle and each electronic state.
  !
  ! The initial guess for the CS coefficients is taken from psi_cs_coef if
  ! it is non-zero, otherwise from the real CI coefficients psi_coef.
  END_DOC

  complex*16, allocatable :: psi_cs(:,:), energy(:), l_energy(:,:), corr(:), l_corr(:,:)
  double precision :: write_e(2,N_states), theta_cs_save, theta
  integer :: i,j

  if (do_cs) then
    if (theta_cs >= 0d0 .and. n_steps_cs >= 1) then

      if (n_steps_cs > 1 .and. theta_step_size <= 0d0) then
        print*,'Please set theta_step_size > 0d0 if you want to screen a range of theta values'
        call abort()
      endif

      allocate(psi_cs(N_det,N_states),energy(N_states), l_energy(N_states,n_steps_cs))
      allocate(corr(N_states), l_corr(N_states,n_steps_cs))

      if (sum(cdabs(psi_cs_coef)) < 1d-6) then
        psi_cs = dcmplx(psi_coef,0d0)
      else
        psi_cs = psi_cs_coef
      endif

      theta_cs_save = theta_cs
      do i = 1, n_steps_cs

        call diagonalize_ci_cs(psi_cs,energy,corr)
        l_energy(1:N_states,i) = energy(1:N_states)
        l_corr(1:N_states,i) = corr(1:N_states)

        if (n_steps_cs > 1) then
          theta_cs += theta_step_size
          touch theta_cs
        endif
      enddo
      theta_cs = theta_cs_save
      touch theta_cs

      ! Energy
      write(*,*) ''
      write(*,*) ' #####################'
      write(*,*) ' ###  CS energies  ###'
      write(*,*) ' #####################'
      write(*,*) ''
      write(*,'(A,I8)') ' Number of theta values: ', n_steps_cs
      write(*,'(A,I8)') ' Number of states: ', N_states
      write(*,*) ''
      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      write(*,'(A15)',advance='no') ' '
      do i = 1, N_states-1
        write(*,'(A17,I4,A12)',advance='no') 'State ',i,' '
      enddo
      write(*,'(A17,I4,A12)') '          State ',N_states,'     '

      write(*,'(A15)',advance='no') '   Theta    '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      write(*,'(A15)',advance='no') ' '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') '   Energy (Re)     Energy (Im)  '
      enddo
      write(*,'(A33)') '   Energy (Re)     Energy (Im)  '

      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      theta = theta_cs
      do i = 1, n_steps_cs
        do j = 1, N_states
          write_e(1,j) = dble(l_energy(j,i))
          write_e(2,j) = dimag(l_energy(j,i))
        enddo
        write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') theta, write_e(1:2,1:N_states)
        if (n_steps_cs > 1) then
          theta += theta_step_size
        endif
      enddo

      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '
      write(*,*) ''

      ! Energy correction
      write(*,*) ''
      write(*,*) ' ###########################'
      write(*,*) ' #### Energy corrections ###'
      write(*,*) ' ###########################'
      write(*,*) ''
      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      write(*,'(A15)',advance='no') ' '
      do i = 1, N_states-1
        write(*,'(A17,I4,A12)',advance='no') 'State ',i,' '
      enddo
      write(*,'(A17,I4,A12)') '          State ',N_states,'     '

      write(*,'(A15)',advance='no') '   Theta    '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      write(*,'(A15)',advance='no') ' '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') '   Energy (Re)     Energy (Im)  '
      enddo
      write(*,'(A33)') '   Energy (Re)     Energy (Im)  '

      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '

      theta = theta_cs
      do i = 1, n_steps_cs
        do j = 1, N_states
          write_e(1,j) = dble(l_corr(j,i))
          write_e(2,j) = dimag(l_corr(j,i))
        enddo
        write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') theta, write_e(1:2,1:N_states)
        if (n_steps_cs > 1) then
          theta += theta_step_size
        endif
      enddo

      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '
      write(*,*) ''

      deallocate(psi_cs,energy,l_energy,l_corr)

    else

      print*,''
      print*,'Parameters for CS calculation not set correctly.'
      if (theta_cs < 0d0) then
        print*,'theta_cs must be >= 0d0'
      endif
      if (n_steps_cs < 1) then
        print*,'n_steps_cs must be >= 1'
      endif
      print*,''
      call abort()

    endif
  endif
end
