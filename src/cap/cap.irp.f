subroutine cap()

  implicit none
  BEGIN_DOC
  ! Main CAP driver routine.
  !
  ! Loops over n_steps_cap values of eta_cap (starting at eta_cap, incrementing by
  ! eta_step_size), calls diagonalize_ci_cap for each value, and prints a formatted
  ! table of complex CI energies and first-order corrections for all states.
  ! Aborts with an error message if the CAP parameters are not set correctly.
  END_DOC

  complex*16, allocatable :: psi_cap(:,:), energy(:), l_energy(:,:), corr(:), l_corr(:,:)
  double precision :: write_e(2,N_states), eta_cap_save, eta
  integer :: i,j

  if (do_cap) then
    if (cap_onset_x >= 0d0 .and. cap_onset_y >= 0d0 .and. cap_onset_z >= 0d0 &
      .and. eta_cap >= 0d0 .and. n_steps_cap >= 1) then

      if (n_steps_cap > 1 .and. eta_step_size <= 0d0) then
        print*,'Please set eta_step_size > 0d0 if you want screen a range of eta values'
        call abort()
      endif

      allocate(psi_cap(N_det,N_states),energy(N_states), l_energy(N_states,n_steps_cap))
      allocate(corr(N_states), l_corr(N_states,n_steps_cap))

      if (sum(cdabs(psi_cap_coef)) < 1d-6) then
        psi_cap = dcmplx(psi_coef,0d0)
      else
        psi_cap = psi_cap_coef
      endif

      print*,'eta_cap',eta_cap
      eta_cap_save = eta_cap
      do i = 1, n_steps_cap

        call diagonalize_ci_cap(psi_cap,energy,corr)
        l_energy(1:N_states,i) = energy(1:N_states)
        l_corr(1:N_states,i) = corr(1:N_states)

        if (n_steps_cap > 1) then
          eta_cap += eta_step_size
          touch eta_cap
        endif
      enddo
      eta_cap = eta_cap_save
      touch eta_cap

      ! Energy
      write(*,*) ''
      write(*,*) ' ####################'
      write(*,*) ' ### CAP energies ###'
      write(*,*) ' ####################'
      write(*,*) ''
      write(*,'(A,I8)') ' Number of eta values: ', n_steps_cap
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

      write(*,'(A15)',advance='no') '   Eta      '
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

      eta = eta_cap
      do i = 1, n_steps_cap
        do j = 1, N_states
          write_e(1,j) = dble(l_energy(j,i))
          write_e(2,j) = dimag(l_energy(j,i))
        enddo
        write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') eta, write_e(1:2,1:N_states)
        if (n_steps_cap > 1) then
          eta += eta_step_size
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

      write(*,'(A15)',advance='no') '   Eta      '
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

      eta = eta_cap
      do i = 1, n_steps_cap
        do j = 1, N_states
          write_e(1,j) = dble(l_corr(j,i))
          write_e(2,j) = dimag(l_corr(j,i))
        enddo
        write(*,'(ES12.4,3X,100(F16.10,1X,F14.10,2X))') eta, write_e(1:2,1:N_states)
        if (n_steps_cap > 1) then
          eta += eta_step_size
        endif
      enddo

      write(*,'(A15)',advance='no') '============== '
      do i = 1, N_states-1
        write(*,'(A33)',advance='no') ' =============================== '
      enddo
      write(*,'(A33)') ' =============================== '
      write(*,*) ''

      deallocate(psi_cap,energy,l_energy,l_corr)

    else

      print*,''
      print*,'Parameters for cap calculation not set correctly.'
      if (cap_onset_x < 0d0) then
        print*,'cap_onset_x must be >= 0d0'
      endif
      if (cap_onset_y < 0d0) then
        print*,'cap_onset_y must be >= 0d0'
      endif
      if (cap_onset_z < 0d0) then
        print*,'cap_onset_z must be >= 0d0'
      endif
      if (eta_cap < 0d0) then
        print*,'eta_cap must be >= 0d0'
      endif
      if (n_steps_cap < 0d0) then
        print*,'n_steps_cap must be >= 1'
      endif
      print*,''
      call abort()

    endif
  endif
end
