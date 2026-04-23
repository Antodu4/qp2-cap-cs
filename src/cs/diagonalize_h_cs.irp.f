program diagonalize_h_cs

  implicit none

  BEGIN_DOC
  ! Entry point for the Complex Scaling diagonalization program.
  !
  ! Calls cs(), which diagonalizes the CS Hamiltonian
  ! H_CS = e^{-i*theta} * (H + (e^{-i*theta} - 1)*T) using the Davidson algorithm,
  ! then prints the complex energies for each state and theta value.
  END_DOC

  call cs()

end
