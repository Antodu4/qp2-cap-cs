subroutine print_summary_cap(e_,pt2_data,pt2_data_err,n_det_,n_configuration_,n_st,s2_)
  use selection_types
  implicit none
  BEGIN_DOC
! Print the extrapolated energy in the output
  END_DOC

  integer, intent(in)            :: n_det_, n_configuration_, n_st
  double precision, intent(in)   :: e_(n_st), s2_(n_st)
  type(pt2_type)  , intent(in)   :: pt2_data, pt2_data_err
  integer                        :: i, k
  integer                        :: N_states_p
  character*(9)                  :: pt2_string
  character*(512)                :: fmt

  if (do_pt2) then
    pt2_string = '        '
  else
    pt2_string = '(approx)'
  endif

  N_states_p = min(N_det_,n_st)

  print *, ''
  print '(A,I12)',  'Summary at N_det = ', N_det_
  print '(A)',      '-----------------------------------'
  print *, ''

  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  write(fmt,*) '(13X,', N_states_p, '(6X,A7,1X,I6,10X))'
  write(*,fmt) ('State',k, k=1,N_states_p)
  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  write(fmt,*) '(A13,', 2*N_states_p, '(1X,F14.8))'
  write(*,fmt) '# E          ', psi_energy_cap(1:N_states_p) + dcmplx(nuclear_repulsion,0d0)
  write(*,fmt) '# PT2 '//pt2_string, (pt2_data % pt2(k), pt2_data % pt2_im(k), k=1,N_states_p)
  write(*,'(A)') '#'
  write(*,fmt) '# E+PT2      ', (dble(psi_energy_cap(k))+nuclear_repulsion+pt2_data % pt2(k),dimag(psi_energy_cap(k))+pt2_data % pt2_im(k), k=1,N_states_p)
  write(fmt,*) '(''# ============'',', N_states_p, '(1X,''=============================''))'
  write(*,fmt)
  print *,  ''

  do k=1, N_states_p
    print*,'* State ',k
    print *,  '< S^2 >         = ', s2_(k)
    print *,  'E               = ', dble(psi_energy_cap(k))+nuclear_repulsion, dimag(psi_energy_cap(k))
    print *,  'PT2 Re          = ', pt2_data % pt2(k), ' +/- ', pt2_data_err % pt2(k)
    print *,  'PT2 Im          = ', pt2_data % pt2_im(k), ' +/- ', pt2_data_err % pt2_im(k)
    print *,  'E+PT2 '//pt2_string//' = ', dble(psi_energy_cap(k))+nuclear_repulsion+pt2_data % pt2(k),dimag(psi_energy_cap(k))+pt2_data % pt2_im(k), ' +/- ', dcmplx(pt2_data_err % pt2(k), pt2_data_err % pt2_im(k))
    print *,  ''
  enddo

!  call print_energy_components()

end subroutine

