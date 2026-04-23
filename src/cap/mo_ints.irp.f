 BEGIN_PROVIDER [double precision, mo_wx_cap , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_wy_cap , (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_wz_cap , (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * W_x MO_j for CAP
 ! array of the integrals of MO_i * W_y MO_j for CAP
 ! array of the integrals of MO_i * W_z MO_j for CAP
 END_DOC
 implicit none

  call ao_to_mo(                                                     &
      ao_wx_cap,                                                   &
      size(ao_wx_cap,1),                                           &
      mo_wx_cap,                                                   &
      size(mo_wx_cap,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_wy_cap,                                                   &
      size(ao_wy_cap,1),                                           &
      mo_wy_cap,                                                   &
      size(mo_wy_cap,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_wz_cap,                                                   &
      size(ao_wz_cap,1),                                           &
      mo_wz_cap,                                                   &
      size(mo_wz_cap,1)                                            &
      )
END_PROVIDER

 BEGIN_PROVIDER [double precision, mo_spread_x_exact, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_spread_y_exact, (mo_num,mo_num)]
&BEGIN_PROVIDER [double precision, mo_spread_z_exact, (mo_num,mo_num)]
 BEGIN_DOC
 ! array of the integrals of MO_i * x^2 MO_j
 ! array of the integrals of MO_i * y^2 MO_j
 ! array of the integrals of MO_i * z^2 MO_j
 END_DOC
 implicit none
  call ao_to_mo(                                                     &
      ao_spread_x_exact,                                                   &
      size(ao_spread_x_exact,1),                                           &
      mo_spread_x_exact,                                                   &
      size(mo_spread_x_exact,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_spread_y_exact,                                                   &
      size(ao_spread_y_exact,1),                                           &
      mo_spread_y_exact,                                                   &
      size(mo_spread_y_exact,1)                                            &
      )
  call ao_to_mo(                                                     &
      ao_spread_z_exact,                                                   &
      size(ao_spread_z_exact,1),                                           &
      mo_spread_z_exact,                                                   &
      size(mo_spread_z_exact,1)                                            &
      )
END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_ints_cap_alpha,(mo_num,mo_num)]
&BEGIN_PROVIDER [ double precision, mo_ints_cap_beta,(mo_num,mo_num)]
  implicit none
  BEGIN_DOC
  ! Cap integrals in mo basis.
  END_DOC

  integer :: i,j

  mo_ints_cap_alpha = mo_wx_cap + mo_wy_cap + mo_wz_cap
  mo_ints_cap_beta = mo_ints_cap_alpha

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_one_e_integrals_cap,(mo_num,mo_num)]
  implicit none
  BEGIN_DOC
 ! Cap integrals in mo basis.
  END_DOC

  mo_one_e_integrals_cap = 0.5d0 * (mo_ints_cap_alpha + mo_ints_cap_beta)

END_PROVIDER
