 BEGIN_PROVIDER [ double precision, ao_wx_cap, (ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_wy_cap, (ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_wz_cap, (ao_num,ao_num)]
BEGIN_DOC
! * array of the integrals of AO_i * W_x AO_j
!
! * array of the integrals of AO_i * W_y AO_j
!
! * array of the integrals of AO_i * W_z AO_j
END_DOC
 implicit none
  integer :: i,j,n,l
  double precision :: f, tmp
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z
  double precision :: alpha, beta
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: lower_exp_val, dx, c,accu_x,accu_y,accu_z
  dim1=500
  lower_exp_val = 40.d0
  ao_spread_x= 0.d0
  ao_spread_y= 0.d0
  ao_spread_z= 0.d0
  !$OMP PARALLEL DO SCHEDULE(GUIDED) &
  !$OMP DEFAULT(NONE) &
  !$OMP PRIVATE(A_center,B_center,power_A,power_B,&
  !$OMP  overlap_x,overlap_y, overlap_z, overlap, &
  !$OMP  alpha, beta,i,j,dx,tmp,c,accu_x,accu_y,accu_z) &
  !$OMP SHARED(nucl_coord,ao_power,ao_prim_num,center_of_mass, &
  !$OMP  cap_onset_x,cap_onset_y,cap_onset_z, &
  !$OMP  ao_wx_cap,ao_wy_cap,ao_wz_cap,ao_num,ao_coef_normalized_ordered_transp, &
  !$OMP  ao_nucl, ao_expo_ordered_transp,dim1,lower_exp_val)
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   do i= 1,ao_num
    B_center(1) = nucl_coord( ao_nucl(i), 1 )
    B_center(2) = nucl_coord( ao_nucl(i), 2 )
    B_center(3) = nucl_coord( ao_nucl(i), 3 )
    power_B(1)  = ao_power( i, 1 )
    power_B(2)  = ao_power( i, 2 )
    power_B(3)  = ao_power( i, 3 )
    accu_x = 0.d0
    accu_y = 0.d0
    accu_z = 0.d0
    do n = 1,ao_prim_num(j)
     alpha = ao_expo_ordered_transp(n,j)
     do l = 1, ao_prim_num(i)
      c = ao_coef_normalized_ordered_transp(n,j)*ao_coef_normalized_ordered_transp(l,i)
      beta = ao_expo_ordered_transp(l,i)

      call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)

      call overlap_bourrin_spread_cap(A_center(1),B_center(1),alpha,beta,power_A(1),power_B(1),tmp,lower_exp_val,dx,dim1,center_of_mass(1),cap_onset_x)
      accu_x +=  c*tmp*overlap_y*overlap_z
      call overlap_bourrin_spread_cap(A_center(2),B_center(2),alpha,beta,power_A(2),power_B(2),tmp,lower_exp_val,dx,dim1,center_of_mass(2),cap_onset_y)
      accu_y +=  c*tmp*overlap_x*overlap_z
      call overlap_bourrin_spread_cap(A_center(3),B_center(3),alpha,beta,power_A(3),power_B(3),tmp,lower_exp_val,dx,dim1,center_of_mass(3),cap_onset_z)
      accu_z +=  c*tmp*overlap_y*overlap_x
     enddo
    enddo
    ao_wx_cap(i,j) = accu_x
    ao_wy_cap(i,j) = accu_y
    ao_wz_cap(i,j) = accu_z
   enddo
  enddo
  !$OMP END PARALLEL DO
END_PROVIDER

subroutine overlap_bourrin_spread_cap(A_center,B_center,alpha,beta,power_A,power_B,overlap_x,lower_exp_val,dx,nx,barycenter_x,cap_x)
  BEGIN_DOC
  ! Computes the following integral :
  !  int [-infty ; +infty] of [(x-A_center)^(power_A) * (x-B_center)^power_B * exp(-alpha(x-A_center)^2) * exp(-beta(x-B_center)^2) * (|x| - x_center)^2 if  abs(x) > x_center ]
  !  needed for the dipole and those things
  END_DOC
  implicit none
  integer :: i,j,k,l
  integer,intent(in) :: power_A,power_B
  double precision, intent(in) :: lower_exp_val,barycenter_x,cap_x
  double precision,intent(in) :: A_center, B_center,alpha,beta
  double precision, intent(out) :: overlap_x,dx
  integer, intent(in) :: nx
  double precision :: x_min,x_max,domain,x,factor,dist,p,p_inv,rho
  double precision :: P_center,pouet_timy
  if(power_A.lt.0.or.power_B.lt.0)then
    overlap_x = 0.d0
    dx = 0.d0
    return
  endif
  p = alpha + beta
  p_inv= 1.d0/p
  rho = alpha * beta * p_inv
  dist = (A_center - B_center)*(A_center - B_center)
  P_center = (alpha * A_center + beta * B_center) * p_inv
  factor = dexp(-rho * dist)
  if(factor.lt.0.000001d0)then
  ! print*,'factor = ',factor
    dx = 0.d0
    overlap_x = 0.d0
    return
  endif
  pouet_timy = dsqrt(lower_exp_val/p)
  x_min = P_center - pouet_timy
  x_max = P_center + pouet_timy
  domain = x_max-x_min
  dx = domain/dble(nx)
  overlap_x = 0.d0
  x = x_min
  do i = 1, nx
    x += dx
    if (dsqrt((x - barycenter_x)**2) >= cap_x) then
      overlap_x += (x-A_center)**(power_A) * (x-B_center)**(power_B) * dexp(-p * (x-P_center)*(x-P_center)) * (dabs(x-barycenter_x) - cap_x)**2
    endif
  enddo
  overlap_x *= factor * dx

end


  BEGIN_PROVIDER [ double precision, ao_spread_x_exact, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_spread_y_exact, (ao_num,ao_num)]
 &BEGIN_PROVIDER [ double precision, ao_spread_z_exact, (ao_num,ao_num)]
 BEGIN_DOC
 ! * array of the integrals of AO_i * x^2 AO_j
 !
 ! * array of the integrals of AO_i * y^2 AO_j
 !
 ! * array of the integrals of AO_i * z^2 AO_j
 END_DOC
 implicit none
  integer :: i,j,n,l
  double precision :: f, tmp
  integer :: dim1
  double precision :: overlap, overlap_x, overlap_y, overlap_z, overlap_gaussian_x
  double precision :: alpha, beta
  integer :: power_Ax,power_Ay,power_Az,power_Bx,power_By,power_Bz
  double precision :: overlap_x0, overlap_x1, overlap_x2
  double precision :: overlap_y0, overlap_y1, overlap_y2
  double precision :: overlap_z0, overlap_z1, overlap_z2
  double precision :: A_center(3), B_center(3)
  integer :: power_A(3), power_B(3)
  double precision :: lower_exp_val, dx, c,accu_x,accu_y,accu_z
  dim1=500
  lower_exp_val = 40.d0
  ao_spread_x_exact= 0.d0
  ao_spread_y_exact= 0.d0
  ao_spread_z_exact= 0.d0
  do j=1,ao_num
   A_center(1) = nucl_coord( ao_nucl(j), 1 )
   A_center(2) = nucl_coord( ao_nucl(j), 2 )
   A_center(3) = nucl_coord( ao_nucl(j), 3 )
   power_A(1)  = ao_power( j, 1 )
   power_A(2)  = ao_power( j, 2 )
   power_A(3)  = ao_power( j, 3 )
   do i= 1,ao_num
    B_center(1) = nucl_coord( ao_nucl(i), 1 )
    B_center(2) = nucl_coord( ao_nucl(i), 2 )
    B_center(3) = nucl_coord( ao_nucl(i), 3 )
    power_B(1)  = ao_power( i, 1 )
    power_B(2)  = ao_power( i, 2 )
    power_B(3)  = ao_power( i, 3 )
    accu_x = 0.d0
    accu_y = 0.d0
    accu_z = 0.d0
    do n = 1,ao_prim_num(j)
     alpha = ao_expo_ordered_transp(n,j)
     do l = 1, ao_prim_num(i)
      c = ao_coef_normalized_ordered_transp(n,j)*ao_coef_normalized_ordered_transp(l,i)
      beta = ao_expo_ordered_transp(l,i)

      call overlap_gaussian_xyz(A_center,B_center,alpha,beta,power_A,power_B,overlap_x0,overlap_y0,overlap_z0,overlap,dim1)

      power_Ax = power_A(1) + 1
      power_Bx = power_B(1)
      overlap_x1 = overlap_gaussian_x(A_center(1), B_center(1), alpha, beta, power_Ax, power_Bx, dim1)

      power_Ax = power_A(1) + 2
      power_Bx = power_B(1)
      overlap_x2 = overlap_gaussian_x(A_center(1), B_center(1), alpha, beta, power_Ax, power_Bx, dim1)
      accu_x +=  c*(overlap_x2 + 2d0 * A_center(1) * overlap_x1 + A_center(1)**2 * overlap_x0 ) * overlap_y0 * overlap_z0

      power_Ay = power_A(2) + 1
      power_By = power_B(2)
      overlap_y1 = overlap_gaussian_x(A_center(2), B_center(2), alpha, beta, power_Ay, power_By, dim1)

      power_Ay = power_A(2) + 2
      power_By = power_B(2)
      overlap_y2 = overlap_gaussian_x(A_center(2), B_center(2), alpha, beta, power_Ay, power_By, dim1)
      accu_y +=  c*(overlap_y2 + 2d0 * A_center(2) * overlap_y1 + A_center(2)**2 * overlap_y0 ) * overlap_x0 * overlap_z0

      power_Az = power_A(3) + 1
      power_Bz = power_B(3)
      overlap_z1 = overlap_gaussian_x(A_center(3), B_center(3), alpha, beta, power_Az, power_Bz, dim1)

      power_Az = power_A(3) + 2
      power_Bz = power_B(3)
      overlap_z2 = overlap_gaussian_x(A_center(3), B_center(3), alpha, beta, power_Az, power_Bz, dim1)
      accu_z +=  c*(overlap_z2 + 2d0 * A_center(3) * overlap_z1 + A_center(3)**2 * overlap_z0 ) * overlap_x0 * overlap_y0

     enddo
    enddo
    ao_spread_x_exact(i,j) = accu_x
    ao_spread_y_exact(i,j) = accu_y
    ao_spread_z_exact(i,j) = accu_z
   enddo
  enddo
 END_PROVIDER

