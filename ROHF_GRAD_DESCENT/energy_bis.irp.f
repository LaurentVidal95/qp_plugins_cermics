BEGIN_PROVIDER [ double precision, energy_one_e_pd_ps]
 implicit none
 integer :: i,j
 energy_one_e_pd_ps = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
!  energy_one_e_pd_ps += ao_one_e_integrals_ortho(j,i) * ( 2.d0 * density_mat_d_ortho(j,i) + density_mat_s_ortho(j,i) )
   energy_one_e_pd_ps += ao_one_e_integrals(j,i) * ( 2.d0 * density_mat_d(j,i) + density_mat_s(j,i) )
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, energy_two_e_pd_ps]
 implicit none
 integer :: i,j
 energy_two_e_pd_ps = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   energy_two_e_pd_ps += (2.D0 * coulomb_d(j,i) - exchange_d(j,i)) * (density_mat_d(j,i) + density_mat_s(j,i)) & 
                      +  0.5d0 * (coulomb_s(j,i) - exchange_s(j,i)) * density_mat_s(j,i)
  enddo 
 enddo

END_PROVIDER 

BEGIN_PROVIDER [ double precision, gradient_d_ortho, (ao_num, ao_num)]
 implicit none
 double precision :: mat_tmp(ao_num,ao_num)
!mat_tmp = ao_one_e_integrals + 2.d0 * coulomb_d - exchange_d + coulomb_s - 0.5d0 * exchange_s
 integer :: i,j
 mat_tmp = 0.d0
! do i = 1, ao_num
!  do j = 1, ao_num
!   mat_tmp(j,i) += ao_one_e_integrals(j,i) + 2.d0 * coulomb_d(j,i) - exchange_d(j,i) + coulomb_s(j,i) - 0.5d0 * exchange_s(j,i)
!  enddo
! enddo
 mat_tmp += ao_one_e_integrals + 2.d0 * coulomb_d - exchange_d + coulomb_s - 0.5d0 * exchange_s
 call ao_ortho_to_ao(mat_tmp, gradient_d_ortho)
 gradient_d_ortho = 2.d0 * gradient_d_ortho
END_PROVIDER 

BEGIN_PROVIDER [ double precision, gradient_s_ortho, (ao_num, ao_num)]
 implicit none
 double precision :: mat_tmp(ao_num,ao_num)
 mat_tmp = 0.5d0 * (ao_one_e_integrals + 2.d0 * coulomb_d - exchange_d + coulomb_s - exchange_s)
 call ao_ortho_to_ao(mat_tmp, gradient_s_ortho)
 gradient_s_ortho = 2.d0 * gradient_s_ortho
END_PROVIDER 

