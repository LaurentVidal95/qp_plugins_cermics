BEGIN_PROVIDER [ double precision, energy_one_e_pd_ps]
 implicit none
 double precision :: energy_one_e_dm
 energy_one_e_pd_ps = energy_one_e_dm(density_mat_d,density_mat_s)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, energy_two_e_pd_ps]
 implicit none
 double precision :: energy_two_e_dm
 energy_two_e_pd_ps = energy_two_e_dm(density_mat_d,density_mat_s,coulomb_d,exchange_d,coulomb_s,exchange_s)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, energy_tot_pd_ps]
 implicit none
 double precision :: e_one_e,e_two_e
 call energy_tot_dm(density_mat_d,density_mat_s,energy_tot_pd_ps,e_one_e,e_two_e) 
END_PROVIDER 

double precision function energy_one_e_dm(pd,ps)
 implicit none
 double precision, intent(in) :: pd(ao_num, ao_num), ps(ao_num, ao_num)
 integer :: i,j
 energy_one_e_dm = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   energy_one_e_dm += (2.d0 * pd(j,i) + ps(j,i)) * ao_one_e_integrals(j,i)
  enddo
 enddo
end

double precision function energy_two_e_dm(pd,ps,coul_d,exch_d,coul_s,exch_s)
 implicit none
 double precision, intent(in) :: pd(ao_num, ao_num), ps(ao_num, ao_num)
 double precision, intent(in) :: coul_d(ao_num, ao_num), exch_d(ao_num, ao_num)
 double precision, intent(in) :: coul_s(ao_num, ao_num), exch_s(ao_num, ao_num)
 integer :: i,j
 energy_two_e_dm = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   energy_two_e_dm+= (2.D0 * coul_d(j,i) - exch_d(j,i)) * (pd(j,i) + ps(j,i)) & 
                      +  0.5d0 * (coul_s(j,i) - exch_s(j,i)) * ps(j,i)
  enddo 
 enddo
end

subroutine energy_tot_dm(pd,ps,e_tot,e_one_e,e_two_e)
 implicit none
 double precision, intent(in) :: pd(ao_num, ao_num), ps(ao_num, ao_num)
 double precision, intent(out):: e_tot,e_one_e,e_two_e
 double precision :: coul_d(ao_num, ao_num), exch_d(ao_num, ao_num)
 double precision :: coul_s(ao_num, ao_num), exch_s(ao_num, ao_num)
 double precision :: energy_one_e_dm,energy_two_e_dm
 call get_j_k_dm(pd,ps,coul_d,exch_d,coul_s,exch_s)
 e_one_e = energy_one_e_dm(pd,ps)
 e_two_e = energy_two_e_dm(pd,ps,coul_d,exch_d,coul_s,exch_s)
 e_tot = e_one_e + e_two_e
end

BEGIN_PROVIDER [ double precision, gradient_d_ortho, (ao_num, ao_num)]
 implicit none
 double precision :: mat_tmp(ao_num,ao_num)
 integer :: i,j
 mat_tmp = 0.d0
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

