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


subroutine energy_Ftot_phi(phi_d,phi_s,do_grad,e_tot,Fd, Fs)
 implicit none
 double precision, intent(in) :: phi_d(ao_num, n_d_occ), phi_s(ao_num, n_s_occ)
 logical, intent(in) :: do_grad
 double precision, intent(out):: e_tot, Fd(ao_num, n_d_occ), Fs(ao_num, n_d_occ)
 double precision :: pd(ao_num, ao_num), ps(ao_num, ao_num)
 double precision :: coul_d(ao_num, ao_num), exch_d(ao_num, ao_num)
 double precision :: coul_s(ao_num, ao_num), exch_s(ao_num, ao_num)
 double precision :: energy_one_e_dm,energy_two_e_dm,e_one_e,e_two_e

 call mo_to_dm(phi_d, n_d_occ, ao_num, pd)
 call mo_to_dm(phi_s, n_s_occ, ao_num, ps)
 call get_j_k_dm(pd,ps,coul_d,exch_d,coul_s,exch_s)
 e_one_e = energy_one_e_dm(pd,ps)
 e_two_e = energy_two_e_dm(pd,ps,coul_d,exch_d,coul_s,exch_s)
 e_tot = e_one_e + e_two_e
 if(do_grad)then
  call get_Fock_ortho(coul_d,exch_d,coul_s,exch_s, Fd,Fs)
 endif
end

 BEGIN_PROVIDER [ double precision, F_d_ortho, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, F_s_ortho, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, norm_F_d_ortho]
&BEGIN_PROVIDER [ double precision, norm_F_s_ortho]
 implicit none
 double precision :: mat_tmp(ao_num,ao_num)
 call get_Fock_ortho(coulomb_d,exchange_d,coulomb_s,exchange_s, F_d_ortho,F_s_ortho)
 integer :: i,j
 double precision :: accu1, accu2
 accu1 = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   norm_F_d_ortho += F_d_ortho(j,i)**2
   norm_F_s_ortho += F_s_ortho(j,i)**2
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, grad_d_ortho, (ao_num, n_d_occ)]
&BEGIN_PROVIDER [ double precision, grad_s_ortho, (ao_num, n_s_occ)]
&BEGIN_PROVIDER [ double precision, norm_grad_d_ortho]
&BEGIN_PROVIDER [ double precision, norm_grad_s_ortho]
 implicit none
 call sym_rect_mat_mul(F_d_ortho,mo_coef_d_ortho,ao_num,ao_num,n_d_occ,grad_d_ortho)
 call sym_rect_mat_mul(F_s_ortho,mo_coef_s_ortho,ao_num,ao_num,n_s_occ,grad_s_ortho)
 integer :: i,l
 norm_grad_d_ortho = 0.d0
 do l = 1, n_d_occ
  do i = 1, ao_num
   norm_grad_d_ortho += grad_d_ortho(i,l)**2
  enddo
 enddo

 norm_grad_s_ortho = 0.d0
 do l = 1, n_s_occ
  do i = 1, ao_num
   norm_grad_s_ortho += grad_s_ortho(i,l)**2
  enddo
 enddo
END_PROVIDER 

subroutine get_Fock_ortho(coul_d,exch_d,coul_s,exch_s, Fd,Fs)
 implicit none
 double precision, intent(in)  :: coul_d(ao_num, ao_num), exch_d(ao_num, ao_num)
 double precision, intent(in)  :: coul_s(ao_num, ao_num), exch_s(ao_num, ao_num)
 double precision, intent(out) :: Fd(ao_num, ao_num), Fs(ao_num, ao_num)
 double precision, allocatable :: mat_tmp(:,:)
 allocate(mat_tmp(ao_num,ao_num))
 mat_tmp = 0.d0
 mat_tmp += ao_one_e_integrals + 2.d0 * coul_d - exch_d + coul_s - 0.5d0 * exch_s
 call ao_ortho_to_ao(mat_tmp, Fd)
 Fd = 4.d0 * Fd
 
 mat_tmp = 0.d0
 mat_tmp = 0.5d0 * (ao_one_e_integrals + 2.d0 * coul_d - exch_d + coul_s - exch_s)
 call ao_ortho_to_ao(mat_tmp, Fs)
 Fs = 4.d0 * Fs

end
