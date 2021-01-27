program pouet
 implicit none
 integer :: i,j
 double precision :: accu_d, accu_s
 accu_d = 0.d0
 accu_s = 0.d0
 do i = 1, ao_num
  accu_d += density_mat_d_ortho(i,i)
  accu_s += density_mat_s_ortho(i,i)
 enddo
 print*,'accu_d = ',accu_d
 print*,'accu_s = ',accu_s
 double precision, allocatable :: mat_ao(:,:)
 allocate( mat_ao(ao_num, ao_num) )
 call ao_ortho_to_ao(density_mat_d_ortho, mat_ao )
 accu_d = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   accu_d += dabs(mat_ao(j,i) - density_mat_d(j,i))
  enddo
 enddo
 print*,'accu_d = ',accu_d
 print*,''
 print*,''
 print*,''
 print*,'Test on the energy'
 print*,''
 print*,'One e part '
 print*,''
 print*,'One e energy Pd,Ps  = ',energy_one_e_pd_ps
 double precision :: accu_one_e
 accu_one_e = 0.d0
 do i = 1, elec_beta_num
  accu_one_e += 2.d0 * mo_one_e_integrals(i,i) 
 enddo
 do i = elec_beta_num+1, elec_alpha_num
  accu_one_e += 1.d0 * mo_one_e_integrals(i,i) 
 enddo
 print*,'Test usual          = ',accu_one_e
 print*,'psi_energy_h_core   = ',psi_energy_h_core
 print*,''
 print*,'Two e part'
 print*,''
 print*,'energy_two_e_pd_ps  = ',energy_two_e_pd_ps
 print*,'psi_energy_two_e    = ',psi_energy_two_e
 print*,''
 print*,'Test for gradients '
 print*,''
 double precision :: norm_of_mat
 print*,'norm_of_mat(grad_d) = ',norm_of_mat(gradient_d_ortho,ao_num)
 print*,'norm_of_mat(grad_s) = ',norm_of_mat(gradient_s_ortho,ao_num)
 print*,''
 double precision, allocatable :: C_stupid(:,:),C(:,:)
 double precision :: accu_mat_mul, accu_mat_mul_stupid 
 allocate(C_stupid(ao_num,ao_num),C(ao_num,ao_num))
 call sym_square_mat_mul_stupid(sqrt_overlap,inv_sqrt_overlap,ao_num,C_stupid)
 call sym_square_mat_mul(sqrt_overlap,inv_sqrt_overlap,ao_num,C)
 accu_mat_mul_stupid = 0.d0
 accu_mat_mul = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   accu_mat_mul_stupid += dabs(C_stupid(j,i))
   accu_mat_mul        += dabs(C(j,i))
  enddo
 enddo
 print*,'accu_mat_mul_stupid = ',accu_mat_mul_stupid
 print*,'accu_mat_mul        = ',accu_mat_mul
 print*,'ao_num              = ',ao_num
end


subroutine sym_square_mat_mul_stupid(A,B,n,C)
 implicit none
 double precision, intent(in) :: A(n,n),B(n,n)
 double precision, intent(out):: C(n,n)
 integer, intent(in) :: n
 integer :: i,j,k
 C = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   do k = 1, ao_num
    C(j,i) += A(j,k) * B(k,j)
   enddo
  enddo
 enddo
end
