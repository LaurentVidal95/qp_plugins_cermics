program test_rohf
 call test_energy_grad
end
subroutine test_energy_grad
 implicit none
 print*,'energy_tot_pd_ps  =',energy_tot_pd_ps
 print*,''
 print*,''
 print*,'norm_F_d_ortho    =' , norm_F_d_ortho
 print*,'norm_F_s_ortho    =' , norm_F_s_ortho
 print*,''
 print*,''
 print*,'norm_grad_d_ortho =' , norm_grad_d_ortho
 print*,'norm_grad_s_ortho =' , norm_grad_s_ortho
 print*,''
 print*,''
 print*,'norm_X_mat_prov   =',norm_X_mat_prov
 print*,''
end

subroutine test_phi_coef
 implicit none
 integer :: i,j,k
 double precision, allocatable :: mat_tmp1(:,:), mat_tmp2(:,:)
 allocate(mat_tmp1(ao_num,ao_num),mat_tmp2(ao_num,ao_num))
 call mo_to_dm(mo_coef_d, n_d_occ, ao_num, mat_tmp1)
 mat_tmp2 = 0.d0
 do k = 1, n_d_occ
  do j = 1, ao_num
   do i = 1, ao_num
    mat_tmp2(i,j) += mo_coef_d(i,k) * mo_coef_d(j,k)
   enddo
  enddo
 enddo
 double precision :: accu
 accu = 0.d0
 do i = 1, ao_num
  do j = 1, ao_num
   accu += dabs(mat_tmp2(j,i) - mat_tmp1(j,i))
  enddo
 enddo
 print*,'accu = ',accu
end
