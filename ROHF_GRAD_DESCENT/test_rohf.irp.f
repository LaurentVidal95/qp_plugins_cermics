program test_rohf
! call test_energy_grad
! call test_projection
! call test_retraction
! call test_complex_mat
! call test_exp_mat
! provide proj_virt_prov_generic
! provide U_proj_virt_generic
 provide retract_d_generic
end

subroutine test_retraction
 implicit none
 print*,' norm_M_r_generic = ',norm_M_r_generic
 print*,' norm_X_r_generic = ',norm_x_r_generic
end

subroutine test_projection
 implicit none
 print*,'norm_proj_d_generic = ',norm_proj_d_generic
 print*,'norm_proj_s_generic = ',norm_proj_s_generic
 print*,''
 print*,'norm_proj_d_prov    = ',norm_proj_d_prov
 print*,'norm_proj_s_prov    = ',norm_proj_s_prov
 provide proj_d_test
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
 print*,'norm_X_mat_t_prov =',norm_X_mat_t_prov
 print*,''
 print*,'norm_X_mat_generic   =',norm_X_mat_generic
 print*,'norm_X_mat_t_generic =',norm_X_mat_t_generic
 print*,''
 print*,'trace_dm_d_ortho     =',trace_dm_d_ortho
 print*,'n_d_occ              =',n_d_occ
 print*,'trace_dm_s_ortho     =',trace_dm_s_ortho
 print*,'n_s_occ              =',n_s_occ
 print*,''
 print*,'norm_proj_d_prov     =',norm_proj_d_prov
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
