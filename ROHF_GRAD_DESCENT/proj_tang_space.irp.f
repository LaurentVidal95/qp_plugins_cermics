subroutine get_X_proj_tan(phi_d_ortho, phi_s_ortho, psi_d_ortho, psi_s_ortho,X_mat)
 implicit none
 BEGIN_DOC
! X matrix for projection onto tangeant space at (phi_d_ortho, phi_s_ortho) and for a general matrix (psi_d_ortho, psi_s_ortho)
 END_DOC
 double precision, intent(in)  :: phi_d_ortho(ao_num, n_d_occ), phi_s_ortho(ao_num, n_s_occ)
 double precision, intent(in)  :: psi_d_ortho(ao_num, n_d_occ), psi_s_ortho(ao_num, n_s_occ)
 double precision, intent(out) :: X_mat(n_s_occ,n_d_occ)
 double precision, allocatable :: phi_s_transpose(:,:),psi_s_transpose(:,:), X_mat_tmp(:,:)
 allocate(phi_s_transpose(n_s_occ, ao_num),psi_s_transpose(n_s_occ, ao_num), X_mat_tmp(n_s_occ,n_d_occ))
 call dtranspose(phi_s_ortho,ao_num,phi_s_transpose,n_s_occ,ao_num,n_s_occ)
 call dtranspose(psi_s_ortho,ao_num,psi_s_transpose,n_s_occ,ao_num,n_s_occ)
 X_mat = 0.d0
 call sym_rect_mat_mul(phi_s_transpose,psi_d_ortho,n_s_occ,ao_num,n_d_occ,X_mat)
 call sym_rect_mat_mul(psi_s_transpose,phi_d_ortho,n_s_occ,ao_num,n_d_occ,X_mat_tmp)
 X_mat += X_mat_tmp

end


 BEGIN_PROVIDER [ double precision, X_mat_prov, (n_s_occ,n_d_occ)]
&BEGIN_PROVIDER [ double precision, norm_X_mat_prov]
 implicit none
 call get_X_proj_tan(mo_coef_d_ortho, mo_coef_s_ortho, grad_d_ortho, grad_s_ortho,X_mat_prov)
 integer :: i,j
 norm_X_mat_prov = 0.d0
 do i = 1, n_d_occ
  do j = 1, n_s_occ
   norm_X_mat_prov += X_mat_prov(j,i)**2
  enddo
 enddo
END_PROVIDER 
