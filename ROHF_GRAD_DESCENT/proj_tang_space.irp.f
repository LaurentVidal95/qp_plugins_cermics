subroutine get_X_proj_tan(phi_d_ortho, phi_s_ortho, psi_d_ortho, psi_s_ortho,X_mat, X_mat_t)
 implicit none
 BEGIN_DOC
! X matrix for projection onto tangeant space at (phi_d_ortho, phi_s_ortho) and for a general matrix (psi_d_ortho, psi_s_ortho)
 END_DOC
 double precision, intent(in)  :: phi_d_ortho(ao_num, n_d_occ), phi_s_ortho(ao_num, n_s_occ)
 double precision, intent(in)  :: psi_d_ortho(ao_num, n_d_occ), psi_s_ortho(ao_num, n_s_occ)
 double precision, intent(out) :: X_mat(n_s_occ,n_d_occ), X_mat_t(n_d_occ,n_s_occ)
 double precision, allocatable :: phi_s_transpose(:,:),psi_s_transpose(:,:)
 allocate(phi_s_transpose(n_s_occ, ao_num),psi_s_transpose(n_s_occ, ao_num))
 call dtranspose(phi_s_ortho,ao_num,phi_s_transpose,n_s_occ,ao_num,n_s_occ)
 call dtranspose(psi_s_ortho,ao_num,psi_s_transpose,n_s_occ,ao_num,n_s_occ)
 
 X_mat = 0.d0
 call sym_rect_mat_mul(phi_s_transpose,psi_d_ortho,n_s_occ,ao_num,n_d_occ,X_mat)
 call sym_rect_mat_mul(psi_s_transpose,phi_d_ortho,n_s_occ,ao_num,n_d_occ,X_mat)

 call dtranspose(X_mat,n_s_occ,X_mat_t,n_d_occ,n_s_occ,n_d_occ)

end

subroutine get_tang_space(phi_d_ortho, phi_s_ortho, psi_d_ortho, psi_s_ortho, proj_d, proj_s)
 implicit none
 BEGIN_DOC
! projection onto tangeant space at (phi_d_ortho, phi_s_ortho) and for a general matrix (psi_d_ortho, psi_s_ortho)
 END_DOC
 double precision, intent(in)  :: phi_d_ortho(ao_num, n_d_occ), phi_s_ortho(ao_num, n_s_occ)
 double precision, intent(in)  :: psi_d_ortho(ao_num, n_d_occ), psi_s_ortho(ao_num, n_s_occ)
 double precision, intent(out) :: proj_d(ao_num,n_d_occ), proj_s(ao_num,n_s_occ)

 double precision :: X_mat(n_s_occ,n_d_occ), X_mat_t(n_d_occ,n_s_occ)
 double precision, allocatable :: proj_d_tmp(:,:),proj_s_tmp(:,:),phi_d_phi_d_t(:,:)
 integer :: i
 allocate(proj_d_tmp(ao_num, n_d_occ),proj_s_tmp(ao_num, n_s_occ),phi_d_phi_d_t(ao_num,ao_num))
 proj_s = 0.d0
 call get_X_proj_tan(phi_d_ortho, phi_s_ortho, psi_d_ortho, psi_s_ortho,X_mat, X_mat_t)
 print*,'X_mat'
 do i = 1, n_s_occ
  write(*,'(100(F16.10,X))')X_mat(i,:)
 enddo
 proj_d_tmp = 0.d0
 call sym_rect_mat_mul(phi_s_ortho,X_mat,ao_num,n_s_occ,n_d_occ,proj_d_tmp)
 proj_d_tmp *= -0.5d0
 print*,'proj_d_tmp'
 do i = 1, ao_num
  write(*,'(100(F16.10,X))')proj_d_tmp(i,:)
 enddo
 call phi_phi_t(phi_d_ortho,ao_num,n_d_occ,phi_d_phi_d_t)
 print*,'phi_d_phi_d_t'
 do i = 1, ao_num
  write(*,'(100(F16.10,X))')phi_d_phi_d_t(i,:)
 enddo
 phi_d_phi_d_t -=  Id_ao_num
 phi_d_phi_d_t *= -1.d0
 print*,'phi_d_phi_d_t'
 do i = 1, ao_num
  write(*,'(100(F16.10,X))')phi_d_phi_d_t(i,:)
 enddo
 proj_d = 0.d0
 call sym_rect_mat_mul(phi_d_phi_d_t,psi_d_ortho,ao_num,ao_num,n_d_occ,proj_d)
 proj_d += proj_d_tmp
end


 BEGIN_PROVIDER [ double precision, X_mat_prov  , (n_s_occ,n_d_occ)]
&BEGIN_PROVIDER [ double precision, X_mat_t_prov, (n_d_occ,n_s_occ)]
 implicit none
 call get_X_proj_tan(mo_coef_d_ortho, mo_coef_s_ortho, grad_d_ortho, grad_s_ortho,X_mat_prov,X_mat_t_prov)
 integer :: i
 print*,'X_mat_prov'
 do i = 1, n_s_occ
  write(*,'(100(F16.10,X))')X_mat_prov(i,:)
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, psi_d_generic, (ao_num, n_d_occ)]
&BEGIN_PROVIDER [ double precision, psi_s_generic, (ao_num, n_s_occ)]
 implicit none
 integer :: i,j
 do i = 1, n_d_occ
  do j = 1, ao_num
   psi_d_generic(j,i) = 1.d0
  enddo
 enddo

 do i = 1, n_s_occ
  do j = 1, ao_num
   psi_s_generic(j,i) = 1.d0
  enddo
 enddo
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, X_mat_generic  , (n_s_occ,n_d_occ)]
&BEGIN_PROVIDER [ double precision, X_mat_t_generic, (n_d_occ,n_s_occ)]
 implicit none
 call get_X_proj_tan(mo_coef_d_ortho, mo_coef_s_ortho, psi_d_generic, psi_s_generic,X_mat_generic,X_mat_t_generic)
 integer :: i
 print*,'X_mat_generic'
 do i = 1, n_s_occ
  write(*,'(100(F16.10,X))')X_mat_generic(i,:)
 enddo
END_PROVIDER _

 BEGIN_PROVIDER [double precision, proj_d_generic, (ao_num,n_d_occ)]
&BEGIN_PROVIDER [double precision, proj_s_generic, (ao_num,n_s_occ)]
 implicit none
 double precision, allocatable :: proj_d_tmp(:,:),proj_s_tmp(:,:),phi_d_phi_d_t(:,:)
 allocate(proj_d_tmp(ao_num, n_d_occ),proj_s_tmp(ao_num, n_s_occ),phi_d_phi_d_t(ao_num,ao_num))
 proj_d_generic = 0.d0
 proj_d_tmp = 0.d0
 call sym_rect_mat_mul(mo_coef_s_ortho,X_mat_generic,ao_num,n_s_occ,n_d_occ,proj_d_tmp)
 proj_d_tmp *= -0.5d0
 phi_d_phi_d_t = Id_ao_num - density_mat_d_ortho
 call sym_rect_mat_mul(phi_d_phi_d_t,psi_d_generic,ao_num,ao_num,n_d_occ,proj_d_generic)
 proj_d_generic += proj_d_tmp

 proj_s_generic = 0.d0
 proj_s_tmp = 0.d0
 call sym_rect_mat_mul(mo_coef_d_ortho,X_mat_t_generic,ao_num,n_d_occ,n_s_occ,proj_s_tmp)
 proj_s_tmp *= -0.5d0
 phi_d_phi_d_t = Id_ao_num - density_mat_s_ortho
 call sym_rect_mat_mul(phi_d_phi_d_t,psi_s_generic,ao_num,ao_num,n_s_occ,proj_s_generic)
 proj_s_generic += proj_s_tmp
END_PROVIDER 

 BEGIN_PROVIDER [double precision, norm_proj_d_generic]
&BEGIN_PROVIDER [double precision, norm_proj_s_generic]
 implicit none
 integer :: i,j
 norm_proj_d_generic = 0.d0
 do i = 1, n_d_occ
  do j = 1, ao_num
  norm_proj_d_generic += proj_d_generic(j,i)**2 
  enddo
 enddo

 norm_proj_s_generic = 0.d0
 do i = 1, n_s_occ
  do j = 1, ao_num
  norm_proj_s_generic += proj_s_generic(j,i)**2 
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, proj_d_prov, (ao_num,n_d_occ)]
&BEGIN_PROVIDER [double precision, proj_s_prov, (ao_num,n_s_occ)]
 implicit none
 double precision, allocatable :: proj_d_tmp(:,:),proj_s_tmp(:,:),phi_d_phi_d_t(:,:)
 allocate(proj_d_tmp(ao_num, n_d_occ),proj_s_tmp(ao_num, n_s_occ),phi_d_phi_d_t(ao_num,ao_num))
 proj_d_prov = 0.d0
 proj_d_tmp = 0.d0
 call sym_rect_mat_mul(mo_coef_s_ortho,X_mat_prov,ao_num,n_s_occ,n_d_occ,proj_d_tmp)
 proj_d_tmp *= -0.5d0
 phi_d_phi_d_t = Id_ao_num - density_mat_d_ortho
 call sym_rect_mat_mul(phi_d_phi_d_t,grad_d_ortho,ao_num,ao_num,n_d_occ,proj_d_prov)
 proj_d_prov += proj_d_tmp

 proj_s_prov = 0.d0
 proj_s_tmp = 0.d0
 call sym_rect_mat_mul(mo_coef_d_ortho,X_mat_t_prov,ao_num,n_d_occ,n_s_occ,proj_s_tmp)
 proj_s_tmp *= -0.5d0
 phi_d_phi_d_t = Id_ao_num - density_mat_s_ortho
 call sym_rect_mat_mul(phi_d_phi_d_t,grad_s_ortho,ao_num,ao_num,n_s_occ,proj_s_prov)
 proj_s_prov += proj_s_tmp
END_PROVIDER 

 BEGIN_PROVIDER [double precision, norm_proj_d_prov]
&BEGIN_PROVIDER [double precision, norm_proj_s_prov]
 implicit none
 integer :: i,j
 norm_proj_d_prov = 0.d0
 do i = 1, n_d_occ
  do j = 1, ao_num
  norm_proj_d_prov += proj_d_prov(j,i)**2 
  enddo
 enddo

 norm_proj_s_prov = 0.d0
 do i = 1, n_s_occ
  do j = 1, ao_num
  norm_proj_s_prov += proj_s_prov(j,i)**2 
  enddo
 enddo
 END_PROVIDER 

 BEGIN_PROVIDER [ double precision, norm_X_mat_generic]
&BEGIN_PROVIDER [ double precision, norm_X_mat_t_generic]
 implicit none
 integer :: i,j
 norm_X_mat_generic = 0.d0
 norm_X_mat_t_generic = 0.d0
 do i = 1, n_d_occ
  do j = 1, n_s_occ
   norm_X_mat_generic += X_mat_generic(j,i)**2
   norm_X_mat_t_generic += X_mat_t_generic(i,j)**2
  enddo
 enddo
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, norm_X_mat_prov]
&BEGIN_PROVIDER [ double precision, norm_X_mat_t_prov]
 implicit none
 integer :: i,j
 norm_X_mat_prov = 0.d0
 norm_X_mat_t_prov = 0.d0
 do i = 1, n_d_occ
  do j = 1, n_s_occ
   norm_X_mat_prov += X_mat_prov(j,i)**2
   norm_X_mat_t_prov += X_mat_t_prov(i,j)**2
  enddo
 enddo
END_PROVIDER 

