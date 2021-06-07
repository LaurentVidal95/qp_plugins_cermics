subroutine get_X_r(phi_d_ortho, psi_s_ortho, X_r)
 implicit none
 double precision, intent(in) :: phi_d_ortho(ao_num, n_d_occ), psi_s_ortho(ao_num, n_s_occ)
 double precision, intent(out):: X_r(n_d_occ, n_s_occ)
 double precision, allocatable :: phi_d_t(:,:) 
 allocate(phi_d_t(n_d_occ, ao_num))
 call dtranspose(phi_d_ortho,ao_num, phi_d_t,n_d_occ,ao_num,n_d_occ)

 call sym_rect_mat_mul(phi_d_t,psi_s_ortho,n_d_occ,ao_num,n_s_occ,X_r)
 X_r *= -1.d0
end

subroutine get_M_r(X_r,M_r)
 implicit none
 double precision, intent(in) :: X_r(n_d_occ, n_s_occ)
 double precision, intent(out):: M_r(n_occ_prov, n_occ_prov)
 M_r = 0.d0
 integer :: i,j
 ! M_r =  | 0         -X_r |
 !        | X_r^t      0   |
 do j = 1, n_d_occ
  do i = n_d_occ + 1, n_occ_prov
   M_r(i,j) = X_r(j,i-n_d_occ)
  enddo
 enddo

 do j = n_d_occ + 1, n_occ_prov
  do i = 1, n_d_occ
   M_r(i,j) = -1.d0 * X_r(i,j-n_d_occ)
  enddo
 enddo

end

subroutine get_prov_virt_psi(psi_d, psi_s, proj_psi_s_psi_d)
 implicit none
 double precision, intent(in) :: psi_d(ao_num, n_d_occ), psi_s(ao_num, n_s_occ)
 double precision, intent(out):: proj_psi_s_psi_d(ao_num, n_occ_prov)
 double precision, allocatable :: tmp_d(:,:), tmp_s(:,:)
 integer :: i,j
 allocate(tmp_d(ao_num, n_d_occ), tmp_s(ao_num,n_s_occ))
 tmp_d = 0.d0
 call sym_rect_mat_mul(proj_virt_prov,psi_d,ao_num,ao_num,n_d_occ,tmp_d)
 tmp_s = 0.d0
 call sym_rect_mat_mul(proj_virt_prov,psi_s,ao_num,ao_num,n_s_occ,tmp_s)
 proj_psi_s_psi_d = 0.d0
 call concat_mat(tmp_d,tmp_s,proj_psi_s_psi_d,n_d_occ,n_s_occ,ao_num)
end

subroutine get_svd_proj_virt(psi_d, psi_s,U,sigma,Vt)
 implicit none
 double precision, intent(in) :: psi_d(ao_num, n_d_occ), psi_s(ao_num, n_s_occ)
 double precision, intent(out):: U(ao_num, n_occ_prov), sigma(n_occ_prov), Vt(n_occ_prov,n_occ_prov)
 double precision, allocatable  :: proj_psi_s_psi_d(:,:)
 allocate(proj_psi_s_psi_d(ao_num, n_occ_prov))
 call get_prov_virt_psi(psi_d, psi_s, proj_psi_s_psi_d) 
 call svd(proj_psi_s_psi_d,ao_num,U,ao_num,sigma,Vt,n_occ_prov,ao_num,n_occ_prov)
end

BEGIN_PROVIDER [ double precision, X_r_generic, (n_d_occ, n_s_occ)]
 implicit none
 call get_X_r(mo_coef_d_ortho, psi_s_generic, X_r_generic)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, norm_X_r_generic]
 implicit none
 integer :: i,j
 norm_X_r_generic = 0.d0
 do j = 1, n_s_occ
  do i = 1, n_d_occ
   norm_X_r_generic += X_r_generic(i,j)**2
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, proj_virt_prov_generic, (ao_num, n_occ_prov)]
 implicit none
 call get_prov_virt_psi(psi_d_generic, psi_s_generic, proj_virt_prov_generic)
 integer :: i,j
 double precision :: accu
 accu = 0.d0
 do i = 1, n_occ_prov
  do j = 1, ao_num
   accu += proj_virt_prov_generic(j,i)**2
  enddo
 enddo
 print*,'proj_virt_prov_generic'
 do i = 1, ao_num
  write(*,'(100(F16.10,X))')proj_virt_prov_generic(i,:)
 enddo
 print*,'accu =',accu
END_PROVIDER 

 BEGIN_PROVIDER [ double precision, U_proj_virt_generic,(ao_num, n_occ_prov) ]
&BEGIN_PROVIDER [ double precision, sigma_generic, (n_occ_prov) ]
&BEGIN_PROVIDER [ double precision, Vt_generic, (n_occ_prov,n_occ_prov) ]
 implicit none
 call get_svd_proj_virt(psi_d_generic, psi_s_generic,U_proj_virt_generic,sigma_generic,Vt_generic)
 integer :: i,j
 double precision :: accu
 accu = 0.d0
 do i = 1, n_occ_prov
  accu += sigma_generic(i)
 enddo
 print*,'accu sigma = ',accu
 accu = 0.d0
 do j = 1, n_occ_prov
  do i = 1, ao_num
   accu += U_proj_virt_generic(i,j)**2
  enddo
 enddo
 print*,'accu U     = ',accu

 accu = 0.d0
 do j = 1, n_occ_prov
  do i = 1, n_occ_prov
   accu += Vt_generic(i,j)**2
  enddo
 enddo
 print*,'accu Vt    = ',accu
END_PROVIDER 

BEGIN_PROVIDER [ double precision, M_r_generic, (n_occ_prov, n_occ_prov)]
 implicit none
 call get_M_r(X_r_generic,M_r_generic)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, norm_M_r_generic]
 implicit none
 integer :: i,j
 norm_M_r_generic = 0.d0
 do j = 1, n_occ_prov
  do i = 1, n_occ_prov 
   norm_M_r_generic += M_r_generic(i,j)**2
  enddo
 enddo
END_PROVIDER 
