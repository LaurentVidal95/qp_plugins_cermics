!BEGIN_PROVIDER [ double precision, density_mat_s_ortho, (ao_num, ao_num)]
! implicit none
! call ao_to_ao_ortho(density_mat_s, density_mat_s_ortho)
!END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_d_ortho, (ao_num, ao_num)]
 implicit none
 density_mat_d_ortho = 0.d0
 call sym_rect_mat_mul(mo_coef_d_ortho,mo_coef_d_ortho_t,ao_num,n_d_occ,ao_num,density_mat_d_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, trace_dm_d_ortho]
 implicit none
 integer ::i 
 trace_dm_d_ortho = 0.d0
 do i = 1, ao_num
  trace_dm_d_ortho += density_mat_d_ortho(i,i)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_s_ortho, (ao_num, ao_num)]
 implicit none
 density_mat_s_ortho = 0.d0
 call sym_rect_mat_mul(mo_coef_s_ortho,mo_coef_s_ortho_t,ao_num,n_s_occ,ao_num,density_mat_s_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, trace_dm_s_ortho]
 implicit none
 integer ::i 
 trace_dm_s_ortho = 0.d0
 do i = 1, ao_num
  trace_dm_s_ortho += density_mat_s_ortho(i,i)
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_d, (ao_num, ao_num)]
 implicit none
 density_mat_d = SCF_density_matrix_ao_beta
END_PROVIDER 


BEGIN_PROVIDER [ double precision, density_mat_s, (ao_num, ao_num)]
 implicit none
 density_mat_s = SCF_density_matrix_ao_alpha - SCF_density_matrix_ao_beta
END_PROVIDER 

BEGIN_PROVIDER [ double precision, epsilon_delta_mo]
 implicit none
 BEGIN_DOC
! epsilon for variation of MOs
 END_DOC
 epsilon_delta_mo = 1.d-5
END_PROVIDER 

BEGIN_PROVIDER [ double precision, delta_mo_d, (ao_num, n_d_occ)]
 implicit none
 BEGIN_DOC
! epsilonesque variation of the doubly occupied MOs
 END_DOC
 integer :: i,j
 delta_mo_d = 0.d0
 do i = 1, n_d_occ
  do j = 1, ao_num
   delta_mo_d(j,i) = epsilon_delta_mo
  enddo
 enddo
 delta_mo_d = delta_mo_d / dble(n_d_occ * ao_num)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, delta_mo_s, (ao_num, n_s_occ)]
 implicit none
 BEGIN_DOC
! epsilonesque variation of the singly occupied MOs
 END_DOC
 integer :: i,j
 delta_mo_s = 0.d0
 do i = 1, n_d_occ
  do j = 1, ao_num
   delta_mo_s(j,i) = epsilon_delta_mo
  enddo
 enddo
 delta_mo_s = delta_mo_s / dble(n_s_occ * ao_num)
END_PROVIDER 

BEGIN_PROVIDER [ integer, n_s_occ]
 implicit none
 n_s_occ = elec_alpha_num - n_d_occ
END_PROVIDER 

BEGIN_PROVIDER [ integer, list_n_s_occ, (n_s_occ)]
 implicit none
 integer :: i,j
 do i = 1, n_s_occ
  j = i + n_d_occ 
  list_n_s_occ(i) = j
 enddo
 print*,'list_n_s_occ ',list_n_s_occ
END_PROVIDER 


BEGIN_PROVIDER [ integer, n_d_occ]
 implicit none
 n_d_occ =  elec_beta_num
END_PROVIDER 


subroutine mo_to_dm(phi_mo, n_mo, n_ao, dm_ao)
 implicit none
 BEGIN_DOC
! takes a matrix phi_mo(n_ao, n_mo) and gives back 
!
! dm_ao(n_ao, n_ao) which is the density matrix associated 
 END_DOC
 integer :: n_mo, n_ao
 double precision, intent(in) :: phi_mo(n_ao, n_mo)
 double precision, intent(out):: dm_ao(n_ao, n_ao)
 double precision, allocatable :: phi_t(:,:)
 allocate(phi_t(n_mo, n_ao))
 call dtranspose(phi_mo,n_ao,phi_t,n_mo,n_ao,n_mo)
 dm_ao = 0.d0
 call sym_rect_mat_mul(phi_mo,phi_t,n_ao,n_mo,n_ao,dm_ao)

end

BEGIN_PROVIDER [ double precision, mo_coef_d, (ao_num, n_d_occ)]
 implicit none
 integer :: i,j
 mo_coef_d = 0.d0
 do i = 1, n_d_occ
  do j = 1, ao_num
   mo_coef_d(j,i) = mo_coef(j,i)
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_coef_d_ortho, (ao_num, n_d_occ)]
 implicit none
 mo_coef_d_ortho = 0.d0
 call sym_rect_mat_mul(sqrt_overlap,mo_coef_d,ao_num,ao_num,n_d_occ,mo_coef_d_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_coef_s_ortho, (ao_num, n_s_occ)]
 implicit none
 mo_coef_s_ortho = 0.d0
 call sym_rect_mat_mul(sqrt_overlap,mo_coef_s,ao_num,ao_num,n_s_occ,mo_coef_s_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_coef_d_ortho_t, (n_d_occ, ao_num)]
 implicit none
 call dtranspose(mo_coef_d_ortho,ao_num, mo_coef_d_ortho_t,n_d_occ,ao_num,n_d_occ)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, mo_coef_s_ortho_t, (n_s_occ, ao_num)]
 implicit none
 call dtranspose(mo_coef_s_ortho,ao_num, mo_coef_s_ortho_t,n_s_occ,ao_num,n_s_occ)
END_PROVIDER 
BEGIN_PROVIDER [ double precision, mo_coef_s, (ao_num, n_s_occ)]
 implicit none
 integer :: i,j
 mo_coef_s = 0.d0
 do i = 1, n_s_occ
  do j = 1, ao_num
   mo_coef_s(j,i) = mo_coef(j,list_n_s_occ(i))
  enddo
 enddo
END_PROVIDER 


BEGIN_PROVIDER [ double precision, Id_ao_num, (ao_num, ao_num)]
 implicit none
 integer :: i
 Id_ao_num = 0.d0
 do i = 1, ao_num
  Id_ao_num(i,i) = 1.d0
 enddo

END_PROVIDER 
