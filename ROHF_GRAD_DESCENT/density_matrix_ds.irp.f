!BEGIN_PROVIDER [ double precision, proj_basis, (ao_num, ao_num)]
! implicit none
! call dgemm('N','T',ao_num,ao_num,ao_num,1.d0, &                                                                                   
!      mo_coef, size(mo_coef,1), &
!      mo_coef, size(mo_coef,1), 0.d0, &
!      proj_basis, size(proj_basis,1))
!END_PROVIDER 

!BEGIN_PROVIDER [ double precision, proj_basis_ortho, (ao_num, ao_num) ]
! implicit none
! call ao_to_ao_ortho(proj_basis, proj_basis_ortho)
!END_PROVIDER 
!
!BEGIN_PROVIDER [ double precision, proj_virt, (ao_num, ao_num) ]
! implicit none
! proj_virt = proj_basis - density_mat_d - density_mat_s 
!END_PROVIDER 
!
!BEGIN_PROVIDER [ double precision, proj_virt_ortho, (ao_num, ao_num) ]
! implicit none
! call ao_to_ao_ortho(proj_virt, proj_virt_ortho)
!END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_s_ortho, (ao_num, ao_num)]
 implicit none
 call ao_to_ao_ortho(density_mat_s, density_mat_s_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_d_ortho, (ao_num, ao_num)]
 implicit none
 call ao_to_ao_ortho(density_mat_d, density_mat_d_ortho)
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

BEGIN_PROVIDER [ double precision, delta_mo_d, (ao_num, elec_beta_num)]
 implicit none
 BEGIN_DOC
! epsilonesque variation of the doubly occupied MOs
 END_DOC
 integer :: i,j
 delta_mo_d = 0.d0
 do i = 1, elec_beta_num
  do j = 1, ao_num
   delta_mo_d(j,i) = epsilon_delta_mo
  enddo
 enddo
 delta_mo_d = delta_mo_d / dble(elec_beta_num * ao_num)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, delta_mo_s, (ao_num, n_single_occupied)]
 implicit none
 BEGIN_DOC
! epsilonesque variation of the singly occupied MOs
 END_DOC
 integer :: i,j
 delta_mo_s = 0.d0
 do i = 1, elec_beta_num
  do j = 1, ao_num
   delta_mo_s(j,i) = epsilon_delta_mo
  enddo
 enddo
 delta_mo_s = delta_mo_s / dble(n_single_occupied * ao_num)
END_PROVIDER 

BEGIN_PROVIDER [ integer, n_single_occupied]
implicit none
 n_single_occupied = elec_alpha_num - elec_beta_num
END_PROVIDER 

