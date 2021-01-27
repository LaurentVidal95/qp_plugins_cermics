BEGIN_PROVIDER [ double precision, proj_basis, (ao_num, ao_num)]
 implicit none
 call dgemm('N','T',ao_num,ao_num,ao_num,1.d0, &                                                                                   
      mo_coef, size(mo_coef,1), &
      mo_coef, size(mo_coef,1), 0.d0, &
      proj_basis, size(proj_basis,1))
END_PROVIDER 

BEGIN_PROVIDER [ double precision, proj_basis_ortho, (ao_num, ao_num) ]
 implicit none
 call ao_to_ao_ortho(proj_basis, proj_basis_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, proj_virt, (ao_num, ao_num) ]
 implicit none
 proj_virt = proj_basis - density_mat_d - density_mat_s 
END_PROVIDER 

BEGIN_PROVIDER [ double precision, proj_virt_ortho, (ao_num, ao_num) ]
 implicit none
 call ao_to_ao_ortho(proj_virt, proj_virt_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_d, (ao_num, ao_num)]
 implicit none
 density_mat_d = SCF_density_matrix_ao_beta
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_d_ortho, (ao_num, ao_num)]
 implicit none
 call ao_to_ao_ortho(density_mat_d, density_mat_d_ortho)
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_s, (ao_num, ao_num)]
 implicit none
 density_mat_s = SCF_density_matrix_ao_alpha - SCF_density_matrix_ao_beta
END_PROVIDER 

BEGIN_PROVIDER [ double precision, density_mat_s_ortho, (ao_num, ao_num)]
 implicit none
 call ao_to_ao_ortho(density_mat_s, density_mat_s_ortho)
END_PROVIDER 
