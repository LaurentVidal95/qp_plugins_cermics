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
