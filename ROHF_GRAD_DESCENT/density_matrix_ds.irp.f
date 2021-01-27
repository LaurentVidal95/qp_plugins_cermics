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


subroutine ao_to_ao_ortho(mat_ao, mat_ao_ortho)
 implicit none
 BEGIN_DOC
! routines that takes mat_ao in the NON ORTHO AO BASIS  
!
! and performs S^1/2 mat_ao S^1/2 to give mat_ao_ortho in the ORHTO AO BASIS
 END_DOC
 double precision, intent(in)  :: mat_ao(ao_num, ao_num)
 double precision, intent(out) :: mat_ao_ortho(ao_num, ao_num)

 double precision :: mat_ao_tmp(ao_num, ao_num)

 mat_ao_ortho = 0.d0
 mat_ao_tmp = 0.d0

 ! mat_ao_tmp = S^1/2 * mat_ao 
 call dgemm('N','N',ao_num,ao_num,ao_num,1.d0, &
      sqrt_overlap, size(sqrt_overlap,1), &
      mat_ao, size(mat_ao,1), 0.d0, &
      mat_ao_tmp, size(mat_ao_tmp,1))

 ! mat_ao_ortho = mat_ao_tmp * S^1/2
 call dgemm('N','N',ao_num,ao_num,ao_num,1.d0, &
      mat_ao_tmp, size(mat_ao_tmp,1), &
      sqrt_overlap, size(sqrt_overlap,1), 0.d0 , &
      mat_ao_ortho, size(mat_ao_ortho,1))
end

subroutine ao_ortho_to_ao(mat_ao_ortho, mat_ao )
 implicit none
 BEGIN_DOC
! routines that takes mat_ao_ortho in the ORTHO AO BASIS  
!
! and performs S^-1/2 mat_ao_ortho S^-1/2 to give mat_ao in the AO BASIS
 END_DOC
 double precision, intent(in) :: mat_ao_ortho(ao_num, ao_num)
 double precision, intent(out)  :: mat_ao(ao_num, ao_num)

 double precision :: mat_ao_tmp(ao_num, ao_num)

 mat_ao = 0.d0
 mat_ao_tmp = 0.d0

 ! mat_ao_tmp = S^-1/2 * mat_ao_ortho
 call dgemm('N','N',ao_num,ao_num,ao_num,1.d0, &
      inv_sqrt_overlap, size(inv_sqrt_overlap,1), &
      mat_ao_ortho, size(mat_ao_ortho,1), 0.d0, &
      mat_ao_tmp, size(mat_ao_tmp,1))

 ! mat_ao = mat_ao_tmp * S^-1/2
 call dgemm('N','N',ao_num,ao_num,ao_num,1.d0, &
      mat_ao_tmp, size(mat_ao_tmp,1), &
      inv_sqrt_overlap, size(inv_sqrt_overlap,1), 0.d0 , &
      mat_ao, size(mat_ao,1))
end

