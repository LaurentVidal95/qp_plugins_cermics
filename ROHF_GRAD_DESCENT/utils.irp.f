


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

double precision function norm_of_mat(a,n)
 implicit none
 integer, intent(in) :: n
 double precision, intent(in) :: a(n,n)
 integer :: i,j
 norm_of_mat = 0.d0
 do i = 1, n
  do j = 1, n
   norm_of_mat += dabs(a(j,i))
  enddo
 enddo

end


subroutine sym_mat(A,n,A_sym)
 implicit none
 double precision, intent(in) :: A(n,n)
 double precision, intent(out):: A_sym(n,n)
 integer, intent(in) :: n
 integer :: i,j
 do i = 1, n
  do j = 1, n
   A_sym(j,i) = (A(j,i) + A(i,j)) * 0.5d0 
  enddo
 enddo
end

subroutine project_tangent_space(Pd,Ps, Pv, M, N, M_proj, N_proj)
 implicit none
 BEGIN_DOC
! You have two symetric general matrix M and N, and you want to project them on the 
!
! hyperplane tangent to the manifold in (Pd, Ps)
 END_DOC
 double precision, intent(in) :: Pd(ao_num, ao_num), Ps(ao_num, ao_num), Pv(ao_num, ao_num)
 double precision, intent(in) :: M(ao_num, ao_num), N(ao_num, ao_num)
 double precision, intent(out):: M_proj(ao_num, ao_num), N_proj(ao_num, ao_num)
 double precision :: M_minus_N(ao_num, ao_num)
 double precision :: mat_tmp(ao_num, ao_num), mapt_tmp2(ao_num, ao_num)
 M_minus_N = M - N

 ! M_proj = Pd(M-N)Ps
 call sym_square_mat_mul(Pd,M_minus_N,ao_num,mat_tmp)
 call sym_square_mat_mul(mat_tmp,Ps,ao_num,M_proj)

 ! mat_tmp2 = 2 Pd M Pv
 call sym_square_mat_mul(Pd,M,ao_num, mat_tmp)
 call sym_square_mat_mul(map_tmp, Pv, ao_num, mat_tmp2)
 mat_tmp2 = 2.d0 * mat_tmp2

 ! M_proj = sym(Pd(M-N)Ps + 2 Pd M Pv)
 mat_tmp = M_proj + mat_tmp2
 call sym_mat(mat_tmp,ao_num,M_proj)

 ! N_proj = Pd(N-M) Ps
 call sym_square_mat_mul(Pd,M_minus_N,ao_num,mat_tmp)
 mat_tmp = - mat_tmp
 call sym_square_mat_mul(mat_tmp,Ps,ao_num,N_proj)

 ! mat_tmp2 = 2 Ps N Pv
 call sym_square_mat_mul(Ps,N,ao_num, mat_tmp)
 call sym_square_mat_mul(map_tmp, Pv, ao_num, mat_tmp2)
 mat_tmp2 = 2.d0 * mat_tmp2

 ! N_proj = sym(Pd(N-M)Ps + 2 Ps M Pv)
 mat_tmp = N_proj + mat_tmp2
 call sym_mat(mat_tmp,ao_num,N_proj)

  
end

subroutine sym_square_mat_mul(A,B,n,C)
 implicit none
 BEGIN_DOC
! C = A * B for symmetric square matrix A, B
 END_DOC
 double precision, intent(in) :: A(n, n), B(n, n)
 integer, intent(in)  :: n
 double precision, intent(out) :: C(n,n)

 call dgemm('N','N',n,n,n,1.d0, &
      A, size(A,1), &
      B, size(B,1), 0.d0, &
      C, size(C,1))
end