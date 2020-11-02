!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !        OVERLAP, SQRT_OVERLAP, INV_SQRT_OVERLAP         ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER[double precision, overlap_matrix, (ao_num,ao_num)]
  implicit none
  integer   :: i,j
  
  BEGIN_DOC
  ! Simply provides the overlap matrix S
  END_DOC

  do i=1,ao_num
     do j=1,ao_num
        overlap_matrix(i,j) = ao_overlap(i,j)
     enddo
  enddo
  
  END_PROVIDER

  
  BEGIN_PROVIDER[double precision,  sqrt_overlap, (ao_num, ao_num)]
&BEGIN_PROVIDER[double precision,  inv_sqrt_overlap, (ao_num,ao_num)]   
  implicit none
  BEGIN_DOC
  !Provides sqare root and inverse square root of S.
  END_DOC

  integer                                    :: i,j,k,l
  double precision, dimension(ao_num,ao_num) :: eigvectors
  double precision, dimension(ao_num)        :: eigvalues
  double precision                           :: accu, accu_inv

  !Diagonalize S
  call lapack_diagd(eigvalues,eigvectors,overlap_matrix,ao_num,ao_num)

  !Compute the sqare root and inverse square root
  do j=1,ao_num
     do i=1,ao_num
        accu = 0.d0
        accu_inv = 0.d0
        do k =1,ao_num
           accu += eigvectors(i,k)*dsqrt(eigvalues(k))*eigvectors(j,k)
           accu_inv += eigvectors(i,k)*(1/dsqrt(eigvalues(k)))*eigvectors(j,k)
        enddo
        sqrt_overlap(i,j) = accu
        inv_sqrt_overlap(i,j) = accu_inv
     enddo
  enddo

  END_PROVIDER

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !                    INITIAL DENSITIES                   ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
subroutine init_guess(Pd_init,Ps_init)
    
  implicit none
  BEGIN_DOC
  !Provides initial orthogonal projectors on doubly occupied MOs
  !and singly occupied MOs based on the guess provided by QP2.
  !
  ! We recall that :
  ! Pd_init = S^{1/2}*(C_dC_d^T)*S^{1/2}
  ! Ps_init = S^{1/2}*(C_sC_s^T)*S^{1/2}
  !
  ! where Cd and Cs are the coefficient on the AO basis of doubly occupied and singly occupied
  ! MOs respectively.
  END_DOC
  integer                                                 :: i,j,k,l,m
  double precision                                        :: accu_d,accu_s
  double precision, dimension(ao_num,ao_num)              :: Pd_tmp,Ps_tmp
  double precision, dimension(ao_num,ao_num), intent(out) :: Pd_init,Ps_init

  double precision, dimension(ao_num,ao_num)              :: H,tmp_mat,eigvectors
  double precision, dimension(ao_num)                     :: eigvalues

  !Compute Pd Ps for non-orthogonal basis given by mo_coefs

  call huckel_guess()
  
  Pd_tmp = 0.d0
  Ps_tmp = 0.d0

  do j=1,ao_num
     do i=1,ao_num
        accu_d = 0.d0
        accu_s = 0.d0
        do k=1,elec_beta_num
           accu_d += mo_coef_transp(k,i)*mo_coef_transp(k,j)
        enddo
        do k=elec_beta_num+1,elec_alpha_num
           accu_s +=  mo_coef_transp(k,i)*mo_coef_transp(k,j)
        enddo
        Pd_tmp(i,j) = accu_d
        Ps_tmp(i,j) = accu_s
     enddo
  enddo

  ! ! ###### OPTIONAL : Import guess written in files "Pd_core.dat" and "Ps_core.dat"
  ! open(1,file = 'Pd_core.dat')
  ! read(1,*) Pd_tmp
  ! close(1)

  ! open(2,file = 'Ps_core.dat')
  ! read(2,*) Ps_tmp
  ! close(2)
  ! ! ################

  !Orthogonalization
  Pd_init = matmul(sqrt_overlap,matmul(Pd_tmp,sqrt_overlap))
  Ps_init = matmul(sqrt_overlap,matmul(Ps_tmp,sqrt_overlap))


  
        
end subroutine init_guess



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !        TEST THE PROPERTIES OF PROJECTORS               ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


double precision function test_projs(Pd,Ps) result(result)
  implicit none
  BEGIN_DOC
  ! Test the properties of Pd and Ps. The test must be very close to 0.
  END_DOC

  integer                                     :: i,j,k,l
  double precision                            :: accu_d,accu_s, out
  double precision, dimension(ao_num,ao_num)  :: temp_matrix

  double precision, dimension(ao_num,ao_num), intent(in) :: Pd,Ps

  out = 0.d0

  ! Norm of Pd*Pd - Pd
  temp_matrix = matmul(Pd,Pd) - Pd
  out += Norm2(temp_matrix)
  
  ! Norm of Ps*Ps - Ps
  temp_matrix= matmul(Ps,Ps) - Ps
  out += Norm2(temp_matrix)

  !Norm of Ps*Pd
  temp_matrix = matmul(Pd,Ps)
  out  += Norm2(temp_matrix)

  !Test traces
  accu_d = 0.d0
  accu_s = 0.d0
  do i=1,ao_num
     accu_d += Pd(i,i)
  enddo
  
  do i=1,ao_num
     accu_s += Ps(i,i)
  enddo
  out += accu_d - elec_beta_num
  out -= accu_s - (elec_alpha_num-elec_beta_num)

  result = out
  

end function test_projs

  
