use map_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !                      H,  J and K                       ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER[double precision, H_core_ortho, (ao_num,ao_num)]
  implicit none
  integer                                    :: i,j,k,l
  double precision                           :: accu
  double precision, dimension(ao_num,ao_num) :: H_core

  BEGIN_DOC
  ! Provides S^{-1/2}*H_core*S^{-1/2}
  END_DOC

  do j=1,ao_num
     do i=1,ao_num
        accu=0.d0
        do l=1,ao_num
           do k=1,ao_num
              accu += inv_sqrt_overlap(i,k)*ao_one_e_integrals(k,l)*inv_sqrt_overlap(j,l)
           enddo
        enddo
        H_core_ortho(i,j) = accu
     enddo
  enddo
  
END_PROVIDER

  
subroutine get_J_K(P,J_P,K_P)
  use map_module
  implicit none
  
  integer                                                 :: i,j,k,l
  double precision                                        :: accu_J,accu_K, get_ao_two_e_integral
  double precision, dimension(ao_num,ao_num), intent(in)  :: P
  double precision, dimension(ao_num,ao_num), intent(out) :: J_P,K_P
  
  BEGIN_DOC
  ! Compute matrix J(P) and K(P)  s.t. [J(P)]_{i,j} = \sum\limits_{k,l} (ij|kl)P_{kl}
  !                                    [K(P)]_{i,j} = \sum\limits_{k,l} (il|kj)P_{kl}
  ! for a given projector P
  END_DOC

  do j=1,ao_num
     do i=1,ao_num
        accu_J = 0.d0
        accu_K = 0.d0
        do l = 1,ao_num
           do k = 1,ao_num
              ! Beware : (ij|kl) = get_ao_two_e_integral(i,k,j,l,ao_integrals_map)
              accu_J += get_ao_two_e_integral(i,k,j,l,ao_integrals_map)*P(k,l)
              accu_K += get_ao_two_e_integral(i,k,l,j,ao_integrals_map)*P(k,l)
           enddo
        enddo
        J_P(i,j) = accu_J
        K_P(i,j) = accu_K
     enddo
  enddo

end subroutine get_J_K


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !                       ENERGY                           ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  double precision function get_rohf_energy(Pd,Ps) result(result)
    implicit none
    
    double precision, dimension(ao_num,ao_num), intent(in) :: Pd,Ps
    integer                                                :: i,j
    double precision                                       :: energy
    double precision, dimension(ao_num,ao_num)             :: Jd,Kd,Js,Ks, True_Pd,True_Ps       
    
    BEGIN_DOC
    ! Compute the energy
    !
    ! E(Pd,Ps) = tr[h(2*Pd+Ps)] + tr[(2J_d + K_d)(Pd+Ps)] + 1/2*tr[(J_s+K_s)Ps]
    !
    !
    END_DOC

    ! Retrieve true density matrices (ie for non orthogonal basis)
    True_Pd = matmul(inv_sqrt_overlap,matmul(Pd,inv_sqrt_overlap))
    True_Ps = matmul(inv_sqrt_overlap,matmul(Ps,inv_sqrt_overlap))

    ! Compute associated J and K
    call get_J_K(True_Pd,Jd,Kd)
    call get_J_K(True_Ps,Js,Ks)
    
    
    !Compute the trace
    energy =0.d0
    do j=1,ao_num
       do i=1,ao_num
          ! One electron energy
          energy += ao_one_e_integrals(i,j)*(2.d0*True_Pd(j,i) + True_Ps(j,i))
          !Two electron energy
          energy += (2.d0*Jd(i,j) - Kd(i,j))*(True_Pd(j,i) + True_Ps(j,i))&
               + 0.5d0*(Js(i,j) - Ks(i,j))*True_Ps(j,i)
       enddo
    enddo
    
    result = energy

  end function get_rohf_energy



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !                       GRADIENT                         ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_grad_projs(Pd,Ps,Gd,Gs)
  implicit none
  BEGIN_DOC
  !Compute gradient of the energy at point (Pd,Ps) projected on the tangent space
  ! T(Pd,Ps)
  END_DOC

  integer                                                :: i,j,k,l
  double precision                                       :: accu_d,accu_s

  double precision, dimension(ao_num,ao_num),intent(in)  :: Pd,Ps
  double precision, dimension(ao_num,ao_num),intent(out) :: Gd,Gs
  
  double precision, dimension(ao_num,ao_num)             :: True_Pd, True_Ps, Jd,Kd,Js,Ks,Id
  double precision, dimension(ao_num,ao_num)             :: Fd,Fs,Pv,tmp_mat


  !############# 1 : Preliminary computations

  ! Retrieve true density matrices (ie for non orthogonal basis)
  ! S^{-1/2}*Pd*S^{-1/2} and S^{-1/2}*Ps*S^{-1/2}
  
  True_Pd = matmul(inv_sqrt_overlap,matmul(Pd,inv_sqrt_overlap))
  True_Ps = matmul(inv_sqrt_overlap,matmul(Ps,inv_sqrt_overlap))
  !Associated J and K
  call get_J_K(True_Pd,Jd,Kd)
  call get_J_K(True_Ps,Js,Ks)

  !############## 2 : Compute Fock operators Fd,Fs with
  ! Fd = (h + 2Jd - Kd + Ks - 1/2 Ks)
  ! Fs = 1/2(h + 2Jd - Kd + Js - Ks)
  ! Beware : these are the Fock operators for the orthogonalized basis. See formula (12) and (13) in the report.
 
  tmp_mat = 2.d0*Jd - Kd + Js - 0.5d0*Ks
  Fd = H_core_ortho + matmul(inv_sqrt_overlap,matmul(tmp_mat,inv_sqrt_overlap))
  tmp_mat = Jd - 0.5d0*Kd + 0.5d0*Js - 0.5d0*Ks
  Fs = 0.5d0*H_core_ortho + matmul(inv_sqrt_overlap,matmul(tmp_mat,inv_sqrt_overlap))

  !############## 3 : Compute the projected gradients
  !Pv
  Id = 0.d0
  do i=1,ao_num
     Id(i,i)=1.d0
  enddo
  Pv = Id - Ps - Pd
  
  tmp_mat = Fd-Fs  
  Gd = matmul(Pd,matmul(tmp_mat,Ps)) + matmul(Ps,matmul(tmp_mat,Pd))&
       + 2.d0*matmul(Pd,matmul(Fd,Pv)) + 2.d0*matmul(Pv,matmul(Fd,Pd))
  Gs = -matmul(Pd,matmul(tmp_mat,Ps)) - matmul(Ps,matmul(tmp_mat,Pd))&
       + 2.d0*matmul(Pv,matmul(Fs,Ps)) + 2.d0*matmul(Ps,matmul(Fs,Pv))

end subroutine get_grad_projs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!!!!! !!! !                      RETRACTION                        ! !!! !!!!!
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine retraction(P,U,R_P)
  implicit none
  BEGIN_DOC
  ! Apply retraction R to P and retrive R_P, given eigenvectors U
  END_DOC

  integer                          :: i,j,k,N
  double precision                 :: Trace, accu

  double precision, dimension(ao_num,ao_num),intent(in)  :: P,U
  double precision, dimension(ao_num,ao_num),intent(out) :: R_P

  !First compute N
  Trace = 0.d0
  do i=1,ao_num
     Trace += P(i,i)
  enddo
  N = NINT(Trace)
  
  do i=1,ao_num
     do j=1,ao_num
        accu = 0.d0
        do k = 0,N-1
           accu += U(i,ao_num-k)*U(j,ao_num-k)
        enddo
        R_P(i,j) = accu
     enddo
  enddo

end subroutine retraction
