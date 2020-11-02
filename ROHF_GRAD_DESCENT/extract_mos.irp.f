! subroutine extract_occupied_mos(Pd,Ps,orbitals)
!   implicit none
!   BEGIN_DOC
!   ! Diagonalize Pd and Ps and extract MOs with non zero occupation number
!   !
!   ! Do we extract from True densities or Pd Ps for orthogonal basis ?
!   END_DOC


!   integer                                                  :: i,j,k,l
!   double precision, dimension(ao_num,ao_num), intent(in)   :: Pd,Ps
!   double precision, dimension(ao_num,elec_alpha_num), intent(out)  :: orbitals
!   double precision, allocatable                            :: eigvalues(:),eigvectors(:,:),orbitals_tmp(:,:)

!   double precision, dimension(ao_num,ao_num)               :: True_Pd,True_Ps
  
!   allocate(eigvectors(ao_num,ao_num), eigvalues(ao_num),orbitals_tmp(ao_num,ao_num))

!   ! Retrive true density matrices 
!   True_Pd = matmul(inv_sqrt_overlap,matmul(Pd,inv_sqrt_overlap))
!   True_Ps = matmul(inv_sqrt_overlap,matmul(Ps,inv_sqrt_overlap))

!   orbitals = 0.d0

!   ! Extract Mos from Pd
!   ! "-True_Pd" in order to have the eigvectors in decreasing eigvalues order.
!   call lapack_diagd(eigvalues,eigvectors,-True_Pd,ao_num,ao_num)
!   k = 0
!   do i = 1,ao_num
!      if (NINT(eigvalues(i)) == -1) then
!         k+=1
!         do j = 1,ao_num
!            orbitals(j,i) = eigvectors(j,i)
!         enddo
!      endif
!   enddo

  
!   !Extract Mos from Ps
!   call lapack_diagd(eigvalues,eigvectors,-True_Ps,ao_num,ao_num)
!   do i = 1,ao_num
!      if (NINT(eigvalues(i)) == -1) then
!         k+=1
!         do j = 1,ao_num
!            orbitals(j,k) = eigvectors(j,i)
!         enddo
!      endif
!   enddo

!   do j = elec_alpha_num+1,ao_num
!      do i = 1,ao_num
!         ! orbitals(j,j) = 1.d0
!         orbitals(i,j) = eigvectors(i,j)
!      enddo
!   enddo
  
! end subroutine extract_occupied_mos


! subroutine Gramm_Schmidt_ortho(orbitals)
!   implicit none
!   BEGIN_DOC
!   ! Generate virtual orbitals by Gramm-Schmidt orthogonalization and normalize
!   ! d, s et v orbitals.
!   !
!   ! Here, "orbitals_tmp" contains the orthogonal orbitals still to be normalized.
!   ! Let e_i be one AO. Then the orthogonalized orbital e_prime_i is obtained with :
!   !
!   ! for any AO e_k,
!   !
!   ! <e_k,e_prime_i> = S_{k,i}
!   !                   - \sum_{j<i}\sum_{mu,nu = 1,ao_num} C_{mu,j}C_{nu,j}S_{mu,i}S_{nu,k}
!   !
!   END_DOC

!   integer                                                   :: i,j,k,mu,nu
!   double precision, dimension(ao_num,ao_num), intent(inout) :: orbitals
!   double precision                                          :: accu
  
  
!   ! Gramm-schmidt to generate virtual orbitals from AOs
!   do i = elec_alpha_num+1,ao_num ! i = virtual orbital to be orthogonalized
!      do k = 1,ao_num ! k runs through AOs
!         accu = ao_overlap(i,k)
!         do j = 1,i-1 !j runs through already orthogonal orbitals
!            ! add \sum_{mu,nu} C_{nu,j}C_{mu,j} S_{mu,j}S{nu,k} for every j
!            do mu = 1,ao_num
!               do nu = 1,ao_num
!                  accu += -orbitals_tmp(mu,j)*orbitals_tmp(nu,j)*ao_overlap(mu,i)*ao_overlap(nu,k)
!               enddo
!            enddo
!         enddo
!         orbitals_tmp(k,i) = accu
!      enddo
!   enddo

!   orbitals = orbitals_tmp


  
  !Normalize orbitals
  ! do i=1,ao_num
  !    accu = 0.d0
  !    do j=1,ao_num
  !       do k=1,ao_num
           
  !          accu += 
  !       enddo
  !    enddo
     
  ! enddo
  
! end subroutine Gramm_Schmidt_ortho
