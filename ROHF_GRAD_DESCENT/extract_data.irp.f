program extract_data
  implicit none
  use map_module
  
  BEGIN_DOC
  ! Extract the overlap matrix, the one-e integrals and the number of  alpha and beta orbitals.
  ! May also extract the four-index tensor but this part has to be uncommented beforehand.
  ! It is not needed to canonicalize MOs.
  END_DOC

  integer                               :: i,j,k,l
  double precision                      :: accu
  integer                               :: mod
  double precision                      :: get_ao_two_e_integral



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    Write directly in julia array format                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  ! Overlap matrix
  open(10, file = "data_dir/overlap_matrix.dat")
  do i = 1, ao_num
     write(10,'(100(F16.10,X))')overlap_matrix(i,:)
  enddo
  close(10)

  ! Nums_orbitals
  open(2, file = "data_dir/nums_orbitals.dat")
  write(2, '(3(I4,X))')ao_num,elec_beta_num,elec_alpha_num-elec_beta_num
  close(2)


  ! H_core
  double precision, allocatable :: H_core(:,:)
  allocate(H_core(ao_num,ao_num))
  do j=1,ao_num
     do i = 1,ao_num
        H_core(i,j) = ao_one_e_integrals(i,j)
     enddo
  enddo
  
  open(3, file = "data_dir/H_core.dat")
  do i = 1, ao_num
     write(3,'(100(F16.10,X))')H_core(i,:)
  enddo
  close(3)


  !Four index tensor [TIME CONSUMING]
  ! open(4, file = "data_dir/Four_index_tensor.dat")
  ! do i=1,ao_num
  !    do j=1,ao_num
  !       do k=1,ao_num
  !          do l=1,ao_num
  !             accu = get_ao_two_e_integral(i,k,j,l,ao_integrals_map)
  !             write(4,*)i,j,k,l,accu
  !          enddo
  !       enddo
  !    enddo
  ! enddo
  ! close(4)

  print*,"############ EXTRACTION DONE  ############"



  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            Write 1 element and corresponding indices  per line             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  ! !Overlap matrix
  ! open(1,file = 'data_dir/overlap_matrix.dat')
  ! do i=1,ao_num
  !    do j=1,ao_num
  !       write(1,*)i,j,ao_overlap(i,j)
  !    enddo
  ! enddo
  ! close(1)
  
  ! ! Four-index tensor
  ! open(2,file = 'data_dir/Four_index_tensor.dat')
  ! do i=1,ao_num
  !    do j=1,ao_num
  !       do k=1,ao_num
  !          do l=1,ao_num
  !             accu = get_ao_two_e_integral(i,k,j,l,ao_integrals_map)
  !             write(2,*)i,j,k,l,accu
  !          enddo
  !       enddo
  !    enddo
  ! enddo


  ! close(2)

  ! !Core hamiltonian
  ! open(3,file = 'data_dir/H_core.dat')
  ! do i=1,ao_num
  !    do j=1,ao_num
  !       accu = ao_one_e_integrals(i,j)
  !       write(3,*)i,j,accu
  !    enddo
  ! enddo

  ! close(3)

  ! !Nums Orbitals
  ! open(4,file = 'data_dir/nums_orbitals.dat')
  ! accu = INT(elec_alpha_num-elec_beta_num)
  ! write(4,*)ao_num,elec_beta_num,accu
  ! close(4)

  ! print*,"############ EXTRACTION DONE ############"

end program extract_data
