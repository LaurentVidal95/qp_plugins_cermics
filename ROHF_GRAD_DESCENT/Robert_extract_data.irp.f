program Robert_extract_data
  implicit none
  use map_module
  ! use huckel_guess
  ! use test_projs
  
  BEGIN_DOC
  ! Extract the overlap matrix, the four-index tensor, the one-e integrals and the number of
  ! alpha and beta orbitals
  END_DOC

  integer                               :: i,j,k,l
  double precision                      :: accu
  integer                               :: mod
  double precision                      :: get_ao_two_e_integral

  double precision, dimension(ao_num,ao_num) :: Pd_tmp, Ps_tmp
  !double precision, allocatable    :: Pd_tmp(:,:),Ps_tmp(:,:)

  double precision                      :: accu_d,accu_s

  double precision                       :: test_projs
  
  !####################################################### Write 1 element and corresponding indices  per line

  ! allocate(Pd_tmp(ao_num,ao_num),Ps_tmp(ao_num,ao_num))

  
  !Overlap matrix
  open(1,file = 'overlap_matrix.dat')
  do i=1,ao_num
     do j=1,ao_num
        write(1,*)i,j,ao_overlap(i,j)
     enddo
  enddo
  close(1)
  
  ! Four-index tensor
  open(2,file = 'Four_index_tensor.dat')
  do i=1,ao_num
     do j=1,ao_num
        do k=1,ao_num
           do l=1,ao_num
              accu = get_ao_two_e_integral(i,k,j,l,ao_integrals_map)
              write(2,*)i,j,k,l,accu
           enddo
        enddo
     enddo
  enddo


  close(2)

  !Core hamiltonian
  open(3,file = 'H_core.dat')
  do i=1,ao_num
     do j=1,ao_num
        accu = ao_one_e_integrals(i,j)
        write(3,*)i,j,accu
     enddo
  enddo

  close(3)

  !Nums Orbitals
  open(4,file = 'num_orbitals.dat')
  accu = elec_alpha_num-elec_beta_num
  write(4,*)ao_num,elec_beta_num,accu
  close(4)

  ! Call 'init_guess' subroutine of initialization.irp.f ?

  ! ######################################################################
  ! #Added 03/11/20, Robert, copied from initialization.irp.f (ROHF_GRAD_DESCENT plugin)

  call huckel_guess()
  ! mo_coef = variable globale de QP ; actualisée par la routine huckel_guess
  ! ==> après appel de la routine huckel_guess(), tout appel à mo_coef ou mo_coef_transp 
  ! donne les coefficients donnés par le guess de Huckel

  Pd_tmp = 0.d0
  Ps_tmp = 0.d0

  do j=1,ao_num
     do i=1,ao_num
        accu_d = 0.d0
        accu_s = 0.d0
        do k=1,elec_beta_num
           accu_d += mo_coef_transp(k,i)*mo_coef_transp(k,j)
           !print*, "mo_coef_transp(k,i)*mo_coef_transp(k,j)"
           !print*, mo_coef_transp(k,i)*mo_coef_transp(k,j)
        enddo
        do k=elec_beta_num+1,elec_alpha_num
           accu_s +=  mo_coef_transp(k,i)*mo_coef_transp(k,j)
           !print*, "mo_coef_transp(k,i)*mo_coef_transp(k,j)"
           !print*, mo_coef_transp(k,i)*mo_coef_transp(k,j)
        enddo
        Pd_tmp(i,j) = accu_d
        Ps_tmp(i,j) = accu_s
     enddo
  enddo
  print*, "Test projectors Pd, Ps : norm(Pd*S*Pd-Pd)+norm(Ps*S*Ps-Ps)+(Tr(S*Pd)-Nd)-(Tr(S*Ps)-Ns)"
  print*, test_projs(Pd_tmp,Ps_tmp)

  !###### OPTIONAL : Import guess written in files "Pd_huckel.dat" and "Ps_huckel.dat"
  open(5,file = 'Pd_huckel.dat')
  do i=1,ao_num
     do j=1,ao_num
        write(5,*)i,j,Pd_tmp(i,j)
     enddo
  enddo
  close(5)
  
  open(6,file = 'Ps_huckel.dat')
  do i=1,ao_num
     do j=1,ao_num
        write(6,*)i,j,Ps_tmp(i,j)
     enddo
  enddo
  close(6)

  ! Print the Huckel initial MOs in a file 
  ! for comparison with GAMESS Huckel initial orbitals :
  open(7,file = 'mo_huckel.dat')
  do i=1,ao_num
     do k=1,ao_num
        write(7,*)i,k,mo_coef_transp(k,i)
     enddo
  enddo
  close(7)

  ! ################

  print*,"############ EXTRACTION DONE ############"

  !######################################### Write directly in an Julia array format

  
  ! ! Overlap matrix
  ! open(10, file = "working_dir/overlap_matrix.dat")
  ! do i = 1, ao_num
  !    write(10,'(100(F16.10,X))')overlap_matrix(i,:)
  ! enddo
  ! close(10)

  ! ! Nums_orbitals
  ! open(2, file = "working_dir/nums_orbitals.dat")
  ! write(2, '(3(I4,X))')ao_num,elec_beta_num,elec_alpha_num-elec_beta_num
  ! close(2)


  ! ! H_core
  ! double precision, allocatable :: H_core(:,:)
  ! allocate(H_core(ao_num,ao_num))
  ! do j=1,ao_num
  !    do i = 1,ao_num
  !       H_core(i,j) = ao_one_e_integrals(i,j)
  !    enddo
  ! enddo
  
  ! open(3, file = "working_dir/H_core.dat")
  ! do i = 1, ao_num
  !    write(3,'(100(F16.10,X))')H_core(i,:)
  ! enddo
  ! close(3)


  ! !Four index tensor
  ! open(4, file = "working_dir/Four_index_tensor.dat")
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

  ! print*,"############ THE JOB IS DONE ############"

  mo_label = "Guess"
  call save_mos
  
end program Robert_extract_data

double precision function test_projs(Pd,Ps) result(result)
  implicit none
  BEGIN_DOC
  ! Test the properties of Pd and Ps. The test must be very close to 0.
  END_DOC

  integer                                     :: i,j,k,l
  double precision                            :: accu_d,accu_s, out
  double precision, dimension(ao_num,ao_num)  :: temp_matrix

  double precision, dimension(ao_num,ao_num), intent(in) :: Pd,Ps

  double precision, dimension(ao_num,ao_num)  :: S_Pd,S_Ps

  out = 0.d0

  ! Norm of Pd*S*Pd - Pd (S=overlap matrix)
  temp_matrix = matmul(matmul(Pd,ao_overlap),Pd) - Pd
  out += Norm2(temp_matrix)
  
  ! Norm of Ps*S*Ps - Ps
  temp_matrix= matmul(matmul(Ps,ao_overlap),Ps) - Ps
  out += Norm2(temp_matrix)

  !Norm of Ps*Pd
  temp_matrix = matmul(Pd,Ps)
  out  += Norm2(temp_matrix)

  !Test traces
  S_Pd=matmul(ao_overlap,Pd)
  accu_d = 0.d0
  accu_s = 0.d0
  do i=1,ao_num
     accu_d += S_Pd(i,i)
  enddo
  
  S_Ps=matmul(ao_overlap,Ps)
  do i=1,ao_num
     accu_s += S_Ps(i,i)
  enddo
  out += accu_d - elec_beta_num
  out -= accu_s - (elec_alpha_num-elec_beta_num)

  result = out
  

end function test_projs
