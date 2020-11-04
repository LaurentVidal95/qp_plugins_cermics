program ROHF_GRAD_DESCENT
  implicit none
  BEGIN_DOC
  ! Fixed step (or fixed by part) gradient descent algorithm for minimizing the
  ! ROHF energy on the manifold of projectors (Pd,Ps) where :
  !
  ! Pd is the projector on doubly occupied orbitals
  ! Ps is the projector on singly occupied orbitals
  !
  ! TODO : add asserts !
  !
  END_DOC

  integer                                       :: iter, i,j,k,l
  double precision                              :: get_rohf_energy, test_projs !Functions
  double precision                              :: energy_prev, energy, delta_energy, test, grad_norm
  double precision, dimension(ao_num,ao_num)    :: Pd,Ps,P,tmp_mat   
  double precision, dimension(ao_num,ao_num)    :: Gd,Gs, eigvectors, U
  double precision, dimension(ao_num)           :: eigvalues
  
  double precision, dimension(ao_num,ao_num)    :: hd,hs,Qd,Qs,test_Qd,test_Qs,Id,Pv !TEST
  double precision                              :: accu_1,accu_2,accu_test ! TEST
   !TEST



  ! Initial data
  call init_guess(Pd,Ps)

  PROVIDE max_iter threshold step !Found in EZFIO.cfg
  iter = 0
  ! threshold = 1e-15
  delta_energy = 1.d0
  step = 0.008d0 !Adapt the initial step to the system. Greater than 0.05d0 affects convergence.
  
  energy = get_rohf_energy(Pd,Ps)
  test = test_projs(Pd,Ps)
  
  call write_time(6)
  
  call write_double(6,energy, 'Guess energy = ')
  call write_double(6, step,'Step  = ')

  call write_double(6, abs(test), 'Initial test = ')
  
  write(6,'(A6, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    '======','================','================','================','================'
  write(6,'(A6, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
    ' N ', 'Energy  ', 'Energy diff  ',' Test','  Gradient Norm'
  write(6,'(A6, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
       '======','================','================','================','================'

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !           Main Loop             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  do while ( (iter < max_iter) .and. (abs(delta_energy) > threshold) )

     iter +=1
     energy_prev = energy

     call get_grad_projs(Pd,Ps,Gd,Gs)

     grad_norm = NORM2(Gd) + NORM2(Gs)

     !P
     tmp_mat = Pd + Ps - step*Gd - step*Gs
     call lapack_diagd(eigvalues,U,tmp_mat,ao_num,ao_num)
     call retraction(tmp_mat,U,P)

     !Pd
     tmp_mat = matmul(P,matmul(Pd - step*Gd,P))
     call lapack_diagd(eigvalues,U,tmp_mat,ao_num,ao_num)
     call retraction(tmp_mat,U,Pd) 

     !Ps
     Ps = P - Pd

     test = test_projs(Pd,Ps)
     
     ! Compute new delta__energy
     energy  = get_rohf_energy(Pd,Ps) 
     delta_energy = energy -  energy_prev

     !Adapt step if the energy is increasing
     if(energy > energy_prev) then
        step = 0.5d0*step
        print*,'HALF STEP'
     endif

     ! call write_double(6,test,'Test')
     write(6,'(I6, 1X, F16.10, 1X, F16.12, 1X, F16.12, 1X, F16.12)')&
          iter, energy, delta_energy, test, grad_norm
     
  enddo


    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !         END MAIN LOOP           !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  write(6,'(A6, 1X, A16, 1X, A16, 1X, A16, 1X, A16)')  &
       '======','================','================','================','================'
  write(6,*)

  call write_double(6, energy, 'ROHF energy')

  call write_time(6)

  if (iter < max_iter) then
     print*,'CONVERGED at iteration',iter
  else
     print*,'NOT CONVERGED'
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !        WRITE FINAL DMs          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(10, file = "data_dir/Pd.dat")
  do i=1,ao_num
     write(10,'(100(F16.10,X))')Pd(i,:)
  enddo
  close(10)

  open(11, file = "data_dir/Ps.dat")
  do i=1,ao_num
     write(11,'(100(F16.10,X))')Ps(i,:)
  enddo
  close(11)
  

  mo_label = "Guess"
  call save_mos
  
   
end program ROHF_GRAD_DESCENT

