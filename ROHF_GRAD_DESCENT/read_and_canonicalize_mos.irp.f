program read_and_canonicalize_mos
  implicit none
  
  double precision                :: energy
  double precision, allocatable   :: new_mos(:,:)
  
  allocate(new_mos(ao_num,ao_num))
  open(39,file = 'data_dir/new_orbitals.dat')
  read(39,*) new_mos
  close(39) 
  
  mo_coef = transpose(new_mos)
  touch mo_coef
  mo_label = "Natural"
  call save_mos

  energy = scf_energy   ! Time consuming ..

  
  ! print*,"Energy for the new MOs : ",energy

  call roothaan_hall_scf()
  
end program read_and_canonicalize_mos
