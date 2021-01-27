 BEGIN_PROVIDER [ double precision, coulomb_d, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, coulomb_s, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, exchange_d, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, exchange_s, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
! J operator for P_d :: J(i,j) = \sum_{k,l} (ij|kl) P_d(k,l)
 END_DOC 
 use map_module
 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer*8                      :: p,q
 double precision               :: integral, c0, c1, c2
 double precision               :: ao_two_e_integral, local_threshold
 double precision, allocatable  :: coulomb_s_tmp(:,:)
 double precision, allocatable  :: coulomb_d_tmp(:,:)
 double precision, allocatable  :: exchange_s_tmp(:,:)
 double precision, allocatable  :: exchange_d_tmp(:,:)

 PROVIDE ao_two_e_integrals_in_map

 integer(omp_lock_kind) :: lck(ao_num)
 integer(map_size_kind)     :: i8
 integer                        :: ii(8), jj(8), kk(8), ll(8), k2
 integer(cache_map_size_kind)   :: n_elements_max, n_elements
 integer(key_kind), allocatable :: keys(:)
 double precision, allocatable  :: values(:)
 coulomb_d = 0.d0
 coulomb_s = 0.d0
 exchange_d = 0.d0
 exchange_s = 0.d0

!!$OMP PARALLEL DEFAULT(NONE)                                      &
!    !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
!    !$OMP  n_elements,coulomb_s_tmp,coulomb_d_tmp,exchange_s_tmp,exchange_d_tmp)&
!    !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
!    !$OMP  ao_integrals_map, coulomb_s, coulomb_d, density_mat_d, density_mat_s)

  call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
  allocate(keys(n_elements_max), values(n_elements_max))
  allocate(coulomb_s_tmp(ao_num,ao_num), &
           coulomb_d_tmp(ao_num,ao_num))
  allocate(exchange_s_tmp(ao_num,ao_num), &
           exchange_d_tmp(ao_num,ao_num))
  coulomb_s_tmp = 0.d0
  coulomb_d_tmp  = 0.d0
  exchange_s_tmp = 0.d0
  exchange_d_tmp  = 0.d0

  !$OMP DO SCHEDULE(static,1)
  do i8=0_8,ao_integrals_map%map_size
    n_elements = n_elements_max
    call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
    ! keys(n_elements) = i,j,k,l ; values(n_elements) = integral
    do k1=1,n_elements
      ! get the i,j,k,l with all permutations for the key == keys(k1)
      call two_e_integrals_index_reverse(kk,ii,ll,jj,keys(k1))

      do k2=1,8 ! loop over all 8 permutations for i,j,k,l
        if (kk(k2)==0) then
          cycle
        endif
        i = ii(k2)
        j = jj(k2)
        k = kk(k2)
        l = ll(k2)
!       integral = (SCF_density_matrix_ao_alpha(k,l)+SCF_density_matrix_ao_beta(k,l)) * values(k1)
        integral = values(k1) 
        coulomb_d_tmp(i,j)   += density_mat_d(k,l) * integral 
        coulomb_s_tmp(i,j)   += density_mat_s(k,l) * integral 
        exchange_s_tmp(l,j)  += density_mat_s(k,i) * integral
        exchange_d_tmp(l,j)  += density_mat_d(k,i) * integral
      enddo
    enddo
  enddo
! !$OMP END DO NOWAIT
! !$OMP CRITICAL
  coulomb_s += coulomb_s_tmp
  coulomb_d  += coulomb_d_tmp
  exchange_s += exchange_s_tmp
  exchange_d  += exchange_d_tmp
! !$OMP END CRITICAL
  deallocate(keys,values,coulomb_s_tmp,coulomb_d_tmp)
  deallocate(exchange_s_tmp,exchange_d_tmp)
! !$OMP END PARALLEL

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_one_e_integrals_ortho, (ao_num, ao_num)]
 implicit none
 call ao_ortho_to_ao(ao_one_e_integrals, ao_one_e_integrals_ortho)
END_PROVIDER 
