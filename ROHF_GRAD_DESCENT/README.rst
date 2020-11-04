=================
ROHF_GRAD_DESCENT
=================

* WORK IN PROGRESS

This modules contains the basics to perform a fixed step (or fixed by part) gradient descent algorithm for minimizing the ROHF energy.

The details of the algorithms are given in the report.

-----------------------------

There are three main files :

1) :file:`initialization.irp.f` contains the providers for :
   * the overlap matrice, its square root and inverse square root.
   * the routine |init_guess| providing initial densities. By defaut : Core Guess
   * the function |test_projs| that tests whether the projectors are in the right space.
   ie if : Tr(Pd) = Nd, Tr(Ps) = Ns, Pd^2 = Pd, Ps^2 = Ps, Pd*Ps = 0

2) :file:`energy_and_gradients.irp.f` contains the providers and subroutine to compute
   the energy and the projected gradient at each loop.
   It also contains the subroutine |retraction|.

3) :file:`ROHF_GRAD_DESCENT.irp.f` is the main file.

-----------------------------
   
There are two other files, needed to obtain canonical MOs from final densities.

The procedure is such :

1) call `$qp run extract_data` to obtain the |H_core| matrix, |overlap matrix|, and the number of total, d and s orbitals.
   The last part of the file generates the |four-index tensor| but is commented since it is not needed for this particular task.

   These files are writen in a `/data_dir` directory.

2) Copy your file :file:`Pd.dat` and :file:`Ps.dat` in `/data_dir`, containing the converged density matrices.
                          [ a b c;             a b c 
   The format is for M =    d e f;     M.dat = d e f
			    g h i]             g h i

   !!! By default, the routine `ROHF_GRAD_DESCENT` generate those files when it converges. !!!

3) run the julia code :file:`extract_mos.jl` (outside `/data dir`). It should generate the file :file:`new_orbitals.dat` in `/data_dir`

4) call `$ qp run read_and_canonicalize_mos`. Orbitals should be canonical.
  
			           
			           
                           
