=================
ROHF_GRAD_DESCENT
=================

* WORK IN PROGRESS

This modules contains the basics to perform a fixed step (or fixed by part) gradient descent algorithm for minimizing the ROHF energy.

The details of the algorithms are given in the report.

-----------------------------

There are three files :

1) :file:`initialization.irp.f` contains the providers for :
   * the overlap matrice, its square root and inverse square root.
   * the routine |init_guess| providing initial densities. By defaut : Core Guess
   * the function |test_projs| that tests whether the projectors are in the right space.
   ie if : Tr(Pd) = Nd, Tr(Ps) = Ns, Pd^2 = Pd, Ps^2 = Ps, Pd*Ps = 0

2) :file:`energy_and_gradients.irp.f` contains the providers and subroutine to compute
   the energy and the projected gradient at each loop.
   It also contains the subroutine |retraction|.

3) :file:`ROHF_GRAD_DESCENT.irp.f` is the main file.
