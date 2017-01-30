Solver packages
^^^^^^^^^^^^^^^
One of the following solvers must be specified for a model run. All arguments for the solver are non-spatial values.

PCG package
~~~~~~~~~~~
The preconditioned conjugate-gradient package can be enabled with

.. code-block:: c

   res = mf::setPCG(MXITER, ITERI, NPCOND, HCLOSE, RCLOSE, RELAX, NBPOL, DAMP);

where

MXITER
   is the maximum number of outer iterations;

ITERI
   is the number of inner iterations;

NPCOND
   1 - Modified Incomplete Cholesky, 2 - Polynomial matrix conditioning method;

HCLOSE
   is the head change criterion for convergence;

RCLOSE
   is the residual criterion for convergence;

RELAX
   is the relaxation parameter used with NPCOND = 1;

NBPOL
   indicates whether the estimate of the upper bound on the maximum eigenvalue is 2.0 and

DAMP
   is the damping factor.

SOR package
~~~~~~~~~~~
The slice-successive overrelaxation package can be enabled with

.. code-block:: c

   res = mf::setSOR(MXITER, ACCL, HCLOSE);

where

MXITER
   is the maximum number of iterations allowed in a time step;

ACCL
   is the acceleration variable and

HCLOSE
   the head change criterion for convergence.

SIP package
~~~~~~~~~~~
The strongly implicit procedure package can be enabled with

.. code-block:: c

   res = mf::setSIP(MXITER, NPARAM, ACCL, HCLOSE, IPCALC, WSEED);

where

MXITER
   is the maximum number of times through the iteration loop in one time step;

NPARAM
   is the number of iteration variables to be used;

ACCL
   is the acceleration variable;

IPCALC
   0 - the seed entered by the user will be used, 1 - the seed will be calculated at the start of the simulation, and

WSEED
   is the seed for calculating iteration variables.

DSP package
~~~~~~~~~~~
The direct solver package can be enabled with

.. code-block:: c

   res = mf::setDSP(ITMX, MXUP, MXLOW, MXBW, IFREQ, ACCL, HCLOSE);

where

ITMX
   is the maximum number of iterations each time step;

MXUP
   is the maximum number of equations in the upper part of the equations to be solved;

MXLOW
   is the maximum number of equations in the lower part of equations to be solved;

MXBW
   is the maximum band width plus 1 of the [AL] matrix;

IFREQ
   is flag indicating the frequency at which coefficients in [A] change;
ACCL
   is a multiplier for the computed head change for each iteration and
HCLOSE
   is the head change closure criterion.
