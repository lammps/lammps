**SELM Library:** Fluctuating hydrodynamics simulation methods.
Provides integrators based on stochastic eulerian lagrangian methods 
(SELMs) for performing simulations of hydrodynamic coupling subject 
to thermal fluctuations by using stochastic continuum fields.  
Methods include stochastic immersed boundary methods, stochastic 
eulerian lagrangian methods, and implicit-solvent coarse-grained 
methods.  Details can be found in the papers

* Fluctuating Hydrodynamics Methods for Dynamic Coarse-Grained 
Implicit-Solvent Simulations in LAMMPS, Y. Wang, J. K. Sigurdsson, 
and P.J. Atzberger, SIAM J. Sci. Comput. , 38(5), S62–S77, (2016), 
https://doi.org/10.1137/15M1026390

* Stochastic Reductions for Inertial Fluid-Structure Interactions 
Subject to Thermal Fluctuations, G. Tabak and P.J. Atzberger, 
SIAM J. Appl. Math., 75(4), 1884–1914, (2015),
http://dx.doi.org/10.1016/j.jtbi.2010.11.023

* Incorporating Shear into Stochastic Eulerian Lagrangian Methods for 
Rheological Studies of Complex Fluids and Soft Materials., 
P.J. Atzberger, Physica D, Vol. 265, pg. 57–70, (2013),
http://dx.doi.org/10.1016/j.physd.2013.09.002

For example scripts, python jupyter notebooks, pre-compiled binaries, 
and additional information see
http://mango-selm.org

These codes provide default version for library ``libselm``.  The 
USER-SELM provides ``fix selm`` for interfacing this library. The 
latest version of codes and binaries can be downloaded from 
the selm website http://mango-selm.org

NOTE: The current version of the codes should be run in serial mode.
A head node is used to handle the fluid mechanics.

This package was developed and is maintained by

Paul J. Atzberger
University of California Santa Barbara
atzberg@gmail.com

--------------------------------------------------------------------------
Example compilation of the codes is given by Makefile.serial.

Build ``libselm`` using for example

  ``make -f Makefile.serial lib_shared``

  or 

  ``make -f Makefile.serial lib_static``

  Note, the shared library ``lib_shared`` build is preferred.

--------------------------------------------------------------------------

For more information, example scripts, and python jupyter notebooks, see 
http://mango-selm.org

