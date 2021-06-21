**USER-SELM Package:** Fluctuating hydrodynamics simulation methods.
Provides integrators based on stochastic eulerian lagrangian methods
(SELMs) for performing simulations of hydrodynamic coupling subject
to thermal fluctuations by using stochastic continuum fields.
Methods include stochastic immersed boundary methods, stochastic
eulerian lagrangian methods, and implicit-solvent coarse-grained
methods.

This package was developed and is maintained by

Paul J. Atzberger
University of California Santa Barbara
http://atzberger.org
atzberg@gmail.com

Details can be found in the papers

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

NOTE: This package currently should be run in serial mode.
The head node handles the continuum fluid mechanics.

NOTE: The USER-SELM interfaces codes organized into a library and
requires building the SELM library ``libselm`` for operation.
A default version of the library codes can be found in ``/lib/selm``
along with build scripts and instructions.  The latest
version of codes and binaries can be downloaded from
the selm website http://mango-selm.org/

NOTE: This package depends on MOLECULE which must be installed,
for example ``make yes-molecule``.

NOTE: While not strictly required, many of the example scripts make
use of USER-VTK, which ideally should be installed.  For example,
``make yes-user-vtk``.

--------------------------------------------------------------------------
Fixes provided by this package:

fix_selm.cpp: fluctuating hydrodynamics simulations using continuum
              stochastic fields.  Parameters specified in a related
              input XML file (see examples and docs for details).

Example usage:

``fix s1 all selm parameters.xml``

--------------------------------------------------------------------------

To build the package copy files from USER-SELM/MAKE to src/MAKE/MINE
and editing as needed.  Then use

``make yes-molecule``
``make yes-user-selm``
``make atz_selm_serial mode=shlib LMP_INC="-DLAMMPS_EXCEPTIONS" JPG_LIB="-lpng -ljpeg"``

(optional) For VTK support in examples suggested also to use before build
``make yes-user-vtk``

NOTE: For the shared library, you make need to adjust library search path
that includes the ``libselm``.  For example, on linux you may need to set
LD_LIBRARY_PATH to include location of ``libselm.so``, such as

``LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../lib/selm``

or

``LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(LAMMPS_DIR)/lib/selm``

--------------------------------------------------------------------------

For more information, example scripts, python jupyter notebooks
see http://mango-selm.org


