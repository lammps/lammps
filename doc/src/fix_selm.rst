.. index:: fix selm

fix selm command
==================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID selm file_input ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* selm = style name of this fix command

  .. parsed-literal::

       *file_input* arg = name of input parameter file to use for selm codes (default: NULL)

Examples
""""""""

.. code-block:: LAMMPS

   fix s1 all selm parameters.xml

Description
"""""""""""
This fix provides methods for fluctuating hydrodynamics simulations using
continuum stochastic fluid equations.  This includes integrators based on
stochastic eulerian lagrangian methods (SELMs) for performing simulations
of hydrodynamic coupling subject to thermal fluctuations by using
continuum fields.  Methods include stochastic immersed boundary methods,
stochastic eulerian lagrangian methods, and implicit-solvent coarse-grained
methods.  More details can be found in the papers :ref:`(SELM LAMMPS) <SELM_LAMMPS>`,
:ref:`(SELM Reduc) <SELM_Reduc>`, :ref:`(SELM Shear) <SELM_Shear>`

The package is organized to call routines in the :ref:`SELM Library <user-selm>`.
A version of the library codes can be found in the directory ``lib/selm``.
The latest library version and development is hosted at https://mango-selm.org

Examples LAMMPS scripts can also be found in directory
``examples/USER/selm``

For addition information including Jupyter python notebooks, example
scripts, and pre-compiled binaries, see http://mango-selm.org

----------

The package is controlled by an XML format parameter file as discussed in
the examples and on-line documentation.

Currently, the *group-ID* entry is ignored. LAMMPS will always pass
all the atoms to SELM and there should only be one instance of the
selm fix at a time.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When performing a restart of a calculation that involves SELM you must
specify information in the parameter file. No part of the SELM restart
data is included in the LAMMPS restart files.

The :doc:`fix_modify <fix_modify>` option is not supported by this fix.

This fix does not provide quantities that can be accessed by the
:doc:`output commands <Howto_output>` within LAMMPS.  The data is
instead output to formats such as VTK and XML by methods native to
SELM, such as the hydrodynamic fields and other information.

Restrictions
""""""""""""

Fix selm can only be used in serial mode. It can be compiled
with MPI enabled, but then LAMMPS must only be run with a single
MPI process.  LAMMPS will stop with an error, if this condition
is not met.

There may only be one instance of fix selm at any time.  LAMMPS
will stop with an error, if this condition is not met.

The fix is part of the USER-SELM package.  It is only enabled if
LAMMPS was built with this package.  For more information see the
:doc:`Build package <Build_package>` document page.

Related commands
""""""""""""""""

Default
"""""""

The default options are params_file = NULL

----------

.. _SELM_LAMMPS:

**(SELM LAMMPS)** *Fluctuating Hydrodynamics Methods for Dynamic Coarse-Grained Implicit-Solvent Simulations in LAMMPS,* Y. Wang, J. K. Sigurdsson, and P.J. Atzberger, SIAM J. Sci. Comput. , 38(5), S62-S77, (2016), https://doi.org/10.1137/15M1026390

.. _SELM_Reduc:

**(SELM Reduc)** *Stochastic Reductions for Inertial Fluid-Structure Interactions Subject to Thermal Fluctuations,* G. Tabak and P.J. Atzberger, SIAM J. Appl. Math., 75(4), 1884-1914, (2015), http://dx.doi.org/10.1137/15M1019088

.. _SELM_Shear:

**(SELM Shear)** *Incorporating Shear into Stochastic Eulerian Lagrangian Methods for Rheological Studies of Complex Fluids and Soft Materials,* P.J. Atzberger, Physica D, Vol. 265, pg. 57-70, (2013), http://dx.doi.org/10.1016/j.physd.2013.09.002

