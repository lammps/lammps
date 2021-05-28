.. index:: fix polarize/bem/icc
.. index:: fix polarize/bem/gmres
.. index:: fix polarize/functional

fix polarize command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID polarize nevery tolerance ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* polarize/bem/icc, polarize/bem/gmres, or polarize/functional  = style name of this fix command
* Nevery = this fixed is invoked every this many timesteps
* tolerance = the tolerance for the iterative solver to stop


Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all polarize/bem/icc 1 0.0001
   fix 2 all polarize/bem/gmres 5 0.0001
   fix 3 all polarize/bem/functional 1 0.0001

Description
"""""""""""

The three fix polarize in the USER-DIELECTRIC package compute the induced charges
at the interface between two impermeable media with different dielectric
constants.

There are some example scripts for using this package with LAMMPS in the
examples/USER/dielectric directory.

----------

The charges of the atoms in the specified group will be computed by the solver.
fix polarize bem/icc computes the induced charges at the boundary elements
(i.e. interface vertices) using the successive overrelaxation as described
in (Tyagi). fix polarize bem/gmres computes the induced charges at
the interface vertices using the successive overrelaxation
as described in (Barros).

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`fix_modify <fix_modify>` 

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the Colvars
energy mentioned above.  The scalar value calculated by this fix is
"extensive".

Restrictions
""""""""""""

This fix is part of the USER-DIELECTRIC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` doc page for more info.

There can only be one colvars fix active at a time. Since the interface
communicates only the minimum amount of information and colvars module
itself can handle an arbitrary number of collective variables, this is
not a limitation of functionality.

Related commands
""""""""""""""""

:doc:`fix smd <fix_smd>`, :doc:`fix spring <fix_spring>`,

Default
"""""""

None.

----------

.. _NguyenTD:

**(NguyenTD)** Nguyen, Li, Bagchi, Solis, Olvera de la Cruz, Mol. Phys., DOI:10.1016/j.cpc.2019.03.006

