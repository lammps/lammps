.. index:: fix polarize/functional

fix polarize/functional command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID polarize nevery tolerance ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* polarize/functional  = style name of this fix command
* Nevery = this fixed is invoked every this many timesteps
* tolerance = the tolerance for the iterative solver to stop


Examples
""""""""

.. code-block:: LAMMPS

   fix 3 all polarize/functional 1 0.001

Description
"""""""""""

The three fix polarize in the USER-DIELECTRIC package compute the induced charges
at the interface between two impermeable media with different dielectric
constants.

There are some example scripts for using this package with LAMMPS in the
examples/USER/dielectric directory.

----------


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

:doc:`fix polarize/bem/icc <fix_polarize_bem_icc>`, :doc:`fix polarize/functional <fix_polarize_bem_gmres>`

Default
"""""""

None.

----------

.. _Jadhao:

**(Jadhao)** Jadhao, Solis, Olvera de la Cruz, J Chem Phys, 138, 054119 (2013)

.. _NguyenTD:

**(NguyenTD)** Nguyen, Li, Bagchi, Solis, Olvera de la Cruz, Comput Phys Commun 241, 80-19 (2019)

