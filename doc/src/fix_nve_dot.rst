.. index:: fix nve/dot

fix nve/dot command
===================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/dot

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/dot = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/dot

Description
"""""""""""

Apply a rigid-body integrator as described in :ref:`(Davidchack) <Davidchack4>`
to a group of atoms, but without Langevin dynamics.
This command performs Molecular dynamics (MD)
via a velocity-Verlet algorithm and an evolution operator that rotates
the quaternion degrees of freedom, similar to the scheme outlined in :ref:`(Miller) <Miller4>`.

This command is the equivalent of the :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`
without damping and noise and can be used to determine the stability range
in a NVE ensemble prior to using the Langevin-type DOTC-integrator
(see also :doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`).
The command is equivalent to the :doc:`fix nve <fix_nve>`.
The particles are always considered to have a finite size.

An example input file can be found in /examples/PACKAGES/cgdna/examples/duplex1/.
Further details of the implementation and stability of the integrator are contained in :ref:`(Henrich) <Henrich4>`.
The preprint version of the article can be found `here <PDF/CG-DNA.pdf>`_.

----------

Restrictions
""""""""""""

These pair styles can only be used if LAMMPS was built with the
:ref:`CG-DNA <PKG-CG-DNA>` package and the MOLECULE and ASPHERE package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix nve/dotc/langevin <fix_nve_dotc_langevin>`, :doc:`fix nve <fix_nve>`

Default
"""""""

none

----------

.. _Davidchack4:

**(Davidchack)** R.L Davidchack, T.E. Ouldridge, and M.V. Tretyakov. J. Chem. Phys. 142, 144114 (2015).

.. _Miller4:

**(Miller)** T. F. Miller III, M. Eleftheriou, P. Pattnaik, A. Ndirango, G. J. Martyna, J. Chem. Phys., 116, 8649-8659 (2002).

.. _Henrich4:

**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E. Ouldridge, Eur. Phys. J. E 41, 57 (2018).
