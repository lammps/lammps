.. index:: dihedral_style helix
.. index:: dihedral_style helix/omp

dihedral_style helix command
============================

Accelerator Variants: *helix/omp*

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style helix

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style helix
   dihedral_coeff 1 80.0 100.0 40.0

Description
"""""""""""

The *helix* dihedral style uses the potential

.. math::

   E = A [1 - \cos(\theta)] + B [1 + \cos(3 \theta)] +
       C [1 + \cos(\theta + \frac{\pi}{4})]

This coarse-grain dihedral potential is described in :ref:`(Guo) <Guo>`.
For dihedral angles in the helical region, the energy function is
represented by a standard potential consisting of three minima, one
corresponding to the trans (t) state and the other to gauche states
(g+ and g-).  The paper describes how the :math:`A`, :math:`B` and,
:math:`C` parameters are chosen so as to balance secondary (largely
driven by local interactions) and
tertiary structure (driven by long-range interactions).

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`A` (energy)
* :math:`B` (energy)
* :math:`C` (energy)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This dihedral style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`

Default
"""""""

none

----------

.. _Guo:

**(Guo)** Guo and Thirumalai, Journal of Molecular Biology, 263, 323-43 (1996).
