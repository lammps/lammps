.. index:: bond_style mesocnt

bond_style mesocnt command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style mesocnt

Examples
""""""""

.. code-block:: LAMMPS

   bond_style mesocnt
   bond_coeff 1 C 10 10 20.0
   bond_coeff 4 custom 800.0 10.0

Description
"""""""""""

.. versionadded:: 15Sep2022

The *mesocnt* bond style is a wrapper for the :doc:`harmonic
<bond_harmonic>` style, and uses the potential

.. math::

   E = K (r - r_0)^2

where :math:`r_0` is the equilibrium bond distance.  Note that the
usual 1/2 factor is included in :math:`K`.  The style implements
parameterization presets of :math:`K` for mesoscopic simulations of
carbon nanotubes based on the atomistic simulations of
:ref:`(Srivastava) <Srivastava_1>`.

Other presets can be readily implemented in the future.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data
<read_data>` or :doc:`read_restart <read_restart>` commands:

* preset = *C* or *custom*
* additional parameters depending on preset

Preset *C* is for carbon nanotubes, and the additional parameters are:

* chiral index :math:`n` (unitless)
* chiral index :math:`m` (unitless)
* :math:`r_0` (distance)

Preset *custom* is simply a direct wrapper for the :doc:`harmonic
<bond_harmonic>` style, and the additional parameters are:

* :math:`K` (energy/distance\^2)
* :math:`r_0` (distance)

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the MOLECULE
and MESONT packages.  See the :doc:`Build package <Build_package>`
page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none

----------

.. _Srivastava_1:

**(Srivastava)** Zhigilei, Wei and Srivastava, Phys. Rev. B 71, 165417
(2005).
