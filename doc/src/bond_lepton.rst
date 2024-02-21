.. index:: bond_style lepton
.. index:: bond_style lepton/omp

bond_style lepton command
=========================

Accelerator Variants: *lepton/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style style args

* style = *lepton*
* args = optional arguments

.. parsed-literal::

   args = *auto_offset* or *no_offset*
     *auto_offset* = offset the potential energy so that the value at r0 is 0.0 (default)
     *no_offset* = do not offset the potential energy

Examples
""""""""

.. code-block:: LAMMPS

   bond_style lepton
   bond_style lepton no_offset

   bond_coeff  1  1.5 "k*r^2; k=250.0"
   bond_coeff  2  1.1 "k2*r^2 + k3*r^3 + k4*r^4; k2=300.0; k3=-100.0; k4=50.0"
   bond_coeff  3  1.3 "k*r^2; k=350.0"

Description
"""""""""""

.. versionadded:: 8Feb2023

Bond style *lepton* computes bonded interactions between two atoms with
a custom function.  The potential function must be provided as an
expression string using "r" as the distance variable relative to the
reference distance :math:`r_0` which is provided as a bond coefficient.
For example `"200.0*r^2"` represents a harmonic potential with a force
constant *K* of 200.0 energy units:

.. math::

   U_{bond,i} = K (r_i - r_0)^2 = K r^2 \qquad r = r_i - r_0

.. versionchanged:: 7Feb2024

By default the potential energy U is shifted so that he value U is 0.0
for $r = r_0$.  This is equivalent to using the optional keyword
*auto_offset*.  When using the keyword *no_offset* instead, the
potential energy is not shifted.

The `Lepton library <https://simtk.org/projects/lepton>`_, that the
*lepton* bond style interfaces with, evaluates this expression string at
run time to compute the pairwise energy.  It also creates an analytical
representation of the first derivative of this expression with respect to
"r" and then uses that to compute the force between the atom pairs forming
bonds as defined by the topology data.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the examples above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* Lepton expression (energy units)
* :math:`r_0` (distance)

The Lepton expression must be either enclosed in quotes or must not
contain any whitespace so that LAMMPS recognizes it as a single keyword.
More on valid Lepton expressions below.  The :math:`r_0` is the
"equilibrium distance".  The potential energy function in the Lepton
expression is shifted in such a way, that the potential energy is 0 for
a bond length :math:`r_i == r_0`.

----------

.. include:: lepton_expression.rst

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style is part of the LEPTON package and only enabled if LAMMPS
was built with this package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`bond_style table <bond_table>`,
:doc:`bond_write <bond_write>`, :doc:`angle_style lepton <angle_lepton>`,
:doc:`dihedral_style lepton <dihedral_lepton>`

Default
"""""""

none
