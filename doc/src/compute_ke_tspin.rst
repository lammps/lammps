.. index:: compute ke/tspin

compute ke/tspin command
========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID ke/tspin

* ID, group-ID are documented in :doc:`compute <compute>` command
* ke/tspin = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute ske all ke/tspin

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the kinetic energy stored in the
spin degrees of freedom of inertial spin dynamics for all atoms in the
group,

.. math::

   E^{s}_{kin} = \sum_i \frac{1}{2} m_s \left| \vec{v}^{s}_i \right|^2

where :math:`m_s` is the spin mass of atom i and :math:`\vec{v}^{s}_i`
its spin velocity.

This energy is not included in the thermodynamic keyword *ke*, which
counts only the translational kinetic energy of the atoms.  Reference
this compute with a *c_ID* keyword in :doc:`thermo_style custom
<thermo_style>` to output it, and add it to the potential energy to
check the energy conservation of a :doc:`fix nve/tspin
<fix_nve_tspin>` run:

.. code-block:: LAMMPS

   compute      ske all ke/tspin
   variable     etot equal pe+c_ske
   thermo_style custom step pe c_ske v_etot

Output info
"""""""""""

This compute calculates a global scalar (the spin kinetic energy).
This value can be used by any command that uses a global scalar value
from a compute as input.  See the :doc:`Howto output <Howto_output>`
page for an overview of LAMMPS output options.

The scalar value calculated by this compute is "extensive".  The scalar
value will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

The *ke/tspin* compute is part of the SPIN package.  This style is only
enabled if LAMMPS was built with this package.  See the :doc:`Build
package <Build_package>` page for more info.

This compute requires :doc:`atom_style spin <atom_style>`.

Related commands
""""""""""""""""

:doc:`compute ke <compute_ke>`,
:doc:`compute spin <compute_spin>`,
:doc:`fix nve/tspin <fix_nve_tspin>`

Default
""""""""

none
