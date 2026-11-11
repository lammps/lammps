.. index:: fix spring/tspin

fix spring/tspin command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID spring/tspin K S0

* ID, group-ID are documented in :doc:`fix <fix>` command
* spring/tspin = style name of this fix command
* K = spring constant (energy per squared Bohr magneton,
  :math:`\mathrm{eV}/\mu_B^2` in metal units)
* S0 = equilibrium spin modulus (Bohr magnetons)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all spring/tspin 1.0 2.2
   fix_modify 1 energy yes

Description
"""""""""""

.. versionadded:: TBD

Apply a harmonic restoring potential to the magnitude of each
magnetic spin in the group,

.. math::

   U_{L,i} = \frac{1}{2}K(|\vec{S}_i|-S_0)^2,

which produces the radial magnetic force

.. math::

   \vec{F}^{m}_{i,\mathrm{spring}} =
   -K(|\vec{S}_i|-S_0)\hat{s}_i.

The spring force is purely radial and produces no direct spin torque.
This fix is useful with :doc:`TSPIN integrators <fix_nve_tspin>` when
the magnetic potential does not provide a longitudinal restoring force.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

Enable :doc:`fix_modify energy yes <fix_modify>` to include the spring
energy in the global potential energy *pe*; it is excluded by default.
This fix also computes a global scalar containing the total spring
energy, accessible as *f_ID*.  This scalar is "extensive" and has
energy units.

The spring force is applied during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix requires :doc:`atom_style spin <atom_style>` and the SPIN
package; see :doc:`Build package <Build_package>`.

Related commands
""""""""""""""""

:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`

Default
""""""""

none
