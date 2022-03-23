.. index:: fix viscous/sphere

fix viscous/sphere command
==========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID viscous/sphere gamma keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* viscous/sphere = style name of this fix command
* gamma = damping coefficient (torque/angular velocity units)
* zero or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *scale*
       *scale* values = type ratio
         type = atom type (1-N)
         ratio = factor to scale the damping coefficient by

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 flow viscous/sphere 0.1
   fix 1 damp viscous/sphere 0.5 scale 3 2.5

Description
"""""""""""

Add a viscous damping torque to finite-size spherical particles in the group
that is proportional to the angular velocity of the atom.  In granular
simulations this can be useful for draining the rotational kinetic energy from
the system in a controlled fashion.  If used without additional thermostatting
(to add kinetic energy to the system), it has the effect of slowly (or rapidly)
freezing the system; hence it can also be used as a simple energy minimization
technique.

The damping torque :math:`T_i` is given by :math:`T_i = - \gamma \omega_i`.
The larger the coefficient, the faster the rotational kinetic energy is reduced.
If the optional keyword *scale* is used, :math:`\gamma` can be scaled up or
down by the specified factor for atoms of that type.  It can be used
multiple times to adjust :math:`\gamma` for several atom types.

.. note::

   You should specify gamma in torque/angular velocity units.  This is not
   the same as mass/time units, at least for some of the LAMMPS
   :doc:`units <units>` options like "real" or "metal" that are not
   self-consistent.

In the current implementation, rather than have the user specify a viscosity,
:math:`\gamma` is specified directly in torque/angular velocity units.
If needed, :math:`\gamma` can be adjusted for atoms of different sizes
(i.e. :math:`\sigma`) by using the *scale* keyword.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is modifying torques. Default is the outermost level.

The torques due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  This fix should only
be used with damped dynamics minimizers that allow for
non-conservative forces.  See the :doc:`min_style <min_style>` command
for details.

Restrictions
""""""""""""

This fix requires that atoms store torque and angular velocity (omega)
and a radius as defined by the :doc:`atom_style sphere <atom_style>`
command.

All particles in the group must be finite-size spheres.  They cannot
be point particles.

Related commands
""""""""""""""""

:doc:`fix viscous <fix_viscous>`, :doc:`fix damping/cundall <fix_damping_cundall>`

Default
"""""""

none
