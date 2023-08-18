.. index:: fix damping/cundall

fix damping/cundall command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID damping/cundall gamma_l gamma_a keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* damping/cundall = style name of this fix command
* gamma_l = linear damping coefficient (dimensionless)
* gamma_a = angular damping coefficient (dimensionless)
* zero or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *scale*
       *scale* values = *type ratio* or *v_name*
         type = atom type (1-N)
         ratio = factor to scale the damping coefficients by
         v_name = reference to atom style variable *name*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all damping/cundall 0.8 0.8
   fix 1 all damping/cundall 0.8 0.5 scale 3 2.5
   fix a all damping/cundall 0.8 0.5 scale v_radscale

Description
"""""""""""

Add damping force and torque to finite-size spherical particles in the group
following the model of :ref:`Cundall, 1987 <Cundall1987>`, as implemented in other
granular physics code (e.g., :ref:`Yade-DEM <YadeDEM>`, :ref:`PFC <PFC>`).

The damping is constructed to always have negative mechanical power with respect
to the current velocity/angular velocity to ensure dissipation of kinetic energy.
If used without additional thermostatting (to add kinetic energy to the system),
it has the effect of slowly (or rapidly) freezing the system; hence it can also
be used as a simple energy minimization technique.

The magnitude of the damping force/torque :math:`F_d`/:math:`T_d` is a fraction
:math:`\gamma \in [0;1]` of the current force/torque :math:`F`/:math:`T` on the
particle. Damping is applied component-by-component in each direction
:math:`k\in\{x, y, z\}`:

.. math::

   {F_d}_k = - \gamma_l \, F_k \, \mathrm{sign}(F_k v_k)

.. math::

   {T_d}_k = - \gamma_a \, T_k \, \mathrm{sign}(T_k \omega_k)

The larger the coefficients, the faster the kinetic energy is reduced.

If the optional keyword *scale* is used, :math:`\gamma_l` and :math:`\gamma_a`
can be scaled up or down by the specified factor for atoms.  This factor can be
set for different atom types and thus the *scale* keyword used multiple times
followed by the atom type and the associated scale factor.  Alternately the
scaling factor can be computed for each atom (e.g. based on its radius) by
using an :doc:`atom-style variable <variable>`.

.. Note::

  The damping force/torque is computed based on the force/torque at the moment
  this fix is invoked. Any force/torque added after this fix, e.g., by
  :doc:`fix addforce <fix_addforce>` or :doc:`fix addtorque <fix_addtorque>`
  will not be damped. When performing simulations with gravity, invoking
  :doc:`fix gravity <fix_gravity>` after this fix will maintain the specified
  gravitational acceleration.

.. Note::

  This scheme is dependent on the coordinates system and does not correspond to
  realistic physical processes. It is constructed for numerical convenience and
  efficacy.

This non-viscous damping presents the following advantages:

1. damping is independent of velocity, equally damping regions with distinct natural frequencies,
2. damping affects acceleration and vanishes for steady uniform motion of the particles,
3. damping parameter :math:`\gamma` is dimensionless and does not require scaling.

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
integrator the fix is modifying forces/torques. Default is the outermost level.

The forces/torques due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  This fix should only
be used with damped dynamics minimizers that allow for
non-conservative forces.  See the :doc:`min_style <min_style>` command
for details.

Restrictions
""""""""""""

This fix is part of the GRANULAR package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

This fix requires that atoms store torque and a radius as defined by the
:doc:`atom_style sphere <atom_style>` command.

Related commands
""""""""""""""""

:doc:`fix viscous <fix_viscous>`, :doc:`fix viscous/sphere <fix_viscous_sphere>`

Default
"""""""

none

References
""""""""""

.. _Cundall1987:

**(Cundall, 1987)** Cundall, P. A. Distinct Element Models of Rock and Soil
Structure, in Analytical and Computational Methods in Engineering Rock
Mechanics, Ch. 4, pp. 129-163. E. T. Brown, ed. London: Allen & Unwin., 1987.

.. _PFC:

**(PFC)** PFC Particle Flow Code 6.0 Documentation. Itasca Consulting Group.

.. _YadeDEM:

**(Yade-DEM)** V. Smilauer et al. (2021), Yade Documentation 3rd ed.
The Yade Project. DOI:10.5281/zenodo.5705394 (https://yade-dem.org/doc/)
