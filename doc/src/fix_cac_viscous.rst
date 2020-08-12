.. index:: fix cac/viscous

fix cac/viscous command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID cac/viscous gamma keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* viscous = style name of this fix command
* gamma = damping coefficient (force/velocity units)
* zero or more keyword/value pairs may be appended
  
  .. parsed-literal::
  
     keyword = *scale*
       *scale* values = type ratio
         type = atom type (1-N)
         ratio = factor to scale the damping coefficient by

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 flow cac/viscous 0.1
   fix 1 damp cac/viscous 0.5 scale 3 2.5

Description
"""""""""""

Add a viscous damping force to atoms and finite elements in the group 
that is proportional to the velocity of the atom or the finite element's nodal velocities. 
The added force can be thought of as a
frictional interaction with implicit solvent, i.e. the no-slip Stokes
drag on a spherical particle.  This can be useful for draining the
kinetic energy from the system in a controlled fashion.
If used without additional thermostatting (to add kinetic
energy to the system), it has the effect of slowly (or rapidly)
freezing the system; hence it can also be used as a simple energy
minimization technique.

The damping force F is given by F = - gamma * velocity.  The larger
the coefficient, the faster the kinetic energy is reduced.  If the
optional keyword *scale* is used, gamma can scaled up or down by the
specified factor for atoms of that type.  It can be used multiple
times to adjust gamma for several atom types.

.. note::

   You should specify gamma in force/velocity units.  This is not
   the same as mass/time units, at least for some of the LAMMPS
   :doc:`units <units>` options like "real" or "metal" that are not
   self-consistent.

In a Brownian dynamics context, gamma = Kb T / D, where Kb =
Boltzmann's constant, T = temperature, and D = particle diffusion
coefficient.  D can be written as Kb T / (3 pi eta d), where eta =
dynamic viscosity of the frictional fluid and d = diameter of
particle.  This means gamma = 3 pi eta d, and thus is proportional to
the viscosity of the fluid and the particle diameter.

In the current implementation, rather than have the user specify a
viscosity, gamma is specified directly in force/velocity units.  If
needed, gamma can be adjusted for atoms of different sizes
(i.e. sigma) by using the *scale* keyword.


----------


**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is modifying forces. Default is the outermost level.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  This fix should only
be used with damped dynamics minimizers that allow for
non-conservative forces.  See the :doc:`min_style <min_style>` command
for details.

Restrictions
""""""""""""

This fix requires a CAC atom style

**Default:** none
