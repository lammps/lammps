.. index:: fix nvk

fix nvk command
===============

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nvk

* ID, group-ID are documented in :doc:`fix <fix>` command
* nvk = style name of this fix command

Examples
""""""""

.. parsed-literal::

   fix 1 all nvk

Description
"""""""""""

Perform constant kinetic energy integration using the Gaussian
thermostat to update position and velocity for atoms in the group each
timestep.  V is volume; K is kinetic energy. This creates a system
trajectory consistent with the isokinetic ensemble.

The equations of motion used are those of Minary et al in
:ref:`(Minary) <nvk-Minary>`, a variant of those initially given by Zhang in
:ref:`(Zhang) <nvk-Zhang>`.

The kinetic energy will be held constant at its value given when fix
nvk is initiated. If a different kinetic energy is desired, the
:doc:`velocity <velocity>` command should be used to change the kinetic
energy prior to this fix.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

The Gaussian thermostat only works when it is applied to all atoms in
the simulation box. Therefore, the group must be set to all.

This fix has not yet been implemented to work with the RESPA integrator.

This fix is part of the USER-MISC package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

**Related commands:** none

**Default:** none

----------

.. _nvk-Minary:

**(Minary)** Minary, Martyna, and Tuckerman, J Chem Phys, 18, 2510 (2003).

.. _nvk-Zhang:

**(Zhang)** Zhang, J Chem Phys, 106, 6102 (1997).
