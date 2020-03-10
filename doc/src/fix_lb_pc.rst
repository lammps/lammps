.. index:: fix lb/pc

fix lb/pc command
=================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID lb/pc

* ID, group-ID are documented in the :doc:`fix <fix>` command
* lb/pc = style name of this fix command

Examples
""""""""

.. parsed-literal::

   fix 1 all lb/pc

Description
"""""""""""

Update the positions and velocities of the individual particles
described by *group-ID*\ , experiencing velocity-dependent hydrodynamic
forces, using the integration algorithm described in :ref:`Mackay et al. <Mackay1>`.  This integration algorithm should only be used if a
user-specified value for the force-coupling constant used in :doc:`fix lb/fluid <fix_lb_fluid>` has been set; do not use this integration
algorithm if the force coupling constant has been set by default.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the USER-LB package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Can only be used if a lattice-Boltzmann fluid has been created via the
:doc:`fix lb/fluid <fix_lb_fluid>` command, and must come after this
command.

Related commands
""""""""""""""""

:doc:`fix lb/fluid <fix_lb_fluid>` :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>`

**Default:** None.

----------

.. _Mackay1:

**(Mackay et al.)** Mackay, F. E., Ollila, S.T.T., and Denniston, C., Hydrodynamic Forces Implemented into LAMMPS through a lattice-Boltzmann fluid, Computer Physics Communications 184 (2013) 2021-2031.
