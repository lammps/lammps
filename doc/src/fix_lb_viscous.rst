.. index:: fix lb/viscous

fix lb/viscous command
======================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID lb/viscous

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/viscous = style name of this fix command

Examples
""""""""

fix 1 flow lb/viscous

Description
"""""""""""

This fix is similar to the :doc:`fix viscous <fix_viscous>` command, and
is to be used in place of that command when a lattice-Boltzmann fluid
is present, and the user wishes to integrate the particle motion using
one of the built in LAMMPS integrators.

This fix adds a force, F = - Gamma\*(velocity-fluid\_velocity), to each
atom, where Gamma is the force coupling constant described in the :doc:`fix lb/fluid <fix_lb_fluid>` command (which applies an equal and
opposite force to the fluid).

.. note::

   This fix should only be used in conjunction with one of the
   built in LAMMPS integrators; it should not be used with the :doc:`fix lb/pc <fix_lb_pc>` or :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>` integrators, which
   already include the hydrodynamic forces.  These latter fixes should
   only be used if the force coupling constant has been set by the user
   (instead of using the default value); if the default force coupling
   value is used, then this fix provides the only method for adding the
   hydrodynamic forces to the particles.


----------


For further details, as well as descriptions and results of several
test runs, see :ref:`Mackay et al. <Mackay3>`.  Please include a citation to
this paper if this fix is used in work contributing to published
research.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

As described in the :doc:`fix viscous <fix_viscous>` documentation:

"No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  This fix should only
be used with damped dynamics minimizers that allow for
non-conservative forces.  See the :doc:`min_style <min_style>` command
for details."

Restrictions
""""""""""""


This fix is part of the USER-LB package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Can only be used if a lattice-Boltzmann fluid has been created via the
:doc:`fix lb/fluid <fix_lb_fluid>` command, and must come after this
command.

This fix should not be used if either the :doc:`fix lb/pc <fix_lb_pc>`
or :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>` integrator is
used.

Related commands
""""""""""""""""

:doc:`fix lb/fluid <fix_lb_fluid>`, :doc:`fix lb/pc <fix_lb_pc>`, :doc:`fix lb/rigid/pc/sphere <fix_lb_rigid_pc_sphere>`

**Default:** none


----------


.. _Mackay3:



**(Mackay et al.)** Mackay, F. E., Ollila, S.T.T., and Denniston, C., Hydrodynamic Forces Implemented into LAMMPS through a lattice-Boltzmann fluid, Computer Physics Communications 184 (2013) 2021-2031.


