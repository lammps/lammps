.. index:: fix lb/rigid/pc/sphere

fix lb/rigid/pc/sphere command
==============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID lb/rigid/pc/sphere bodystyle args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* lb/rigid/pc/sphere = style name of this fix command
* bodystyle = *single* or *molecule* or *group*
  
  .. parsed-literal::
  
       *single* args = none
       *molecule* args = none
       *group* args = N groupID1 groupID2 ...
         N = # of groups

* zero or more keyword/value pairs may be appended
* keyword = *force* or *torque* or *innerNodes*
  
  .. parsed-literal::
  
       *force* values = M xflag yflag zflag
         M = which rigid body from 1-Nbody (see asterisk form below)
         xflag,yflag,zflag = off/on if component of center-of-mass force is active
       *torque* values = M xflag yflag zflag
         M = which rigid body from 1-Nbody (see asterisk form below)
         xflag,yflag,zflag = off/on if component of center-of-mass torque is active
       *innerNodes* values = innergroup-ID
         innergroup-ID = ID of the atom group which does not experience a hydrodynamic force from the lattice-Boltzmann fluid



Examples
""""""""


.. parsed-literal::

   fix 1 spheres lb/rigid/pc/sphere
   fix 1 all lb/rigid/pc/sphere force 1 0 0 innerNodes ForceAtoms

Description
"""""""""""

This fix is based on the :doc:`fix rigid <fix_rigid>` command, and was
created to be used in place of that fix, to integrate the equations of
motion of spherical rigid bodies when a lattice-Boltzmann fluid is
present with a user-specified value of the force-coupling constant.
The fix uses the integration algorithm described in :ref:`Mackay et al. <Mackay>` to update the positions, velocities, and orientations of
a set of spherical rigid bodies experiencing velocity dependent
hydrodynamic forces.  The spherical bodies are assumed to rotate as
solid, uniform density spheres, with moments of inertia calculated
using the combined sum of the masses of all the constituent particles
(which are assumed to be point particles).


----------


By default, all of the atoms that this fix acts on experience a
hydrodynamic force due to the presence of the lattice-Boltzmann fluid.
However, the *innerNodes* keyword allows the user to specify atoms
belonging to a rigid object which do not interact with the
lattice-Boltzmann fluid (i.e. these atoms do not feel a hydrodynamic
force from the lattice-Boltzmann fluid).  This can be used to
distinguish between atoms on the surface of a non-porous object, and
those on the inside.

This feature can be used, for example, when implementing a hard sphere
interaction between two spherical objects.  Instead of interactions
occurring between the particles on the surfaces of the two spheres, it
is desirable simply to place an atom at the center of each sphere,
which does not contribute to the hydrodynamic force, and have these
central atoms interact with one another.


----------


Apart from the features described above, this fix is very similar to
the rigid fix (although it includes fewer optional arguments, and
assumes the constituent atoms are point particles); see
:doc:`fix rigid <fix_rigid>` for a complete documentation.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about the *rigid* and *rigid/nve* fixes are written to
:doc:`binary restart files <restart>`.

Similar to the :doc:`fix rigid <fix_rigid>` command: The rigid fix
computes a global scalar which can be accessed by various :doc:`output commands <Howto_output>`.  The scalar value calculated by these
fixes is "intensive".  The scalar is the current temperature of the
collection of rigid bodies.  This is averaged over all rigid bodies
and their translational and rotational degrees of freedom.  The
translational energy of a rigid body is 1/2 m v\^2, where m = total
mass of the body and v = the velocity of its center of mass.  The
rotational energy of a rigid body is 1/2 I w\^2, where I = the moment
of inertia tensor of the body and w = its angular velocity.  Degrees
of freedom constrained by the *force* and *torque* keywords are
removed from this calculation.

All of these fixes compute a global array of values which can be
accessed by various :doc:`output commands <Howto_output>`.  The number
of rows in the array is equal to the number of rigid bodies.  The
number of columns is 15.  Thus for each rigid body, 15 values are
stored: the xyz coords of the center of mass (COM), the xyz components
of the COM velocity, the xyz components of the force acting on the
COM, the xyz components of the torque acting on the COM, and the xyz
image flags of the COM, which have the same meaning as image flags for
atom positions (see the "dump" command).  The force and torque values
in the array are not affected by the *force* and *torque* keywords in
the fix rigid command; they reflect values before any changes are made
by those keywords.

The ordering of the rigid bodies (by row in the array) is as follows.
For the *single* keyword there is just one rigid body.  For the
*molecule* keyword, the bodies are ordered by ascending molecule ID.
For the *group* keyword, the list of group IDs determines the ordering
of bodies.

The array values calculated by these fixes are "intensive", meaning
they are independent of the number of atoms in the simulation.

No parameter of these fixes can be used with the *start/stop* keywords
of the :doc:`run <run>` command.  These fixes are not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""


This fix is part of the USER-LB package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Can only be used if a lattice-Boltzmann fluid has been created via the
:doc:`fix lb/fluid <fix_lb_fluid>` command, and must come after this
command.  Should only be used if the force coupling constant used in
:doc:`fix lb/fluid <fix_lb_fluid>` has been set by the user; this
integration fix cannot be used if the force coupling constant is set
by default.

Related commands
""""""""""""""""

:doc:`fix lb/fluid <fix_lb_fluid>`, :doc:`fix lb/pc <fix_lb_pc>`

Default
"""""""

The defaults are force \* on on on, and torque \* on on on.


----------


.. _Mackay:



**(Mackay et al.)** Mackay, F. E., Ollila, S.T.T., and Denniston, C., Hydrodynamic Forces Implemented into LAMMPS through a lattice-Boltzmann fluid, Computer Physics Communications 184 (2013) 2021-2031.
