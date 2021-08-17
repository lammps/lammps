.. index:: fix store/state

fix store/state command
=======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID store/state N input1 input2 ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* store/state = style name of this fix command
* N = store atom attributes every N steps, N = 0 for initial store only
* input = one or more atom attributes

  .. parsed-literal::

       possible attributes = id, mol, type, mass,
                             x, y, z, xs, ys, zs, xu, yu, zu, xsu, ysu, zsu, ix, iy, iz,
                             vx, vy, vz, fx, fy, fz,
                             q, mux, muy, muz, mu,
                             radius, diameter, omegax, omegay, omegaz,
                             angmomx, angmomy, angmomz, tqx, tqy, tqz,
                             c_ID, c_ID[I], f_ID, f_ID[I], v_name,
                             d_name, i_name, i2_name[I], d2_name[I],

  .. parsed-literal::

           id = atom ID
           mol = molecule ID
           type = atom type
           mass = atom mass
           x,y,z = unscaled atom coordinates
           xs,ys,zs = scaled atom coordinates
           xu,yu,zu = unwrapped atom coordinates
           xsu,ysu,zsu = scaled unwrapped atom coordinates
           ix,iy,iz = box image that the atom is in
           vx,vy,vz = atom velocities
           fx,fy,fz = forces on atoms
           q = atom charge
           mux,muy,muz = orientation of dipolar atom
           mu = magnitued of dipole moment of atom
           radius,diameter = radius.diameter of spherical particle
           omegax,omegay,omegaz = angular velocity of spherical particle
           angmomx,angmomy,angmomz = angular momentum of aspherical particle
           tqx,tqy,tqz = torque on finite-size particles
           *c_ID* = per-atom vector calculated by a compute with ID
           *c_ID[I]* = Ith column of per-atom array calculated by a compute with ID
           *f_ID* = per-atom vector calculated by a fix with ID
           *f_ID[I]* = Ith column of per-atom array calculated by a fix with ID
           *v_name* = per-atom vector calculated by an atom-style variable with name
           *i_name* = custom integer vector with name
           *d_name* = custom floating point vector with name
           *i2_name[I]* = Ith column of custom integer array with name
           *d2_name[I]* = Ith column of custom floating-point array with name

* zero or more keyword/value pairs may be appended
* keyword = *com*

  .. parsed-literal::

       *com* value = *yes* or *no*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all store/state 0 x y z
   fix 1 all store/state 0 xu yu zu com yes
   fix 2 all store/state 1000 vx vy vz

Description
"""""""""""

Define a fix that stores attributes for each atom in the group at the
time the fix is defined.  If *N* is 0, then the values are never
updated, so this is a way of archiving an atom attribute at a given
time for future use in a calculation or output.  See the discussion of
:doc:`output commands <Howto_output>` that take fixes as inputs.

If *N* is not zero, then the attributes will be updated every *N*
steps.

.. note::

   Actually, only atom attributes specified by keywords like *xu*
   or *vy* or *radius* are initially stored immediately at the point in
   your input script when the fix is defined.  Attributes specified by a
   compute, fix, or variable are not initially stored until the first run
   following the fix definition begins.  This is because calculating
   those attributes may require quantities that are not defined in
   between runs.

The list of possible attributes is the same as that used by the
:doc:`dump custom <dump>` command, which describes their meaning.

If the *com* keyword is set to *yes* then the *xu*, *yu*, and *zu*
inputs store the position of each atom relative to the center-of-mass
of the group of atoms, instead of storing the absolute position.

The requested values are stored in a per-atom vector or array as
discussed below.  Zeroes are stored for atoms not in the specified
group.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the per-atom values it stores to :doc:`binary restart
files <restart>`, so that the values can be restored when a simulation
is restarted.  See the :doc:`read_restart <read_restart>` command for
info on how to re-specify a fix in an input script that reads a
restart file, so that the operation of the fix continues in an
uninterrupted fashion.

.. warning::

   When reading data from a restart file, this fix command has to be specified
   **exactly** the same way as before. LAMMPS will only check whether a
   fix is of the same style and has the same fix ID and in case of a match
   will then try to initialize the fix with the data stored in the binary
   restart file.  If the fix store/state command does not match exactly,
   data can be corrupted or LAMMPS may crash.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

If a single input is specified, this fix produces a per-atom vector.
If multiple inputs are specified, a per-atom array is produced where
the number of columns for each atom is the number of inputs.  These
can be accessed by various :doc:`output commands <Howto_output>`.  The
per-atom values be accessed on any timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dump custom <dump>`, :doc:`compute property/atom <compute_property_atom>`,
:doc:`fix property/atom <fix_property_atom>`, :doc:`variable <variable>`

Default
"""""""

The option default is com = no.
