.. index:: fix store/state

fix store/state command
=======================

Syntax
""""""

.. code-block:: LAMMPS

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
* keyword = *com* or *history* or *thresh*

  .. parsed-literal::

       *com* value = *yes* or *no*
       *history* values = Nevery Nrepeat Nfreq
         Nevery = accumulate atom attributes once every this many steps
         Nrepeat = # of times to accumulate atom attributes
         Nfreq = make stored atom attributes (history) available every this many steps
       *thresh* args = variable operator value
         variable = equal style or compatible variable reference
         operator = "<" or "<=" or ">" or ">=" or "==" or "!=" or "\|^"
         value = numeric value to compare to

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all store/state 0 x y z
   fix 1 all store/state 0 xu yu zu com yes
   fix 2 all store/state 1000 vx vy vz
   fix 2 all store/state 0 vx vy vz history 5 100 0
   fix 2 all store/state 0 vx vy vz history 5 20 1000
   fix 1 all store/state 0 x y z thresh v_dist < 5.0

Description
"""""""""""

.. versionchanged:: 4Jul2026

   The *history* keyword was added.

Define a fix that stores attributes for each atom in the group either
once or for multiple recent timesteps by use of the optional *history*
keyword.

If the optional keyword *history* is not used, then atom attributes are
stored once.  If *N* > 0, then the stored attributes will be updated
once every *N* steps.  If *N* = 0, then the attributes are stored when
the fix is defined (see the following Note) and never changed.  The
latter is a way of archiving an atom attribute for future use in a
calculation or output.

.. note::

   Only atom attributes specified by keywords like *xu* or *vy* or
   *radius* are initially stored immediately at the point in your input
   script when the fix is defined.  Attributes specified by a compute,
   fix, or variable are not initially stored until the first run
   following the fix definition begins.  This is because calculating
   those attributes may require quantities that are not defined in
   between runs.

If the optional keyword *history* is used, then the *N* keyword must
be specified as 0, and attributes are stored for multiple timesteps.
This enables calculations involving the recent history of each atom.
The Nevery, Nrepeat, and Nfreq settings determine when the history is
stored and made available for calculations.  There are 2 use cases,
for Nfreq = 0 and Nfreq > 0.  In both cases, attributes are stored
Nrepeat times, once every Nevery steps.  Nevery >= 1, Nrepeat > 1, and
Nfreq = a multiple of Nevery (when Nfreq > 0) are required.

If Nfreq = 0, then the Nrepeat attributes are updated continuously and
the attributes are accessible on any timestep which is a multiple of
Nevery.  For example, this command in the Examples section above:

.. code-block:: LAMMPS

   fix 2 all store/state 0 vx vy vz history 5 100 0

will continuously store attributes once every 5 timesteps repeated 100
times (including the current timestep).  On step 500, each atom stores
attributes for steps 5,10,15, ... 500.  On step 600, each atom stores
attributes for steps 105,105,110, ... 600.  For timesteps < 495, less
than 100 attributes are necessarily stored.

If Nfreq is non-zero, it must be a multiple of Nevery, and the
attributes are accessible only on timesteps which are a multiple of
Nfreq.  This can be useful to limit storage of atom attributes to only
the time windows when the attributes are needed.  For example, this
command in the Examples section above:

.. code-block:: LAMMPS

   fix 2 all store/state 0 vx vy vz history 5 20 1000

will store attributes once every 5 timesteps repeated 20 times
(including the current timestep), but only preceding timesteps which
are multiples of 1000.  On step 1000, each atom stores attributes for
steps 905,910,915, ... 1000.  On step 2000, each atom stores
attributes for steps 1905,1910,1915, ... 2000.  Between timesteps 1000
to 1905, no attributes are stored, which can be more efficient if
those attributes are not needed.

.. warning::

   Specifying a large value for *Nrepeat* may require significant
   extra memory.  Since the attribute values need to persist with each
   atom as it migrates to a new processor, it may also slow down a
   simulation.

----------

The list of possible attributes is the same as that used by the
:doc:`dump custom <dump>` command, which describes their meaning.

The requested values are stored in a per-atom vector or array as
discussed below.  Zeroes are stored for atoms not in the specified
group.

If the *com* keyword is set to *yes* then the *xu*, *yu*, and *zu*
inputs store the position of each atom relative to the center-of-mass
of the group of atoms, instead of storing the absolute position.

.. versionadded:: 10Dec2025

If the *thresh* keyword is used, data is only stored on steps where also
the threshold condition following the keyword is met.  The first
argument *must* be a :doc:`variable <variable>` reference to an
equal-style or compatible variable.  The second argument is a logical
operator and the third argument is a number.  The choice of operators
listed above are the usual comparison operators. The XOR operation
(exclusive or) is also included as "\|\^".  In this context, XOR means
that if either the variable or the value is 0.0 and the other is
non-zero, then the result is "true" and the threshold criterion is met.
Otherwise it is not met.

The *thresh* keyword may be used *exactly once*.  If a more complex
threshold condition is needed or access to data from computes, fixes, or
similar, this can be realized within the expression used when defining
the referenced variable.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

If the optional *history* keyword is not used, this fix writes the
per-atom values it stores to :doc:`binary restart files <restart>`, so
that the values can be restored when a simulation is restarted.  See
the :doc:`read_restart <read_restart>` command for info on how to
re-specify a fix in an input script that reads a restart file, so that
the operation of the fix continues in an uninterrupted fashion.

.. warning::

   To enable the stored per-atom info to be restored from the restart
   file, this fix command has to be specified **exactly** the same way
   as before. LAMMPS will only check whether a fix is of the same
   style and has the same fix ID and in case of a match will then try
   to initialize the fix with the data stored in the binary restart
   file.  If the fix store/state command does not match exactly (for
   example, if the list of values does not match), data can be
   corrupted or LAMMPS may crash.  If the new fix uses the optional
   *history* keyword, the data stored in the restart file will be
   discarded.

If the *history* keyword is used, no information about this fix is
written to :doc:`binary restart files <restart>`.  Instead, you should
simply specify a new fix in the restarted run to begin new
accumulation of per-atom history.

None of the :doc:`fix_modify <fix_modify>` options are relevant to
this fix.

If the optional *history* keyword is not used, this fix produces a
per-atom vector if a single input is specified.  Or a per-atom array
if multiple inputs are specified, where the number of columns for each
atom is the number of inputs.  These can be accessed by various
:doc:`output commands <Howto_output>`.  These per-atom values can be
accessed on any timestep (see the discussion of Noutput below).

If the *history* keyword is used, no per-atom vector or array is
produced.  Instead, the fix provides access to its multiple timesteps
of stored per-atom history via its *extract()* method.
The :doc:`compute property/atom <compute_property_atom>`
command is an example of a command which does this.

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
