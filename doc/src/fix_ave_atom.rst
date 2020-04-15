.. index:: fix ave/atom

fix ave/atom command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID ave/atom Nevery Nrepeat Nfreq value1 value2 ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/atom = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many timesteps
  one or more input values can be listed
* value = x, y, z, vx, vy, vz, fx, fy, fz, c_ID, c_ID[i], f_ID, f_ID[i], v_name

  .. parsed-literal::

       x,y,z,vx,vy,vz,fx,fy,fz = atom attribute (position, velocity, force component)
       c_ID = per-atom vector calculated by a compute with ID
       c_ID[I] = Ith column of per-atom array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = per-atom vector calculated by a fix with ID
       f_ID[I] = Ith column of per-atom array calculated by a fix with ID, I can include wildcard (see below)
       v_name = per-atom vector calculated by an atom-style variable with name

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/atom 1 100 100 vx vy vz
   fix 1 all ave/atom 10 20 1000 c_my_stress[1]
   fix 1 all ave/atom 10 20 1000 c_my_stress[*]

Description
"""""""""""

Use one or more per-atom vectors as inputs every few timesteps, and
average them atom by atom over longer timescales.  The resulting
per-atom averages can be used by other :doc:`output commands <Howto_output>` such as the :doc:`fix ave/chunk <fix_ave_chunk>` or :doc:`dump custom <dump>` commands.

The group specified with the command means only atoms within the group
have their averages computed.  Results are set to 0.0 for atoms not in
the group.

Each input value can be an atom attribute (position, velocity, force
component) or can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an atom-style
:doc:`variable <variable>`.  In the latter cases, the compute, fix, or
variable must produce a per-atom vector, not a global quantity or
local quantity.  If you wish to time-average global quantities from a
compute, fix, or variable, then see the :doc:`fix ave/time <fix_ave_time>` command.

Each per-atom value of each input vector is averaged independently.

:doc:`Computes <compute>` that produce per-atom vectors or arrays are
those which have the word *atom* in their style name.  See the doc
pages for individual :doc:`fixes <fix>` to determine which ones produce
per-atom vectors or arrays.  :doc:`Variables <variable>` of style *atom*
are the only ones that can be used with this fix since they produce
per-atom vectors.

Note that for values from a compute or fix, the bracketed index I can
be specified using a wildcard asterisk with the index to effectively
specify multiple values.  This takes the form "\*" or "\*n" or "n\*" or
"m\*n".  If N = the size of the vector (for *mode* = scalar) or the
number of columns in the array (for *mode* = vector), then an asterisk
with no numeric values means all indices from 1 to N.  A leading
asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to N (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual columns of the array
had been listed one by one.  E.g. these 2 fix ave/atom commands are
equivalent, since the :doc:`compute stress/atom <compute_stress_atom>`
command creates a per-atom array with 6 columns:

.. code-block:: LAMMPS

   compute my_stress all stress/atom NULL
   fix 1 all ave/atom 10 20 1000 c_my_stress[*]
   fix 1 all ave/atom 10 20 1000 c_my_stress[1] c_my_stress[2] &
                                 c_my_stress[3] c_my_stress[4] &
                                 c_my_stress[5] c_my_stress[6]

----------

The *Nevery*\ , *Nrepeat*\ , and *Nfreq* arguments specify on what
timesteps the input values will be used in order to contribute to the
average.  The final averaged quantities are generated on timesteps
that are a multiple of *Nfreq*\ .  The average is over *Nrepeat*
quantities, computed in the preceding portion of the simulation every
*Nevery* timesteps.  *Nfreq* must be a multiple of *Nevery* and
*Nevery* must be non-zero even if *Nrepeat* is 1.  Also, the timesteps
contributing to the average value cannot overlap,
i.e. Nrepeat\*Nevery can not exceed Nfreq.

For example, if Nevery=2, Nrepeat=6, and Nfreq=100, then values on
timesteps 90,92,94,96,98,100 will be used to compute the final average
on timestep 100.  Similarly for timesteps 190,192,194,196,198,200 on
timestep 200, etc.

----------

The atom attribute values (x,y,z,vx,vy,vz,fx,fy,fz) are
self-explanatory.  Note that other atom attributes can be used as
inputs to this fix by using the :doc:`compute property/atom <compute_property_atom>` command and then specifying
an input value from that compute.

.. note::

   The x,y,z attributes are values that are re-wrapped inside the
   periodic box whenever an atom crosses a periodic boundary.  Thus if
   you time average an atom that spends half its time on either side of
   the periodic box, you will get a value in the middle of the box.  If
   this is not what you want, consider averaging unwrapped coordinates,
   which can be provided by the :doc:`compute property/atom <compute_property_atom>` command via its xu,yu,zu
   attributes.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If no bracketed term is
appended, the per-atom vector calculated by the compute is used.  If a
bracketed term containing an index I is appended, the Ith column of
the per-atom array calculated by the compute is used.  Users can also
write code for their own compute styles and :doc:`add them to LAMMPS <Modify>`.  See the discussion above for how I can
be specified with a wildcard asterisk to effectively specify multiple
values.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If no bracketed term is
appended, the per-atom vector calculated by the fix is used.  If a
bracketed term containing an index I is appended, the Ith column of
the per-atom array calculated by the fix is used.  Note that some
fixes only produce their values on certain timesteps, which must be
compatible with *Nevery*\ , else an error will result.  Users can also
write code for their own fix styles and :doc:`add them to LAMMPS <Modify>`.  See the discussion above for how I can be
specified with a wildcard asterisk to effectively specify multiple
values.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script as an :doc:`atom-style variable <variable>` Variables of style *atom* can reference
thermodynamic keywords, or invoke other computes, fixes, or variables
when they are evaluated, so this is a very general means of generating
per-atom quantities to time average.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global scalar or vector quantities are
stored by this fix for access by various :doc:`output commands <Howto_output>`.

This fix produces a per-atom vector or array which can be accessed by
various :doc:`output commands <Howto_output>`.  A vector is produced if
only a single quantity is averaged by this fix.  If two or more
quantities are averaged, then an array of values is produced.  The
per-atom values can only be accessed on timesteps that are multiples
of *Nfreq* since that is when averaging is performed.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix ave/histo <fix_ave_histo>`, :doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`variable <variable>`,

**Default:** none
