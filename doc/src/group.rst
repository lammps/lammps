.. index:: group

group command
=============

Syntax
""""""

.. parsed-literal::

   group ID style args

* ID = user-defined name of the group
* style = *delete* or *clear* or *empty* or *region* or         *type* or *id* or *molecule* or *variable* or         *include* or *subtract* or *union* or *intersect* or         *dynamic* or *static*

  .. parsed-literal::

       *delete* = no args
       *clear* = no args
       *empty* = no args
       *region* args = region-ID
       *type* or *id* or *molecule*
         args = list of one or more atom types, atom IDs, or molecule IDs
           any entry in list can be a sequence formatted as A:B or A:B:C where
           A = starting index, B = ending index,
           C = increment between indices, 1 if not specified
         args = logical value
           logical = "<" or "<=" or ">" or ">=" or "==" or "!="
           value = an atom type or atom ID or molecule ID (depending on *style*\ )
         args = logical value1 value2
           logical = "<>"
           value1,value2 = atom types or atom IDs or molecule IDs (depending on *style*\ )
       *variable* args = variable-name
       *include* args = molecule
         molecule = add atoms to group with same molecule ID as atoms already in group
       *subtract* args = two or more group IDs
       *union* args = one or more group IDs
       *intersect* args = two or more group IDs
       *dynamic* args = parent-ID keyword value ...
         one or more keyword/value pairs may be appended
         keyword = *region* or *var* or *every*
           *region* value = region-ID
           *var* value = name of variable
           *property* value = name of per-atom property
           *every* value = N = update group every this many timesteps
       *static* = no args

Examples
""""""""

.. code-block:: LAMMPS

   group edge region regstrip
   group water type 3 4
   group sub id 10 25 50
   group sub id 10 25 50 500:1000
   group sub id 100:10000:10
   group sub id <= 150
   group polyA molecule <> 50 250
   group hienergy variable eng
   group hienergy include molecule
   group boundary subtract all a2 a3
   group boundary union lower upper
   group boundary intersect upper flow
   group boundary delete
   group mine dynamic all region myRegion every 100

Description
"""""""""""

Identify a collection of atoms as belonging to a group.  The group ID
can then be used in other commands such as :doc:`fix <fix>`,
:doc:`compute <compute>`, :doc:`dump <dump>`, or :doc:`velocity <velocity>`
to act on those atoms together.

If the group ID already exists, the group command adds the specified
atoms to the group.

.. note::

   By default groups are static, meaning the atoms are permanently
   assigned to the group.  For example, if the *region* style is used to
   assign atoms to a group, the atoms will remain in the group even if
   they later move out of the region.  As explained below, the *dynamic*
   style can be used to make a group dynamic so that a periodic
   determination is made as to which atoms are in the group.  Since many
   LAMMPS commands operate on groups of atoms, you should think carefully
   about whether making a group dynamic makes sense for your model.

A group with the ID *all* is predefined.  All atoms belong to this
group.  This group cannot be deleted, or made dynamic.

The *delete* style removes the named group and un-assigns all atoms
that were assigned to that group.  Since there is a restriction (see
below) that no more than 32 groups can be defined at any time, the
*delete* style allows you to remove groups that are no longer needed,
so that more can be specified.  You cannot delete a group if it has
been used to define a current :doc:`fix <fix>` or :doc:`compute <compute>`
or :doc:`dump <dump>`.

The *clear* style un-assigns all atoms that were assigned to that
group.  This may be dangerous to do during a simulation run,
e.g. using the :doc:`run every <run>` command if a fix or compute or
other operation expects the atoms in the group to remain constant, but
LAMMPS does not check for this.

The *empty* style creates an empty group, which is useful for commands
like :doc:`fix gcmc <fix_gcmc>` or with complex scripts that add atoms
to a group.

The *region* style puts all atoms in the region volume into the group.
Note that this is a static one-time assignment.  The atoms remain
assigned (or not assigned) to the group even in they later move out of
the region volume.

The *type*\ , *id*\ , and *molecule* styles put all atoms with the
specified atom types, atom IDs, or molecule IDs into the group.  These
3 styles can use arguments specified in one of two formats.

The first format is a list of values (types or IDs).  For example, the
second command in the examples above puts all atoms of type 3 or 4 into
the group named *water*\ .  Each entry in the list can be a
colon-separated sequence A:B or A:B:C, as in two of the examples
above.  A "sequence" generates a sequence of values (types or IDs),
with an optional increment.  The first example with 500:1000 has the
default increment of 1 and would add all atom IDs from 500 to 1000
(inclusive) to the group sub, along with 10,25,50 since they also
appear in the list of values.  The second example with 100:10000:10
uses an increment of 10 and would thus would add atoms IDs
100,110,120, ... 9990,10000 to the group sub.

The second format is a *logical* followed by one or two values (type
or ID).  The 7 valid logicals are listed above.  All the logicals
except <> take a single argument.  The third example above adds all
atoms with IDs from 1 to 150 to the group named *sub*\ .  The logical <>
means "between" and takes 2 arguments.  The fourth example above adds all
atoms belonging to molecules with IDs from 50 to 250 (inclusive) to
the group named polyA.

The *variable* style evaluates a variable to determine which atoms to
add to the group.  It must be an :doc:`atom-style variable <variable>`
previously defined in the input script.  If the variable evaluates
to a non-zero value for a particular atom, then that atom is added
to the specified group.

Atom-style variables can specify formulas that include thermodynamic
quantities, per-atom values such as atom coordinates, or per-atom
quantities calculated by computes, fixes, or other variables.  They
can also include Boolean logic where 2 numeric values are compared to
yield a 1 or 0 (effectively a true or false).  Thus using the
*variable* style, is a general way to flag specific atoms to include
or exclude from a group.

For example, these lines define a variable "eatom" that calculates the
potential energy of each atom and includes it in the group if its
potential energy is above the threshold value -3.0.

.. code-block:: LAMMPS

   compute         1 all pe/atom
   compute         2 all reduce sum c_1
   thermo_style    custom step temp pe c_2
   run             0

   variable        eatom atom "c_1 > -3.0"
   group           hienergy variable eatom

Note that these lines

.. code-block:: LAMMPS

   compute         2 all reduce sum c_1
   thermo_style    custom step temp pe c_2
   run             0

are necessary to insure that the "eatom" variable is current when the
group command invokes it.  Because the eatom variable computes the
per-atom energy via the pe/atom compute, it will only be current if a
run has been performed which evaluated pairwise energies, and the
pe/atom compute was actually invoked during the run.  Printing the
thermodynamic info for compute 2 insures that this is the case, since
it sums the pe/atom compute values (in the reduce compute) to output
them to the screen.  See the "Variable Accuracy" section of the
:doc:`variable <variable>` doc page for more details on insuring that
variables are current when they are evaluated between runs.

The *include* style with its arg *molecule* adds atoms to a group that
have the same molecule ID as atoms already in the group.  The molecule
ID = 0 is ignored in this operation, since it is assumed to flag
isolated atoms that are not part of molecules.  An example of where
this operation is useful is if the *region* style has been used
previously to add atoms to a group that are within a geometric region.
If molecules straddle the region boundary, then atoms outside the
region that are part of molecules with atoms inside the region will
not be in the group.  Using the group command a second time with *include
molecule* will add those atoms that are outside the region to the
group.

.. note::

   The *include molecule* operation is relatively expensive in a
   parallel sense.  This is because it requires communication of relevant
   molecule IDs between all the processors and each processor to loop
   over its atoms once per processor, to compare its atoms to the list of
   molecule IDs from every other processor.  Hence it scales as N, rather
   than N/P as most of the group operations do, where N is the number of
   atoms, and P is the number of processors.

The *subtract* style takes a list of two or more existing group names
as arguments.  All atoms that belong to the first group, but not to any
of the other groups are added to the specified group.

The *union* style takes a list of one or more existing group names as
arguments.  All atoms that belong to any of the listed groups are
added to the specified group.

The *intersect* style takes a list of two or more existing group names
as arguments.  Atoms that belong to every one of the listed groups are
added to the specified group.

----------

The *dynamic* style flags an existing or new group as dynamic.  This
means atoms will be (re)assigned to the group periodically as a
simulation runs.  This is in contrast to static groups where atoms are
permanently assigned to the group.  The way the assignment occurs is
as follows.  Only atoms in the group specified as the parent group via
the parent-ID are assigned to the dynamic group before the following
conditions are applied.  If the *region* keyword is used, atoms not in
the specified region are removed from the dynamic group.  If the *var*
keyword is used, the variable name must be an atom-style or
atomfile-style variable.  The variable is evaluated and atoms whose
per-atom values are 0.0, are removed from the dynamic group. If the *property*
keyword is used, the per-atom property name must be a previously defined
per-atom property.  The per-atom property is evaluated and atoms whose
values are 0.0 are removed from the dynamic group.

The assignment of atoms to a dynamic group is done at the beginning of
each run and on every timestep that is a multiple of *N*\ , which is the
argument for the *every* keyword (N = 1 is the default).  For an
energy minimization, via the :doc:`minimize <minimize>` command, an
assignment is made at the beginning of the minimization, but not
during the iterations of the minimizer.

The point in the timestep at which atoms are assigned to a dynamic
group is after the initial stage of velocity Verlet time integration
has been performed, and before neighbor lists or forces are computed.
This is the point in the timestep where atom positions have just
changed due to the time integration, so the region criterion should be
accurate, if applied.

.. note::

   If the *region* keyword is used to determine what atoms are in
   the dynamic group, atoms can move outside of the simulation box
   between reneighboring events.  Thus if you want to include all atoms
   on the left side of the simulation box, you probably want to set the
   left boundary of the region to be outside the simulation box by some
   reasonable amount (e.g. up to the cutoff of the potential), else they
   may be excluded from the dynamic region.

Here is an example of using a dynamic group to shrink the set of atoms
being integrated by using a spherical region with a variable radius
(shrinking from 18 to 5 over the course of the run).  This could be
used to model a quench of the system, freezing atoms outside the
shrinking sphere, then converting the remaining atoms to a static
group and running further.

.. code-block:: LAMMPS

   variable        nsteps equal 5000
   variable        rad equal 18-(step/v_nsteps)\*(18-5)
   region          ss sphere 20 20 0 v_rad
   group           mobile dynamic all region ss
   fix             1 mobile nve
   run             ${nsteps}
   group           mobile static
   run             ${nsteps}

.. note::

   All fixes and computes take a group ID as an argument, but they
   do not all allow for use of a dynamic group.  If you get an error
   message that this is not allowed, but feel that it should be for the
   fix or compute in question, then please post your reasoning to the
   LAMMPS mail list and we can change it.

The *static* style removes the setting for a dynamic group, converting
it to a static group (the default).  The atoms in the static group are
those currently in the dynamic group.

----------

Restrictions
""""""""""""

There can be no more than 32 groups defined at one time, including
"all".

The parent group of a dynamic group cannot itself be a dynamic group.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`fix <fix>`, :doc:`region <region>`,
:doc:`velocity <velocity>`

Default
"""""""

All atoms belong to the "all" group.
