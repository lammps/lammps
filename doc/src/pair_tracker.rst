.. index:: pair_style tracker

pair_style tracker command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style tracker fix_ID N keyword values attribute1 attribute2 ...

* fix_ID = ID of associated internal fix to store data
* N = prepare data for output every this many timesteps
* zero or more keywords may be appended
* keyword = *finite* or *time/min* or *type/include*

  .. parsed-literal::

      *finite* value = none
         pair style uses atomic diameters to identify contacts
      *time/min* value = T
         T = minimum number of timesteps of interaction
      *type/include* value = list1 list2
         list1,list2 = separate lists of types (see below)

* one or more attributes may be appended

  .. parsed-literal::

       possible attributes = id1 id2 time/created time/broken time/total
                             r/min r/ave x y z

  .. parsed-literal::

          id1, id2 = IDs of the 2 atoms in each pair interaction
          time/created = the timestep that the 2 atoms began interacting
          time/broken = the timestep that the 2 atoms stopped interacting
          time/total = the total number of timesteps the 2 atoms interacted
          r/min = the minimum radial distance between the 2 atoms during the interaction (distance units)
          r/ave = the average radial distance between the 2 atoms during the interaction (distance units)
          x, y, z = the center of mass position of the 2 atoms when they stopped interacting (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay tracker myfix 1000 id1 id2 type/include 1 * type/include 2 3,4  lj/cut 2.5
   pair_coeff 1 1 tracker 2.0

   pair_style hybrid/overlay tracker myfix 1000 finite x y z time/min 100 granular
   pair_coeff * * tracker

   dump 1 all local 1000 dump.local f_myfix[1] f_myfix[2] f_myfix[3]
   dump_modify 1 write_header no

Description
"""""""""""

Style *tracker* monitors information about pairwise interactions.  It
does not calculate any forces on atoms. :doc:`Pair hybrid/overlay
<pair_hybrid>` can be used to combine this pair style with any other
pair style, as shown in the examples above.

At each timestep, if two neighboring atoms move beyond the interaction
cutoff, pairwise data is processed and transferred to an internal fix
labeled *fix_ID*. This allows the local data to be accessed by other
LAMMPS commands. Additional
filters can be applied using the *time/min* or *type/include* keywords
described below.  Note that this is the interaction cutoff defined by
this pair style, not the short-range cutoff defined by the pair style
that is calculating forces on atoms.

Following any optional keyword/value arguments, a list of one or more
attributes is specified.  These include the IDs of the two atoms in
the pair.  The other attributes for the pair of atoms are the
duration of time they were "interacting" or at the point in time they
started or stopped interacting.  In this context, "interacting" means
the time window during which the two atoms were closer than the
interaction cutoff distance.  The attributes for time/* refer to
timesteps.

Data is continuously accumulated by the internal fix over intervals of *N*
timesteps. At the end of each interval, all of the saved accumulated
data is deleted to make room for new data. Individual datum may
therefore persist anywhere between *1* to *N* timesteps depending on
when they are saved. This data can be accessed using the *fix_ID* and a
:doc:`dump local <dump>` command. To ensure all data is output,
the dump frequency should correspond to the same interval of *N*
timesteps. A dump frequency of an integer multiple of *N* can be used
to regularly output a sample of the accumulated data.

----------

The following optional keywords may be used.

If the *finite* keyword is not used, the following coefficients must
be defined for each pair of atom types via the :doc:`pair_coeff
<pair_coeff>` command as in the examples above, or in the data file or
restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described
below:

* cutoff (distance units)

If the *finite* keyword is used, there are no additional coefficients
to set for each pair of atom types via the
:doc:`pair_coeff <pair_coeff>` command. Interaction cutoffs are
instead calculated based on the diameter of finite particles. However
you must still use the :doc:`pair_coeff <pair_coeff>` for all atom
types. For example the command

.. code-block:: LAMMPS

   pair_coeff * *

should be used.

The *time/min* keyword sets a minimum amount of time that an
interaction must persist to be included.  This setting can be used to
censor short-lived interactions.

The *type/include* keyword filters interactions based on the types of
the two atoms.  Data is only saved for interactions between atoms
whose two atom types appear in *list1* and *list2*.  Atom type 1 must
be in list1 and atom type 2 in list2.  Or vice versa.

Each type list consists of a series of type ranges separated by
commas.  Each range can be specified as a single numeric value, or a
wildcard asterisk can be used to specify a range of values.  This
takes the form "\*" or "\*n" or "n\*" or "m\*n".  For example, if M =
the number of atom types, then an asterisk with no numeric values
means all types from 1 to M.  A leading asterisk means all types from
1 to n (inclusive).  A trailing asterisk means all types from n to M
(inclusive).  A middle asterisk means all types from m to n
(inclusive).  Note that the *type/include* keyword can be specified
multiple times.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the cutoff coefficient and cutoff
distance for this pair style can be mixed.  The cutoff is always mixed
via a *geometric* rule.  The cutoff is mixed according to the
pair_modify mix value.  The default mix value is *geometric*\ .  See
the "pair_modify" command for details.

This pair style writes its information to :doc:`binary restart files
<restart>`, so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

The :doc:`pair_modify <pair_modify>` shift, table, and tail options
are not relevant for this pair style.

The accumulated data is not written to restart files and should be
output before a restart file is written to avoid missing data.

The internal fix calculates a local vector or local array depending on the
number of input values.  The length of the vector or number of rows in
the array is the number of recorded, lost interactions.  If a single
input is specified, a local vector is produced.  If two or more inputs
are specified, a local array is produced where the number of columns =
the number of inputs.  The vector or array can be accessed by any
command that uses local values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS
output options.

The vector or array will be floating point values that correspond to
the specified attribute.

----------

Restrictions
""""""""""""

This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style is currently incompatible with granular pair styles
that extend beyond the contact (e.g. JKR and DMT).

Related commands
""""""""""""""""

Default
"""""""

none
