.. index:: write_data

write_data command
===================

Syntax
""""""

.. code-block:: LAMMPS

   write_data file keyword value ...

* file = name of data file to write out
* zero or more keyword/value pairs may be appended
* keyword = *pair* or *nocoeff* or *nofix* or *nolabelmap*

  .. parsed-literal::

       *nocoeff* = do not write out force field info
       *nofix* = do not write out extra sections read by fixes
       *nolabelmap* = do not write out type labels
       *triclinic/general = write data file in general triclinic format
       *types* value = *numeric* or *labels*
       *pair* value = *ii* or *ij*
         *ii* = write one line of pair coefficient info per atom type
         *ij* = write one line of pair coefficient info per IJ atom type pair

Examples
""""""""

.. code-block:: LAMMPS

   write_data data.polymer
   write_data data.*
   write_data data.solid triclinic/general

Description
"""""""""""

Write a data file in text format of the current state of the
simulation.  Data files can be read by the :doc:`read data <read_data>`
command to begin a simulation.  The :doc:`read_data <read_data>` command
also describes their format.

Similar to :doc:`dump <dump>` files, the data filename can contain a "\*"
wild-card character.  The "\*" is replaced with the current timestep
value.

.. admonition:: Data in Coeff sections
   :class: note

   The write_data command may not always write all coefficient settings
   to the corresponding Coeff sections of the data file.  This can have
   one of multiple reasons. 1) A few styles may be missing the code that
   would write those sections (if you come across one, please notify
   the LAMMPS developers). 2) Some pair styles require a single pair_coeff
   statement and those are not compatible with data files. 3) The
   default for write_data is to write a PairCoeff section, which has
   only entries for atom types i == j. The remaining coefficients would
   be inferred through the currently selected mixing rule.  If there has
   been a pair_coeff command with i != j, this setting would be lost.
   LAMMPS will detect this and print a warning message unless *pair ij*
   is appended to the write_data command.  This will request writing a
   PairIJCoeff section which has information for all pairs of atom types.
   In cases where the coefficient data in the data file is incomplete,
   you will need to re-specify that information in your input script
   that reads the data file.

Because a data file is in text format, if you use a data file written
out by this command to restart a simulation, the initial state of the
new run will be slightly different than the final state of the old run
(when the file was written) which was represented internally by LAMMPS
in binary format.  A new simulation which reads the data file will
thus typically diverge from a simulation that continued in the
original input script.

If you want to do more exact restarts, using binary files, see the
:doc:`restart <restart>`, :doc:`write_restart <write_restart>`, and
:doc:`read_restart <read_restart>` commands.  You can also convert
binary restart files to text data files, after a simulation has run,
using the :doc:`-r command-line switch <Run_options>`.

.. note::

   Only limited information about a simulation is stored in a data
   file.  For example, no information about atom :doc:`groups <group>` and
   :doc:`fixes <fix>` are stored.  :doc:`Binary restart files <read_restart>`
   store more information.

Bond interactions (angle, etc) that have been turned off by the
:doc:`fix shake <fix_shake>` or :doc:`delete_bonds <delete_bonds>`
command will be written to a data file as if they are turned on.  This
means they will need to be turned off again in a new run after the
data file is read.

Bonds that are broken (e.g. by a bond-breaking potential) are not
written to the data file.  Thus these bonds will not exist when the
data file is read.

----------

Use of the *nocoeff* keyword means no force field parameters are
written to the data file. This can be helpful, for example, if you
want to make significant changes to the force field or if the force
field parameters are read in separately, e.g. from an include file.

Use of the *nofix* keyword means no extra sections read by fixes are
written to the data file (see the *fix* option of the :doc:`read_data
<read_data>` command for details). For example, this option excludes
sections for user-created per-atom properties from :doc:`fix
property/atom <fix_property_atom>`.

The *nolabelmap* and *types* keywords refer to type labels that may be
defined for numeric atom types, bond types, angle types, etc.  The
label map can be defined in two ways, either by the :doc:`labelmap
<labelmap>` command or in data files read by the :doc:`read_data
<read_data>` command which have sections for Atom Type Labels, Bond
Type Labels, Angle Type Labels, etc.  See the :doc:`Howto type labels
<Howto_type_labels>` doc page for the allowed syntax of type labels
and a general discussion of how type labels can be used.

Use of the *nolabelmap* keyword means that even if type labels exist
for a given type-kind (Atoms, Bonds, Angles, etc.), type labels are
not written to the data file.  By default, they are written if they
exist.  A type label must be defined for every numeric type (within a
given type-kind) to be written to the data file.

Use of the *triclinic/general* keyword will output a data file which
specifies a general triclinic simulation box as well as per-atom
quantities consistent with the general triclinic box.  The latter means
that per-atom vectors, such as velocities and dipole moments will be
oriented consistent with the 3d rotation implied by the general
triclinic box (relative to the associated restricted triclinic box).

This option can only be requested if the simulation box was initially
defined to be general triclinic.  If if was and the
*triclinic/general* keyword is not used, then the data file will
specify a restricted triclinic box, since that is the internal format
LAMMPS uses for both general and restricted triclinic simulations.
See the :doc:`Howto triclinic <Howto_triclinic>` doc page for more
explanation of how general triclinic simulation boxes are supported by
LAMMPS.  And see the :doc:`read_data <read_data>` doc page for details
of how the format is altered for general triclinic data files.

The *types* keyword determines how atom types, bond types, angle
types, etc are written into these data file sections: Atoms, Bonds,
Angles, etc.  The default is the *numeric* setting, even if type label
maps exist.  If the *labels* setting is used, type labels will be
written to the data file, if the corresponding label map exists.  Note
that when using *types labels*, the *nolabelmap* keyword cannot be
used.

The *pair* keyword lets you specify in what format the pair
coefficient information is written into the data file.  If the value
is specified as *ii*, then one line per atom type is written, to
specify the coefficients for each of the I=J interactions.  This means
that no cross-interactions for I != J will be specified in the data
file and the pair style will apply its mixing rule, as documented on
individual :doc:`pair_style <pair_style>` doc pages.  Of course this
behavior can be overridden in the input script after reading the data
file, by specifying additional :doc:`pair_coeff <pair_coeff>` commands
for any desired I,J pairs.

If the value is specified as *ij*, then one line of coefficients is
written for all I,J pairs where I <= J.  These coefficients will
include any specific settings made in the input script up to that
point.  The presence of these I != J coefficients in the data file
will effectively turn off the default mixing rule for the pair style.
Again, the coefficient values in the data file can be overridden
in the input script after reading the data file, by specifying
additional :doc:`pair_coeff <pair_coeff>` commands for any desired I,J
pairs.

----------

Restrictions
""""""""""""

This command requires inter-processor communication to migrate atoms
before the data file is written.  This means that your system must be
ready to perform a simulation before using this command (force fields
setup, atom masses initialized, etc).

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`write_restart <write_restart>`

Default
"""""""

The option defaults are pair = ii and types_style = numeric.
