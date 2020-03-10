.. index:: fix property/atom

fix property/atom command
=========================

fix property/atom/kk command
============================

Syntax
""""""


.. parsed-literal::

   fix ID group-ID property/atom vec1 vec2 ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* property/atom = style name of this fix command
* vec1,vec2,... = *mol* or *q* or *rmass* or *i\_name* or *d\_name*

  .. parsed-literal::

       *mol* = molecule IDs
       *q* = charge
       *rmass* = per-atom mass
       *i_name* = new integer vector referenced by name
       *d_name* = new floating-point vector referenced by name

* zero of more keyword/value pairs may be appended
* keyword = *ghost*

  .. parsed-literal::

       *ghost* value = *no* or *yes* for whether ghost atom info in communicated



Examples
""""""""


.. parsed-literal::

   fix 1 all property/atom mol
   fix 1 all property/atom i_myflag1 i_myflag2
   fix 1 all property/atom d_sx d_sy d_sz

Description
"""""""""""

Create one or more additional per-atom vectors to store information
about atoms and to use during a simulation.  The specified *group-ID*
is ignored by this fix.

The atom style used for a simulation defines a set of per-atom
properties, as explained on the :doc:`atom_style <atom_style>` and
:doc:`read_data <read_data>` doc pages.  The latter command allows these
properties to be defined for each atom in the system when a data file
is read.  This fix will augment the set of properties with new custom
ones. This can be useful in several scenarios.

If the atom style does not define molecule IDs, per-atom charge,
or per-atom mass, they can be added using the *mol*\ , *q* or *rmass*
keywords.  This can be useful, e.g, to define "molecules" to use as
rigid bodies with the :doc:`fix rigid <fix_rigid>` command, or just to
carry around an extra flag with the atoms (stored as a molecule ID)
that can be used to group atoms without having to use the group
command (which is limited to a total of 32 groups including *all*\ ).

Another application would be to use the *rmass* flag in order to have
per-atom masses instead of per-type masses, for example this can be
useful to study isotope effects with partial isotope substitution.
Please :ref:`see below <isotopes>` for an example of simulating a mixture
of light and heavy water with the TIP4P water potential.

An alternative to using fix *property/atom* in these ways is to
use an atom style that does define molecule IDs or charge or per-atom
mass (indirectly via diameter and density) or to use a hybrid atom
style that combines two or more atom styles
to provide the union of all atom properties. However, this has two
practical drawbacks:  first, it typically necessitates changing the
format of the data file, which can be tedious for large systems;
and second, it may define additional properties that are not needed
such as bond lists, which has some overhead when there are no bonds.

In the future, we may add additional per-atom properties similar to
*mol*\ , *q* or *rmass*\ , which "turn-on" specific properties defined
by some atom styles, so they can be used by atom styles that do not
define them.

More generally, the *i\_name* and *d\_name* vectors allow one or more
new custom per-atom properties to be defined.  Each name must be
unique and can use alphanumeric or underscore characters.  These
vectors can store whatever values you decide are useful in your
simulation.  As explained below there are several ways to initialize
and access and output these values, both via input script commands and
in new code that you add to LAMMPS.

This is effectively a simple way to add per-atom properties to a model
without needing to write code for a new :doc:`atom style <atom_style>`
that defines the properties.  Note however that implementing a new
atom style allows new atom properties to be more tightly and
seamlessly integrated with the rest of the code.

The new atom properties encode values that migrate with atoms to new
processors and are written to restart files.  If you want the new
properties to also be defined for ghost atoms, then use the *ghost*
keyword with a value of *yes*\ .  This will invoke extra communication
when ghost atoms are created (at every re-neighboring) to insure the
new properties are also defined for the ghost atoms.

.. note::

   If you use this command with the *mol*\ , *q* or *rmass* vectors,
   then you most likely want to set *ghost* yes, since these properties
   are stored with ghost atoms if you use an :doc:`atom_style <atom_style>`
   that defines them, and many LAMMPS operations that use molecule IDs or
   charge, such as neighbor lists and pair styles, will expect ghost
   atoms to have these values.  LAMMPS will issue a warning it you define
   those vectors but do not set *ghost* yes.

.. note::

   The properties for ghost atoms are not updated every timestep,
   but only once every few steps when neighbor lists are re-built.  Thus
   the *ghost* keyword is suitable for static properties, like molecule
   IDs, but not for dynamic properties that change every step.  For the
   latter, the code you add to LAMMPS to change the properties will also
   need to communicate their new values to/from ghost atoms, an operation
   that can be invoked from within a :doc:`pair style <pair_style>` or
   :doc:`fix <fix>` or :doc:`compute <compute>` that you write.

.. note::

   If this fix is defined **after** the simulation box is created,
   a 'run 0' command should be issued to properly initialize the storage
   created by this fix.


----------


This fix is one of a small number that can be defined in an input
script before the simulation box is created or atoms are defined.
This is so it can be used with the :doc:`read_data <read_data>` command
as described below.

Per-atom properties that are defined by the :doc:`atom style <atom_style>` are initialized when atoms are created, e.g. by
the :doc:`read_data <read_data>` or :doc:`create_atoms <create_atoms>`
commands.  The per-atom properties defined by this fix are not.  So
you need to initialize them explicitly.  This can be done by the
:doc:`read_data <read_data>` command, using its *fix* keyword and
passing it the fix-ID of this fix.

Thus these commands:


.. parsed-literal::

   fix prop all property/atom mol d_flag
   read_data data.txt fix prop NULL Molecules

would allow a data file to have a section like this:


.. parsed-literal::

   Molecules

   1 4 1.5
   2 4 3.0
   3 10 1.0
   4 10 1.0
   5 10 1.0
   ...
   N 763 4.5

where N is the number of atoms, and the first field on each line is
the atom-ID, followed by a molecule-ID and a floating point value that
will be stored in a new property called "flag".  Note that the list of
per-atom properties can be in any order.

Another way of initializing the new properties is via the
:doc:`set <set>` command.  For example, if you wanted molecules
defined for every set of 10 atoms, based on their atom-IDs,
these commands could be used:


.. parsed-literal::

   fix prop all property/atom mol
   variable cluster atom ((id-1)/10)+1
   set atom \* mol v_cluster

The :doc:`atom-style variable <variable>` will create values for atoms
with IDs 31,32,33,...40 that are 4.0,4.1,4.2,...,4.9.  When the
:doc:`set <set>` commands assigns them to the molecule ID for each atom,
they will be truncated to an integer value, so atoms 31-40 will all be
assigned a molecule ID of 4.

Note that :doc:`atomfile-style variables <variable>` can also be used in
place of atom-style variables, which means in this case that the
molecule IDs could be read-in from a separate file and assigned by the
:doc:`set <set>` command.  This allows you to initialize new per-atom
properties in a completely general fashion.


----------


For new atom properties specified as *i\_name* or *d\_name*, the
:doc:`compute property/atom <compute_property_atom>` command can access
their values.  This means that the values can be output via the :doc:`dump custom <dump>` command, accessed by fixes like :doc:`fix ave/atom <fix_ave_atom>`, accessed by other computes like :doc:`compute reduce <compute_reduce>`, or used in :doc:`atom-style variables <variable>`.

For example, these commands will output two new properties to a custom
dump file:


.. parsed-literal::

   fix prop all property/atom i_flag1 d_flag2
   compute 1 all property/atom i_flag1 d_flag2
   dump 1 all custom 100 tmp.dump id x y z c_1[1] c_1[2]


----------


If you wish to add new :doc:`pair styles <pair_style>`,
:doc:`fixes <fix>`, or :doc:`computes <compute>` that use the per-atom
properties defined by this fix, see the :doc:`Modify atom <Modify_atom>`
doc page which has details on how the properties can be accessed from
added classes.


----------


.. _isotopes:



Example for using per-atom masses with TIP4P water to
study isotope effects. When setting up simulations with the :doc:`TIP4P pair styles <Howto_tip4p>` for water, you have to provide exactly
one atom type each to identify the water oxygen and hydrogen
atoms. Since the atom mass is normally tied to the atom type, this
makes it impossible to study multiple isotopes in the same simulation.
With *fix property/atom rmass* however, the per-type masses are
replaced by per-atom masses. Asumming you have a working input deck
for regular TIP4P water, where water oxygen is atom type 1 and water
hydrogen is atom type 2, the following lines of input script convert
this to using per-atom masses:


.. parsed-literal::

   fix Isotopes all property/atom rmass ghost yes
   set type 1 mass 15.9994
   set type 2 mass 1.008

When writing out the system data with the :doc:`write_data <write_data>`
command, there will be a new section named with the fix-ID
(i.e. *Isotopes* in this case). Alternatively, you can take an
existing data file and just add this *Isotopes* section with
one line per atom containing atom-ID and mass. Either way, the
extended data file can be read back with:


.. parsed-literal::

   fix Isotopes all property/atom rmass ghost yes
   read_data tip4p-isotopes.data fix Isotopes NULL Isotopes

Please note that the first *Isotopes* refers to the fix-ID
and the second to the name of the section. The following input
script code will now change the first 100 water molecules in this
example to heavy water:


.. parsed-literal::

   group hwat id 2:300:3
   group hwat id 3:300:3
   set group hwat mass 2.0141018


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


**Restart, fix\_modify, output, run start/stop, minimize info:**

This fix writes the per-atom values it stores to :doc:`binary restart files <restart>`, so that the values can be restored when a
simulation is restarted.  See the :doc:`read_restart <read_restart>`
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for
access by various :doc:`output commands <Howto_output>`.  No parameter
of this fix can be used with the *start/stop* keywords of the
:doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`read_data <read_data>`, :doc:`set <set>`, :doc:`compute property/atom <compute_property_atom>`

Default
"""""""

The default keyword values are ghost = no.
