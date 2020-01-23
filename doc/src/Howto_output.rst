Output from LAMMPS (thermo, dumps, computes, fixes, variables)
==============================================================

There are four basic kinds of LAMMPS output:

* :doc:`Thermodynamic output <thermo_style>`, which is a list
  of quantities printed every few timesteps to the screen and logfile.
* :doc:`Dump files <dump>`, which contain snapshots of atoms and various
  per-atom values and are written at a specified frequency.
* Certain fixes can output user-specified quantities to files: :doc:`fix ave/time <fix_ave_time>` for time averaging, :doc:`fix ave/chunk <fix_ave_chunk>` for spatial or other averaging, and :doc:`fix print <fix_print>` for single-line output of
  :doc:`variables <variable>`.  Fix print can also output to the
  screen.
* :doc:`Restart files <restart>`.


A simulation prints one set of thermodynamic output and (optionally)
restart files.  It can generate any number of dump files and fix
output files, depending on what :doc:`dump <dump>` and :doc:`fix <fix>`
commands you specify.

As discussed below, LAMMPS gives you a variety of ways to determine
what quantities are computed and printed when the thermodynamics,
dump, or fix commands listed above perform output.  Throughout this
discussion, note that users can also :doc:`add their own computes and fixes to LAMMPS <Modify>` which can then generate values that can then be
output with these commands.

The following sub-sections discuss different LAMMPS command related
to output and the kind of data they operate on and produce:

* :ref:`Global/per-atom/local data <global>`
* :ref:`Scalar/vector/array data <scalar>`
* :ref:`Thermodynamic output <thermo>`
* :ref:`Dump file output <dump>`
* :ref:`Fixes that write output files <fixoutput>`
* :ref:`Computes that process output quantities <computeoutput>`
* :ref:`Fixes that process output quantities <fixprocoutput>`
* :ref:`Computes that generate values to output <compute>`
* :ref:`Fixes that generate values to output <fix>`
* :ref:`Variables that generate values to output <variable>`
* :ref:`Summary table of output options and data flow between commands <table>`

.. _global:

Global/per-atom/local data
---------------------------------------

Various output-related commands work with three different styles of
data: global, per-atom, or local.  A global datum is one or more
system-wide values, e.g. the temperature of the system.  A per-atom
datum is one or more values per atom, e.g. the kinetic energy of each
atom.  Local datums are calculated by each processor based on the
atoms it owns, but there may be zero or more per atom, e.g. a list of
bond distances.

.. _scalar:

Scalar/vector/array data
-------------------------------------

Global, per-atom, and local datums can each come in three kinds: a
single scalar value, a vector of values, or a 2d array of values.  The
doc page for a "compute" or "fix" or "variable" that generates data
will specify both the style and kind of data it produces, e.g. a
per-atom vector.

When a quantity is accessed, as in many of the output commands
discussed below, it can be referenced via the following bracket
notation, where ID in this case is the ID of a compute.  The leading
"c\_" would be replaced by "f\_" for a fix, or "v\_" for a variable:

+-------------+--------------------------------------------+
| c\_ID       | entire scalar, vector, or array            |
+-------------+--------------------------------------------+
| c\_ID[I]    | one element of vector, one column of array |
+-------------+--------------------------------------------+
| c\_ID[I][J] | one element of array                       |
+-------------+--------------------------------------------+

In other words, using one bracket reduces the dimension of the data
once (vector -> scalar, array -> vector).  Using two brackets reduces
the dimension twice (array -> scalar).  Thus a command that uses
scalar values as input can typically also process elements of a vector
or array.

.. _thermo:

Thermodynamic output
---------------------------------

The frequency and format of thermodynamic output is set by the
:doc:`thermo <thermo>`, :doc:`thermo_style <thermo_style>`, and
:doc:`thermo_modify <thermo_modify>` commands.  The
:doc:`thermo_style <thermo_style>` command also specifies what values
are calculated and written out.  Pre-defined keywords can be specified
(e.g. press, etotal, etc).  Three additional kinds of keywords can
also be specified (c\_ID, f\_ID, v\_name), where a :doc:`compute <compute>`
or :doc:`fix <fix>` or :doc:`variable <variable>` provides the value to be
output.  In each case, the compute, fix, or variable must generate
global values for input to the :doc:`thermo_style custom <dump>`
command.

Note that thermodynamic output values can be "extensive" or
"intensive".  The former scale with the number of atoms in the system
(e.g. total energy), the latter do not (e.g. temperature).  The
setting for :doc:`thermo_modify norm <thermo_modify>` determines whether
extensive quantities are normalized or not.  Computes and fixes
produce either extensive or intensive values; see their individual doc
pages for details.  :doc:`Equal-style variables <variable>` produce only
intensive values; you can include a division by "natoms" in the
formula if desired, to make an extensive calculation produce an
intensive result.

.. _dump:

Dump file output
---------------------------

Dump file output is specified by the :doc:`dump <dump>` and
:doc:`dump_modify <dump_modify>` commands.  There are several
pre-defined formats (dump atom, dump xtc, etc).

There is also a :doc:`dump custom <dump>` format where the user
specifies what values are output with each atom.  Pre-defined atom
attributes can be specified (id, x, fx, etc).  Three additional kinds
of keywords can also be specified (c\_ID, f\_ID, v\_name), where a
:doc:`compute <compute>` or :doc:`fix <fix>` or :doc:`variable <variable>`
provides the values to be output.  In each case, the compute, fix, or
variable must generate per-atom values for input to the :doc:`dump custom <dump>` command.

There is also a :doc:`dump local <dump>` format where the user specifies
what local values to output.  A pre-defined index keyword can be
specified to enumerate the local values.  Two additional kinds of
keywords can also be specified (c\_ID, f\_ID), where a
:doc:`compute <compute>` or :doc:`fix <fix>` or :doc:`variable <variable>`
provides the values to be output.  In each case, the compute or fix
must generate local values for input to the :doc:`dump local <dump>`
command.

.. _fixoutput:

Fixes that write output files
---------------------------------------------

Several fixes take various quantities as input and can write output
files: :doc:`fix ave/time <fix_ave_time>`, :doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix ave/histo <fix_ave_histo>`,
:doc:`fix ave/correlate <fix_ave_correlate>`, and :doc:`fix print <fix_print>`.

The :doc:`fix ave/time <fix_ave_time>` command enables direct output to
a file and/or time-averaging of global scalars or vectors.  The user
specifies one or more quantities as input.  These can be global
:doc:`compute <compute>` values, global :doc:`fix <fix>` values, or
:doc:`variables <variable>` of any style except the atom style which
produces per-atom values.  Since a variable can refer to keywords used
by the :doc:`thermo_style custom <thermo_style>` command (like temp or
press) and individual per-atom values, a wide variety of quantities
can be time averaged and/or output in this way.  If the inputs are one
or more scalar values, then the fix generate a global scalar or vector
of output.  If the inputs are one or more vector values, then the fix
generates a global vector or array of output.  The time-averaged
output of this fix can also be used as input to other output commands.

The :doc:`fix ave/chunk <fix_ave_chunk>` command enables direct output
to a file of chunk-averaged per-atom quantities like those output in
dump files.  Chunks can represent spatial bins or other collections of
atoms, e.g. individual molecules.  The per-atom quantities can be atom
density (mass or number) or atom attributes such as position,
velocity, force.  They can also be per-atom quantities calculated by a
:doc:`compute <compute>`, by a :doc:`fix <fix>`, or by an atom-style
:doc:`variable <variable>`.  The chunk-averaged output of this fix can
also be used as input to other output commands.

The :doc:`fix ave/histo <fix_ave_histo>` command enables direct output
to a file of histogrammed quantities, which can be global or per-atom
or local quantities.  The histogram output of this fix can also be
used as input to other output commands.

The :doc:`fix ave/correlate <fix_ave_correlate>` command enables direct
output to a file of time-correlated quantities, which can be global
values.  The correlation matrix output of this fix can also be used as
input to other output commands.

The :doc:`fix print <fix_print>` command can generate a line of output
written to the screen and log file or to a separate file, periodically
during a running simulation.  The line can contain one or more
:doc:`variable <variable>` values for any style variable except the
vector or atom styles).  As explained above, variables themselves can
contain references to global values generated by :doc:`thermodynamic keywords <thermo_style>`, :doc:`computes <compute>`,
:doc:`fixes <fix>`, or other :doc:`variables <variable>`, or to per-atom
values for a specific atom.  Thus the :doc:`fix print <fix_print>`
command is a means to output a wide variety of quantities separate
from normal thermodynamic or dump file output.

.. _computeoutput:

Computes that process output quantities
-----------------------------------------------------------

The :doc:`compute reduce <compute_reduce>` and :doc:`compute reduce/region <compute_reduce>` commands take one or more per-atom
or local vector quantities as inputs and "reduce" them (sum, min, max,
ave) to scalar quantities.  These are produced as output values which
can be used as input to other output commands.

The :doc:`compute slice <compute_slice>` command take one or more global
vector or array quantities as inputs and extracts a subset of their
values to create a new vector or array.  These are produced as output
values which can be used as input to other output commands.

The :doc:`compute property/atom <compute_property_atom>` command takes a
list of one or more pre-defined atom attributes (id, x, fx, etc) and
stores the values in a per-atom vector or array.  These are produced
as output values which can be used as input to other output commands.
The list of atom attributes is the same as for the :doc:`dump custom <dump>` command.

The :doc:`compute property/local <compute_property_local>` command takes
a list of one or more pre-defined local attributes (bond info, angle
info, etc) and stores the values in a local vector or array.  These
are produced as output values which can be used as input to other
output commands.

.. _fixprocoutput:

Fixes that process output quantities
--------------------------------------------------------

The :doc:`fix vector <fix_vector>` command can create global vectors as
output from global scalars as input, accumulating them one element at
a time.

The :doc:`fix ave/atom <fix_ave_atom>` command performs time-averaging
of per-atom vectors.  The per-atom quantities can be atom attributes
such as position, velocity, force.  They can also be per-atom
quantities calculated by a :doc:`compute <compute>`, by a
:doc:`fix <fix>`, or by an atom-style :doc:`variable <variable>`.  The
time-averaged per-atom output of this fix can be used as input to
other output commands.

The :doc:`fix store/state <fix_store_state>` command can archive one or
more per-atom attributes at a particular time, so that the old values
can be used in a future calculation or output.  The list of atom
attributes is the same as for the :doc:`dump custom <dump>` command,
including per-atom quantities calculated by a :doc:`compute <compute>`,
by a :doc:`fix <fix>`, or by an atom-style :doc:`variable <variable>`.
The output of this fix can be used as input to other output commands.

.. _compute:

Computes that generate values to output
-----------------------------------------------------

Every :doc:`compute <compute>` in LAMMPS produces either global or
per-atom or local values.  The values can be scalars or vectors or
arrays of data.  These values can be output using the other commands
described in this section.  The doc page for each compute command
describes what it produces.  Computes that produce per-atom or local
values have the word "atom" or "local" in their style name.  Computes
without the word "atom" or "local" produce global values.

.. _fix:

Fixes that generate values to output
----------------------------------------------

Some :doc:`fixes <fix>` in LAMMPS produces either global or per-atom or
local values which can be accessed by other commands.  The values can
be scalars or vectors or arrays of data.  These values can be output
using the other commands described in this section.  The doc page for
each fix command tells whether it produces any output quantities and
describes them.

.. _variable:

Variables that generate values to output
-------------------------------------------------------

:doc:`Variables <variable>` defined in an input script can store one or
more strings.  But equal-style, vector-style, and atom-style or
atomfile-style variables generate a global scalar value, global vector
or values, or a per-atom vector, respectively, when accessed.  The
formulas used to define these variables can contain references to the
thermodynamic keywords and to global and per-atom data generated by
computes, fixes, and other variables.  The values generated by
variables can be used as input to and thus output by the other
commands described in this section.

.. _table:

Summary table of output options and data flow between commands
--------------------------------------------------------------------------

This table summarizes the various commands that can be used for
generating output from LAMMPS.  Each command produces output data of
some kind and/or writes data to a file.  Most of the commands can take
data from other commands as input.  Thus you can link many of these
commands together in pipeline form, where data produced by one command
is used as input to another command and eventually written to the
screen or to a file.  Note that to hook two commands together the
output and input data types must match, e.g. global/per-atom/local
data and scalar/vector/array data.

Also note that, as described above, when a command takes a scalar as
input, that could be an element of a vector or array.  Likewise a
vector input could be a column of an array.

+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| Command                                                | Input                                        | Output                                    |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`thermo_style custom <thermo_style>`             | global scalars                               | screen, log file                          |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`dump custom <dump>`                              | per-atom vectors                             | dump file                                 |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`dump local <dump>`                               | local vectors                                | dump file                                 |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix print <fix_print>`                           | global scalar from variable                  | screen, file                              |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`print <print>`                                   | global scalar from variable                  | screen                                    |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`computes <compute>`                              | N/A                                          | global/per-atom/local scalar/vector/array |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fixes <fix>`                                     | N/A                                          | global/per-atom/local scalar/vector/array |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`variables <variable>`                            | global scalars and vectors, per-atom vectors | global scalar and vector, per-atom vector |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`compute reduce <compute_reduce>`                 | per-atom/local vectors                       | global scalar/vector                      |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`compute slice <compute_slice>`                   | global vectors/arrays                        | global vector/array                       |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`compute property/atom <compute_property_atom>`   | per-atom vectors                             | per-atom vector/array                     |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`compute property/local <compute_property_local>` | local vectors                                | local vector/array                        |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix vector <fix_vector>`                         | global scalars                               | global vector                             |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix ave/atom <fix_ave_atom>`                     | per-atom vectors                             | per-atom vector/array                     |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix ave/time <fix_ave_time>`                     | global scalars/vectors                       | global scalar/vector/array, file          |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix ave/chunk <fix_ave_chunk>`                   | per-atom vectors                             | global array, file                        |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix ave/histo <fix_ave_histo>`                   | global/per-atom/local scalars and vectors    | global array, file                        |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix ave/correlate <fix_ave_correlate>`           | global scalars                               | global array, file                        |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+
| :doc:`fix store/state <fix_store_state>`               | per-atom vectors                             | per-atom vector/array                     |
+--------------------------------------------------------+----------------------------------------------+-------------------------------------------+


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
