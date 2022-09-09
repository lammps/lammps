.. index:: thermo_modify

thermo_modify command
=====================

Syntax
""""""

.. code-block:: LAMMPS

   thermo_modify keyword value ...

* one or more keyword/value pairs may be listed
* keyword = *lost* or *lost/bond* or *warn* or *norm* or *flush* or *line* or *colname* or *format* or *temp* or *press*

  .. parsed-literal::

       *lost* value = *error* or *warn* or *ignore*
       *lost/bond* value = *error* or *warn* or *ignore*
       *warn* value = *ignore* or *reset* or *default* or a number
       *norm* value = *yes* or *no*
       *flush* value = *yes* or *no*
       *line* value = *one* or *multi* or *yaml*
       *colname* values =  ID string, or *default*
         string = new column header name
         ID = integer from 1 to N, or integer from -1 to -N, where N = # of quantities being output
              *or* a thermo keyword or reference to compute, fix, property or variable.
       *format* values = *line* string, *int* string, *float* string, ID string, or *none*
         string = C-style format string
         ID = integer from 1 to N, or integer from -1 to -N, where N = # of quantities being output
              *or* a thermo keyword or reference to compute, fix, property or variable.
       *temp* value = compute ID that calculates a temperature
       *press* value = compute ID that calculates a pressure

Examples
""""""""

.. code-block:: LAMMPS

   thermo_modify lost ignore flush yes
   thermo_modify temp myTemp format 3 %15.8g
   thermo_modify temp myTemp format line "%ld %g %g %15.8g"
   thermo_modify line multi format float %g
   thermo_modify line yaml format none
   thermo_modify colname 1 Timestep colname -2 Pressure colname f_1[1] AvgDensity

Description
"""""""""""

Set options for how thermodynamic information is computed and printed by
LAMMPS.

.. note::

   These options apply to the *currently defined* thermo style.  When
   you specify a :doc:`thermo_style <thermo_style>` command, all
   thermodynamic settings are restored to their default values,
   including those previously reset by a thermo_modify command.  Thus if
   your input script specifies a thermo_style command, you should use
   the thermo_modify command **after** it.

The *lost* keyword determines whether LAMMPS checks for lost atoms each
time it computes thermodynamics and what it does if atoms are lost.  An
atom can be "lost" if it moves across a non-periodic simulation box
:doc:`boundary <boundary>` or if it moves more than a box length outside
the simulation domain (or more than a processor sub-domain length)
before reneighboring occurs.  The latter case is typically due to bad
dynamics, e.g. too large a timestep or huge forces and velocities.  If
the value is *ignore*, LAMMPS does not check for lost atoms.  If the
value is *error* or *warn*, LAMMPS checks and either issues an error or
warning.  The code will exit with an error and continue with a warning.
A warning will only be issued once, the first time an atom is lost.
This can be a useful debugging option.

The *lost/bond* keyword determines whether LAMMPS throws an error or
not if an atom in a bonded interaction (bond, angle, etc) cannot be
found when it creates bonded neighbor lists.  By default this is a
fatal error.  However in some scenarios it may be desirable to only
issue a warning or ignore it and skip the computation of the missing
bond, angle, etc.  An example would be when gas molecules in a vapor
are drifting out of the box through a fixed boundary condition (see
the :doc:`boundary <boundary>` command).  In this case one atom may be
deleted before the rest of the molecule is, on a later timestep.

The *warn* keyword allows you to control whether LAMMPS will print
warning messages and how many of them.  Most warning messages are only
printed by MPI rank 0.  They are usually pointing out important issues
that should be investigated, but LAMMPS cannot determine for
certain whether they are an indication of an error.

Some warning messages are printed during a run (or immediately before)
each time a specific MPI rank encounters the issue, e.g. bonds that are
stretched too far or dihedrals in extreme configurations. These number
of these can quickly blow up the size of the log file and screen output.
Thus a limit of 100 warning messages is applied by default.  The warning
count is applied to the entire input unless reset with a ``thermo_modify
warn reset`` command.  If there are more warnings than the limit, LAMMPS
will print one final warning that it will not print any additional
warning messages.

.. note::

   The warning limit is enforced on either the per-processor count or
   the total count across all processors. For efficiency reasons,
   however, the total count is only updated at steps with thermodynamic
   output. Thus when running on a large number of processors in
   parallel, the total number of warnings printed can be significantly
   larger than the given limit.

Any number after the keyword *warn* will change the warning limit
accordingly.  With the value *ignore* all warnings will be suppressed,
with the value *always* no limit will be applied and warnings will
always be printed, with the value *reset* the internal warning counter
will be reset to zero, and with the value *default*, the counter is
reset and the limit set to 100.  An example usage of either *reset* or
*default* would be to re-enable warnings that were disabled or have
reached the limit during equilibration, where the warnings would be
acceptable while the system is still adjusting, but then change to all
warnings for the production run, where they would indicate problems that
would require a closer look at what is causing them.

The *norm* keyword determines whether various thermodynamic output
values are normalized by the number of atoms or not, depending on
whether it is set to *yes* or *no*\ .  Different unit styles have
different defaults for this setting (see below).  Even if *norm* is set
to *yes*, a value is only normalized if it is an "extensive" quantity,
meaning that it scales with the number of atoms in the system.  For the
thermo keywords described by the page for the :doc:`thermo_style
<thermo_style>` command, all energy-related keywords are extensive, such
as *pe* or *ebond* or *enthalpy*\ .  Other keywords such as *temp* or
*press* are "intensive" meaning their value is independent (in a
statistical sense) of the number of atoms in the system and thus are
never normalized.  For thermodynamic output values extracted from fixes
and computes in a :doc:`thermo_style custom <thermo_style>` command, the
page for the individual :doc:`fix <fix>` or :doc:`compute <compute>`
lists whether the value is "extensive" or "intensive" and thus whether
it is normalized.  Thermodynamic output values calculated by a variable
formula are assumed to be "intensive" and thus are never normalized.
You can always include a divide by the number of atoms in the variable
formula if this is not the case.

The *flush* keyword invokes a flush operation after thermodynamic info
is written to the screen and log file.  This insures the output is
updated and not buffered (by the application) even if LAMMPS halts
before the simulation completes.  Please note that this does not affect
buffering by the OS or devices, so you may still lose data in case the
simulation stops due to a hardware failure.

The *line* keyword determines whether thermodynamics will be output as a
series of numeric values on one line ("one"), in a multi-line format
with 3 quantities with text strings per line and a dashed-line header
containing the timestep and CPU time ("multi"), or in a YAML format
block ("yaml").  This modify option overrides the *one*, *multi*, or
*yaml* thermo_style settings.

.. versionadded:: 4May2022

The *colname* keyword can be used to change the default header keyword
for a column or field of thermodynamic output.  The setting for *ID
string* replaces the default text with the provided string.  *ID* can be
a positive integer when it represents the column number counting from
the left, a negative integer when it represents the column number from
the right (i.e. -1 is the last column/keyword), or a thermo keyword (or
compute, fix, property, or variable reference) and then it replaces the
string for that specific thermo keyword.

The *colname* keyword can be used multiple times. If multiple *colname*
settings refer to the same keyword, the last setting has precedence.  A
setting of *default* clears all previous settings, reverting all values
to their default values.

The *format* keyword can be used to change the default numeric format of
any of quantities the :doc:`thermo_style <thermo_style>` command
outputs.  All the specified format strings are C-style formats, e.g. as
used by the C/C++ printf() command.  The *line* keyword takes a single
argument which is the format string for the entire line of thermo
output, with N fields, which you must enclose in quotes if it is more
than one field.  The *int* and *float* keywords take a single format
argument and are applied to all integer or floating-point quantities
output.  The setting for *ID string* also takes a single format argument
which is used for the indexed value in each line.  The interpretation is
the same as for *colname*, i.e. a positive integer is the n-th value
corresponding to the n-th thermo keyword, a negative integer is counting
backwards, and a string matches the entry with the thermo keyword.,
e.g. the fifth column is output in high precision for "format 5 %20.15g"
and the pair energy for "format epair %20.15g".

The *format* keyword can be used multiple times.  The precedence is
that for each value in a line of output, the *ID* format (if specified)
is used, else the *int* or *float* setting (if specified) is used,
else the *line* setting (if specified) for that value is used, else
the default setting is used.  A setting of *none* clears all previous
settings, reverting all values to their default format.

.. note::

   The thermo output values *step* and *atoms* are stored internally as
   8-byte signed integers, rather than the usual 4-byte signed integers.
   When specifying the *format int* option you can use a "%d"-style
   format identifier in the format string and LAMMPS will convert this
   to the corresponding 8-byte form when it is applied to those
   keywords.  However, when specifying the *line* option or *format ID
   string* option for *step* and *natoms*, you should specify a format
   string appropriate for an 8-byte signed integer, e.g. one with "%ld"
   or "%lld" depending on the platform.

The *temp* keyword is used to determine how thermodynamic temperature is
calculated, which is used by all thermo quantities that require a
temperature ("temp", "press", "ke", "etotal", "enthalpy", "pxx", etc).
The specified compute ID must have been previously defined by the user
via the :doc:`compute <compute>` command and it must be a style of
compute that calculates a temperature.  As described in the
:doc:`thermo_style <thermo_style>` command, thermo output uses a default
compute for temperature with ID = *thermo_temp*.  This option allows the
user to override the default.

The *press* keyword is used to determine how thermodynamic pressure is
calculated, which is used by all thermo quantities that require a
pressure ("press", "enthalpy", "pxx", etc).  The specified compute ID
must have been previously defined by the user via the :doc:`compute
<compute>` command and it must be a style of compute that calculates a
pressure.  As described in the :doc:`thermo_style <thermo_style>`
command, thermo output uses a default compute for pressure with ID =
*thermo_press*.  This option allows the user to override the default.

.. note::

   If both the *temp* and *press* keywords are used in a single
   thermo_modify command (or in two separate commands), then the order
   in which the keywords are specified is important.  Note that a
   :doc:`pressure compute <compute_pressure>` defines its own
   temperature compute as an argument when it is specified.  The *temp*
   keyword will override this (for the pressure compute being used by
   thermodynamics), but only if the *temp* keyword comes after the
   *press* keyword.  If the *temp* keyword comes before the *press*
   keyword, then the new pressure compute specified by the *press*
   keyword will be unaffected by the *temp* setting.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`thermo <thermo>`, :doc:`thermo_style <thermo_style>`

Default
"""""""

The option defaults are lost = error, warn = 100, norm = yes for unit
style of *lj*, norm = no for unit style of *real* and *metal*,
flush = no, and temp/press = compute IDs defined by thermo_style.

The defaults for the line and format options depend on the thermo style.
For styles "one" and "custom", the line and format defaults are "one",
"%10d", and "%14.8g".  For style "multi", the line and format defaults
are "multi", "%14d", and "%14.4f". For style "yaml", the line and format
defaults are "%d" and "%.15g".
