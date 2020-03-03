.. index:: thermo_modify

thermo_modify command
=====================

Syntax
""""""

 
.. code-block:: LAMMPS

   thermo_modify keyword value ...

* one or more keyword/value pairs may be listed
  
  .. parsed-literal::
  
     keyword = *lost* or *lost/bond* or *norm* or *flush* or *line* or *format* or *temp* or *press*\ :l
       *lost* value = *error* or *warn* or *ignore*
       *lost/bond* value = *error* or *warn* or *ignore*
       *norm* value = *yes* or *no*
       *flush* value = *yes* or *no*
       *line* value = *one* or *multi*
       *format* values = *line* string, *int* string, *float* string, M string, or *none*
         string = C-style format string
         M = integer from 1 to N, where N = # of quantities being output
       *temp* value = compute ID that calculates a temperature
       *press* value = compute ID that calculates a pressure



Examples
""""""""


.. code-block:: LAMMPS

   thermo_modify lost ignore flush yes
   thermo_modify temp myTemp format 3 %15.8g
   thermo_modify temp myTemp format line "%ld %g %g %15.8g"
   thermo_modify line multi format float %g

Description
"""""""""""

Set options for how thermodynamic information is computed and printed
by LAMMPS.

.. note::

   These options apply to the currently defined thermo style.  When
   you specify a :doc:`thermo_style <thermo_style>` command, all
   thermodynamic settings are restored to their default values, including
   those previously reset by a thermo\_modify command.  Thus if your input
   script specifies a thermo\_style command, you should use the
   thermo\_modify command after it.

The *lost* keyword determines whether LAMMPS checks for lost atoms
each time it computes thermodynamics and what it does if atoms are
lost.  An atom can be "lost" if it moves across a non-periodic
simulation box :doc:`boundary <boundary>` or if it moves more than a box
length outside the simulation domain (or more than a processor
sub-domain length) before reneighboring occurs.  The latter case is
typically due to bad dynamics, e.g. too large a timestep or huge
forces and velocities.  If the value is *ignore*\ , LAMMPS does not
check for lost atoms.  If the value is *error* or *warn*\ , LAMMPS
checks and either issues an error or warning.  The code will exit with
an error and continue with a warning.  A warning will only be issued
once, the first time an atom is lost.  This can be a useful debugging
option.

The *lost/bond* keyword determines whether LAMMPS throws an error or
not if an atom in a bonded interaction (bond, angle, etc) cannot be
found when it creates bonded neighbor lists.  By default this is a
fatal error.  However in some scenarios it may be desirable to only
issue a warning or ignore it and skip the computation of the missing
bond, angle, etc.  An example would be when gas molecules in a vapor
are drifting out of the box through a fixed boundary condition (see
the :doc:`boundary <boundary>` command).  In this case one atom may be
deleted before the rest of the molecule is, on a later timestep.

The *norm* keyword determines whether various thermodynamic output
values are normalized by the number of atoms or not, depending on
whether it is set to *yes* or *no*\ .  Different unit styles have
different defaults for this setting (see below).  Even if *norm* is
set to *yes*\ , a value is only normalized if it is an "extensive"
quantity, meaning that it scales with the number of atoms in the
system.  For the thermo keywords described by the doc page for the
:doc:`thermo_style <thermo_style>` command, all energy-related keywords
are extensive, such as *pe* or *ebond* or *enthalpy*\ .  Other keywords
such as *temp* or *press* are "intensive" meaning their value is
independent (in a statistical sense) of the number of atoms in the
system and thus are never normalized.  For thermodynamic output values
extracted from fixes and computes in a :doc:`thermo_style custom <thermo_style>` command, the doc page for the individual
:doc:`fix <fix>` or :doc:`compute <compute>` lists whether the value is
"extensive" or "intensive" and thus whether it is normalized.
Thermodynamic output values calculated by a variable formula are
assumed to be "intensive" and thus are never normalized.  You can
always include a divide by the number of atoms in the variable formula
if this is not the case.

The *flush* keyword invokes a flush operation after thermodynamic info
is written to the log file.  This insures the output in that file is
current (no buffering by the OS), even if LAMMPS halts before the
simulation completes.

The *line* keyword determines whether thermodynamics will be output as
a series of numeric values on one line or in a multi-line format with
3 quantities with text strings per line and a dashed-line header
containing the timestep and CPU time.  This modify option overrides
the *one* and *multi* thermo\_style settings.

The *format* keyword can be used to change the default numeric format
of any of quantities the :doc:`thermo_style <thermo_style>` command
outputs.  All the specified format strings are C-style formats,
e.g. as used by the C/C++ printf() command.  The *line* keyword takes
a single argument which is the format string for the entire line of
thermo output, with N fields, which you must enclose in quotes if it
is more than one field.  The *int* and *float* keywords take a single
format argument and are applied to all integer or floating-point
quantities output.  The setting for *M string* also takes a single
format argument which is used for the Mth value output in each line,
e.g. the 5th column is output in high precision for "format 5
%20.15g".

The *format* keyword can be used multiple times.  The precedence is
that for each value in a line of output, the *M* format (if specified)
is used, else the *int* or *float* setting (if specified) is used,
else the *line* setting (if specified) for that value is used, else
the default setting is used.  A setting of *none* clears all previous
settings, reverting all values to their default format.

.. note::

   The thermo output values *step* and *atoms* are stored
   internally as 8-byte signed integers, rather than the usual 4-byte
   signed integers.  When specifying the *format int* option you can use
   a "%d"-style format identifier in the format string and LAMMPS will
   convert this to the corresponding 8-byte form when it is applied to
   those keywords.  However, when specifying the *line* option or *format
   M string* option for *step* and *natoms*\ , you should specify a format
   string appropriate for an 8-byte signed integer, e.g. one with "%ld".

The *temp* keyword is used to determine how thermodynamic temperature
is calculated, which is used by all thermo quantities that require a
temperature ("temp", "press", "ke", "etotal", "enthalpy", "pxx", etc).
The specified compute ID must have been previously defined by the user
via the :doc:`compute <compute>` command and it must be a style of
compute that calculates a temperature.  As described in the
:doc:`thermo_style <thermo_style>` command, thermo output uses a default
compute for temperature with ID = *thermo\_temp*.  This option allows
the user to override the default.

The *press* keyword is used to determine how thermodynamic pressure is
calculated, which is used by all thermo quantities that require a
pressure ("press", "enthalpy", "pxx", etc).  The specified compute ID
must have been previously defined by the user via the
:doc:`compute <compute>` command and it must be a style of compute that
calculates a pressure.  As described in the
:doc:`thermo_style <thermo_style>` command, thermo output uses a default
compute for pressure with ID = *thermo\_press*.  This option allows the
user to override the default.

.. note::

   If both the *temp* and *press* keywords are used in a single
   thermo\_modify command (or in two separate commands), then the order in
   which the keywords are specified is important.  Note that a :doc:`pressure compute <compute_pressure>` defines its own temperature compute as
   an argument when it is specified.  The *temp* keyword will override
   this (for the pressure compute being used by thermodynamics), but only
   if the *temp* keyword comes after the *press* keyword.  If the *temp*
   keyword comes before the *press* keyword, then the new pressure
   compute specified by the *press* keyword will be unaffected by the
   *temp* setting.

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`thermo <thermo>`, :doc:`thermo_style <thermo_style>`

Default
"""""""

The option defaults are lost = error, norm = yes for unit style of
*lj*\ , norm = no for unit style of *real* and *metal*\ , flush = no,
and temp/press = compute IDs defined by thermo\_style.

The defaults for the line and format options depend on the thermo
style.  For styles "one" and "custom", the line and format defaults
are "one", "%8d", and "%12.8g".  For style "multi", the line and
format defaults are "multi", "%8d", and "%14.4f".
