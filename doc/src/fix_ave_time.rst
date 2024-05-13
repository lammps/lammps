.. index:: fix ave/time

fix ave/time command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ave/time Nevery Nrepeat Nfreq value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/time = style name of this fix command
* Nevery = use input values every this many time steps
* Nrepeat = # of times to use input values for calculating averages
* Nfreq = calculate averages every this many time steps
* one or more input values can be listed
* value = c_ID, c_ID[N], f_ID, f_ID[N], v_name

  .. parsed-literal::

       c_ID = global scalar or vector calculated by a compute with ID
       c_ID[I] = Ith component of global vector or Ith column of global array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = global scalar or vector calculated by a fix with ID
       f_ID[I] = Ith component of global vector or Ith column of global array calculated by a fix with ID, I can include wildcard (see below)
       v_name = value(s) calculated by an equal-style or vector-style variable with name
       v_name[I] = value calculated by a vector-style variable with name, I can include wildcard (see below)

* zero or more keyword/arg pairs may be appended
* keyword = *mode* or *file* or *append* or *ave* or *start* or *off* or *overwrite* or *format* or *title1* or *title2* or *title3*

  .. parsed-literal::

       *mode* arg = *scalar* or *vector*
         scalar = all input values are global scalars
         vector = all input values are global vectors or global arrays
       *ave* args = *one* or *running* or *window M*
         one = output a new average value every Nfreq steps
         running = output cumulative average of all previous Nfreq steps
         window M = output average of M most recent Nfreq steps
       *start* args = Nstart
         Nstart = start averaging on this time step
       *off* arg = M = do not average this value
         M = value # from 1 to Nvalues
       *file* arg = filename
         filename = name of file to output time averages to
       *append* arg = filename
         filename = name of file to append time averages to
       *overwrite* arg = none = overwrite output file with only latest output
       *format* arg = string
         string = C-style format string
       *title1* arg = string
         string = text to print as 1st line of output file
       *title2* arg = string
         string = text to print as 2nd line of output file
       *title3* arg = string
         string = text to print as 3rd line of output file, only for vector mode

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/time 100 5 1000 c_myTemp c_thermo_temp file temp.profile
   fix 1 all ave/time 100 5 1000 c_thermo_press[2] ave window 20 &
                                 title1 "My output values"
   fix 1 all ave/time 100 5 1000 c_thermo_press[*]
   fix 1 all ave/time 1 100 1000 f_indent f_indent[1] file temp.indent off 1

Description
"""""""""""

Use one or more global values as inputs every few time steps, and
average them over longer timescales.  The resulting averages can be
used by other :doc:`output commands <Howto_output>` such as
:doc:`thermo_style custom <thermo_style>`, and can also be written to a
file.  Note that if no time averaging is done, this command can be
used as a convenient way to simply output one or more global values to
a file.

The group specified with this command is ignored.  However, note that
specified values may represent calculations performed by computes and
fixes which store their own "group" definitions.

Each listed value can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an equal-style or vector-style
:doc:`variable <variable>`.  In each case, the compute, fix, or variable
must produce a global quantity, not a per-atom or local quantity.  If
you wish to spatial- or time-average or histogram per-atom quantities
from a compute, fix, or variable, then see the
:doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix ave/atom <fix_ave_atom>`,
or :doc:`fix ave/histo <fix_ave_histo>` commands.  If you wish to sum a
per-atom quantity into a single global quantity, see the
:doc:`compute reduce <compute_reduce>` command.

:doc:`Computes <compute>` that produce global quantities are those which
do not have the word *atom* in their style name.  Only a few
:doc:`fixes <fix>` produce global quantities.  See the doc pages for
individual fixes for info on which ones produce such values.
:doc:`Variables <variable>` of style *equal* and *vector* are the only
ones that can be used with this fix.  Variables of style *atom* cannot
be used, since they produce per-atom values.

The input values must either be all scalars or all vectors depending
on the setting of the *mode* keyword.  In both cases, the averaging is
performed independently on each input value (i.e., each input scalar
is averaged independently or each element of each input vector is
averaged independently).

If *mode* = scalar, then the input values must be scalars, or vectors
with a bracketed term appended, indicating the :math:`I^\text{th}` value of the
vector is used.

If *mode* = vector, then the input values must be vectors, or arrays
with a bracketed term appended, indicating the Ith column of the array
is used.  All vectors must be the same length, which is the length of
the vector or number of rows in the array.

----------

For input values from a compute or fix or variable, the bracketed
index I can be specified using a wildcard asterisk with the index to
effectively specify multiple values.  This takes the form "\*" or
"\*n" or "m\*" or "m\*n".  If :math:`N` is the size of the vector (for *mode* =
scalar) or the number of columns in the array (for *mode* = vector),
then an asterisk with no numeric values means all indices from 1 to :math:`N`.
A leading asterisk means all indices from 1 to n (inclusive).  A trailing
asterisk means all indices from n to :math:`N` (inclusive).  A middle asterisk
means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual elements of the
vector or columns of the array had been listed one by one.  For example, the
following two fix ave/time commands are equivalent, since the :doc:`compute rdf
<compute_rdf>` command creates, in this case, a global array with three
columns, each of length 50:

.. code-block:: LAMMPS

   compute myRDF all rdf 50 1 2
   fix 1 all ave/time 100 1 100 c_myRDF[*] file tmp1.rdf mode vector
   fix 2 all ave/time 100 1 100 c_myRDF[1] c_myRDF[2] c_myRDF[3] file tmp2.rdf mode vector

.. note::

   For a vector-style variable, only the wildcard forms "\*n" or
   "m\*n" are allowed.  You must specify the upper bound, because
   vector-style variable lengths are not determined until the variable
   is evaluated.  If n is specified larger than the vector length
   turns out to be, zeroes are output for missing vector values.

----------

The :math:`N_\text{every}`, :math:`N_\text{repeat}`, and :math:`N_\text{freq}`
arguments specify on what time steps the input values will be used in order to
contribute to the average.  The final averaged quantities are generated on
time steps that are a multiple of :math:`N_\text{freq}`\ .  The average is over
:math:`N_\text{repeat}` quantities, computed in the preceding portion of the
simulation every :math:`N_\text{every}` time steps.  :math:`N_\text{freq}` must
be a multiple of :math:`N_\text{every}` and :math:`N_\text{every}` must be
non-zero even if :math:`N_\text{repeat} = 1`.  Also, the time steps
contributing to the average value cannot overlap,
i.e. Nrepeat\*Nevery can not exceed Nfreq.

For example, if :math:`N_\text{every}=2`, :math:`N_\text{repeat}=6`, and
:math:`N_\text{freq}=100`, then values on time steps 90, 92, 94, 96, 98, and
100 will be used to compute the final average on time step 100.  Similarly for
time steps 190, 192, 194, 196, 198, and 200 on time step 200, etc.
If :math:`N_\text{repeat}=1` and :math:`N_\text{freq} = 100`, then no time
averaging is done; values are simply generated on time steps
100, 200, etc.

----------

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If *mode* = scalar, then if
no bracketed term is appended, the global scalar calculated by the
compute is used.  If a bracketed term is appended, the Ith element of
the global vector calculated by the compute is used.  If *mode* =
vector, then if no bracketed term is appended, the global vector
calculated by the compute is used.  If a bracketed term is appended,
the Ith column of the global array calculated by the compute is used.
See the discussion above for how I can be specified with a wildcard
asterisk to effectively specify multiple values.

Note that there is a :doc:`compute reduce <compute_reduce>` command
that can sum per-atom quantities into a global scalar or vector, which
can then be accessed by fix ave/time.  It can also be a compute defined
not in your input script, but by :doc:`thermodynamic output
<thermo_style>` or other fixes such as :doc:`fix nvt <fix_nh>` or
:doc:`fix temp/rescale <fix_temp_rescale>`.  See the doc pages for
these commands which give the IDs of these computes.  Users can also
write code for their own compute styles and :doc:`add them to LAMMPS
<Modify>`.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If *mode* = scalar, then if
no bracketed term is appended, the global scalar calculated by the fix
is used.  If a bracketed term is appended, the Ith element of the
global vector calculated by the fix is used.  If *mode* = vector, then
if no bracketed term is appended, the global vector calculated by the
fix is used.  If a bracketed term is appended, the Ith column of the
global array calculated by the fix is used.  See the discussion above
for how I can be specified with a wildcard asterisk to effectively
specify multiple values.

Note that some fixes only produce their values on certain time steps,
which must be compatible with *Nevery*, else an error will result.
Users can also write code for their own fix styles and :doc:`add them to LAMMPS <Modify>`.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  If *mode* = scalar, then
only equal-style or vector-style variables can be used, which both
produce global values.  In this mode, a vector-style variable requires
a bracketed term to specify the Ith element of the vector calculated
by the variable.  If *mode* = vector, then only a vector-style
variable can be used, without a bracketed term.  See the
:doc:`variable <variable>` command for details.

Note that variables of style *equal* and *vector* define a formula
which can reference individual atom properties or thermodynamic
keywords, or they can invoke other computes, fixes, or variables when
they are evaluated, so this is a very general means of specifying
quantities to time average.

----------

Additional optional keywords also affect the operation of this fix.

If the *mode* keyword is set to *scalar*, then all input values must
be global scalars, or elements of global vectors.  If the *mode*
keyword is set to *vector*, then all input values must be global
vectors, or columns of global arrays.  They can also be global arrays,
which are converted into a series of global vectors (one per column),
as explained above.

The *ave* keyword determines how the values produced every
:math:`N_\text{freq}` steps are averaged with values produced on previous steps
that were multiples of :math:`N_\text{freq}`, before they are accessed by
another output command or written to a file.

If the *ave* setting is *one*, then the values produced on time steps
that are multiples of :math:`N_\text{freq}` are independent of each other; they
are output as-is without further averaging.

If the *ave* setting is *running*, then the values produced on
time steps that are multiples of :math:`N_\text{freq}` are summed and averaged
in a cumulative sense before being output.  Each output value is thus the
average of the value produced on that time step with all preceding
values.  This running average begins when the fix is defined; it can
only be restarted by deleting the fix via the :doc:`unfix <unfix>`
command, or by re-defining the fix by re-specifying it.

If the *ave* setting is *window*, then the values produced on
time steps that are multiples of *Nfreq* are summed and averaged within
a moving "window" of time, so that the last M values are used to
produce the output.  For example, if :math:`M = 3` and
:math:`N_\text{freq} = 1000`, then the output on step 10000 will be the average
of the individual values on steps 8000, 9000, and 10000.  Outputs on early
steps will average over less than :math:`M` values if they are not available.

The *start* keyword specifies what time step averaging will begin on.
The default is step 0.  Often input values can be 0.0 at time 0, so
setting *start* to a larger value can avoid including a 0.0 in a
running or windowed average.

The *off* keyword can be used to flag any of the input values.  If a
value is flagged, it will not be time averaged.  Instead the most
recent input value will always be stored and output.  This is useful
if one of more of the inputs produced by a compute or fix or variable
are effectively constant or are simply current values (e.g., they are
being written to a file with other time-averaged values for purposes
of creating well-formatted output).

.. versionadded:: 17Apr2024
   new keyword *append*

The *file* or *append* keywords allow a filename to be specified.  If
*file* is used, then the filename is overwritten if it already exists.
If *append* is used, then the filename is appended to if it already
exists, or created if it does not exist.  Every *Nfreq* steps, one
quantity or vector of quantities is written to the file for each input
value specified in the fix ave/time command.  For *mode* = scalar, this
means a single line is written each time output is performed.  Thus the
file ends up to be a series of lines, i.e. one column of numbers for
each input value.  For *mode* = vector, an array of numbers is written
each time output is performed.  The number of rows is the length of the
input vectors, and the number of columns is the number of values.  Thus
the file ends up to be a series of these array sections.

.. versionadded:: 4May2022

If the filename ends in '.yaml' or '.yml' then the output format
conforms to the `YAML standard <https://yaml.org/>`_ which allows
easy import that data into tools and scripts that support reading YAML
files. The :doc:`structured data Howto <Howto_structured_data>` contains
examples for parsing and plotting such data with very little programming
effort in Python using the *pyyaml*, *pandas*, and *matplotlib*
packages.

The *overwrite* keyword will continuously overwrite the output file
with the latest output, so that it only contains one time step worth of
output.  This option can only be used with the *ave running* setting.

The *format* keyword sets the numeric format of each value when it is
printed to a file via the *file* keyword.  Note that all values are
floating point quantities.  The default format is %g.  You can specify
a higher precision if desired (e.g., %20.16g).

The *title1* and *title2* and *title3* keywords allow specification of
the strings that will be printed as the first 2 or 3 lines of the
output file, assuming the *file* keyword was used.  LAMMPS uses
default values for each of these, so they do not need to be specified.

By default, these header lines are as follows for *mode* = scalar:

.. parsed-literal::

   # Time-averaged data for fix ID
   # TimeStep value1 value2 ...

In the first line, ID is replaced with the fix-ID.  In the second line
the values are replaced with the appropriate fields from the fix
ave/time command.  There is no third line in the header of the file,
so the *title3* setting is ignored when *mode* = scalar.

By default, these header lines are as follows for *mode* = vector:

.. parsed-literal::

   # Time-averaged data for fix ID
   # TimeStep Number-of-rows
   # Row value1 value2 ...

In the first line, ID is replaced with the fix-ID.  The second line
describes the two values that are printed at the first of each section
of output.  In the third line the values are replaced with the
appropriate fields from the fix ave/time command.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. versionadded:: 4May2022

No information about this fix is written to :doc:`binary restart files
<restart>`.  The :doc:`fix_modify colname <fix_modify>` option can be
used to change the name of the column in the output file.  When writing
a YAML format file this name will be in the list of keywords.

This fix produces a global scalar or global vector or global array
which can be accessed by various :doc:`output commands <Howto_output>`.
The values can only be accessed on time steps that are multiples of
:math:`N_\text{freq}` since that is when averaging is performed.

A scalar is produced if only a single input value is averaged and
*mode* = scalar.  A vector is produced if multiple input values are
averaged for *mode* = scalar, or a single input value for *mode* =
vector.  In the first case, the length of the vector is the number of
inputs.  In the second case, the length of the vector is the same as
the length of the input vector.  An array is produced if multiple
input values are averaged and *mode* = vector.  The global array has #
of rows = length of the input vectors and # of columns = number of
inputs.

If the fix produces a scalar or vector, then the scalar and each
element of the vector can be either "intensive" or "extensive",
depending on whether the values contributing to the scalar or vector
element are "intensive" or "extensive".  If the fix produces an array,
then all elements in the array must be the same, either "intensive" or
"extensive".  If a compute or fix provides the value being time
averaged, then the compute or fix determines whether the value is
intensive or extensive; see the page for that compute or fix for
further info.  Values produced by a variable are treated as intensive.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix ave/atom <fix_ave_atom>`,
:doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix ave/histo <fix_ave_histo>`,
:doc:`variable <variable>`, :doc:`fix ave/correlate <fix_ave_correlate>`,

Default
"""""""

The option defaults are mode = scalar, ave = one, start = 0, no file
output, format = %g, title 1,2,3 = strings as described above, and no
off settings for any input values.
