.. index:: fix ave/histo
.. index:: fix ave/histo/weight

fix ave/histo command
=====================

fix ave/histo/weight command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style Nevery Nrepeat Nfreq lo hi Nbin value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *ave/histo* or *ave/histo/weight* = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = # of times to use input values for calculating histogram
* Nfreq = calculate histogram every this many timesteps
* lo,hi = lo/hi bounds within which to histogram
* Nbin = # of histogram bins
* one or more input values can be listed
* value = *x*, *y*, *z*, *vx*, *vy*, *vz*, *fx*, *fy*, *fz*, c_ID, c_ID[N], f_ID, f_ID[N], v_name

  .. parsed-literal::

       x,y,z,vx,vy,vz,fx,fy,fz = atom attribute (position, velocity, force component)
       c_ID = scalar or vector calculated by a compute with ID
       c_ID[I] = Ith component of vector or Ith column of array calculated by a compute with ID, I can include wildcard (see below)
       f_ID = scalar or vector calculated by a fix with ID
       f_ID[I] = Ith component of vector or Ith column of array calculated by a fix with ID, I can include wildcard (see below)
       v_name = value(s) calculated by an equal-style or vector-style or atom-style variable with name
       v_name[I] = value calculated by a vector-style variable with name, I can include wildcard (see below)

* zero or more keyword/arg pairs may be appended
* keyword = *mode* or *kind* or *file* or *ave* or *start* or *beyond* or *overwrite* or *title1* or *title2* or *title3*

  .. parsed-literal::

       *mode* arg = *scalar* or *vector*
         scalar = all input values are scalars
         vector = all input values are vectors
       *kind* arg = *global* or *peratom* or *local*
       *file* arg = filename
         filename = name of file to output histogram(s) to
       *ave* args = *one* or *running* or *window*
         one = output a new average value every Nfreq steps
         running = output cumulative average of all previous Nfreq steps
         window M = output average of M most recent Nfreq steps
       *start* args = Nstart
         Nstart = start averaging on this timestep
       *beyond* arg = *ignore* or *end* or *extra*
         ignore = ignore values outside histogram lo/hi bounds
         end = count values outside histogram lo/hi bounds in end bins
         extra = create 2 extra bins for value outside histogram lo/hi bounds
       *overwrite* arg = none = overwrite output file with only latest output
       *title1* arg = string
         string = text to print as 1st line of output file
       *title2* arg = string
         string = text to print as 2nd line of output file
       *title3* arg = string
         string = text to print as 3rd line of output file, only for vector mode

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/histo 100 5 1000 0.5 1.5 50 c_myTemp file temp.histo ave running
   fix 1 all ave/histo 100 5 1000 -5 5 100 c_thermo_press[2] c_thermo_press[3] title1 "My output values"
   fix 1 all ave/histo 100 5 1000 -5 5 100 c_thermo_press[*]
   fix 1 all ave/histo 1 100 1000 -2.0 2.0 18 vx vy vz mode vector ave running beyond extra
   fix 1 all ave/histo/weight 1 1 1 10 100 2000 c_XRD[1] c_XRD[2]

Description
"""""""""""

Use one or more values as inputs every few timesteps to create a
single histogram.  The histogram can then be averaged over longer
timescales.  The resulting histogram can be used by other :doc:`output commands <Howto_output>`, and can also be written to a file.  The
fix ave/histo/weight command has identical syntax to fix ave/histo,
except that exactly two values must be specified.  See details below.

The group specified with this command is ignored for global and local
input values.  For per-atom input values, only atoms in the group
contribute to the histogram.  Note that regardless of the specified
group, specified values may represent calculations performed by
computes and fixes which store their own "group" definition.

A histogram is simply a count of the number of values that fall within
a histogram bin.  *Nbins* are defined, with even spacing between *lo*
and *hi*\ .  Values that fall outside the lo/hi bounds can be treated in
different ways; see the discussion of the *beyond* keyword below.

Each input value can be an atom attribute (position, velocity, force
component) or can be the result of a :doc:`compute <compute>` or
:doc:`fix <fix>` or the evaluation of an equal-style or vector-style or
atom-style :doc:`variable <variable>`.  The set of input values can be
either all global, all per-atom, or all local quantities.  Inputs of
different kinds (e.g. global and per-atom) cannot be mixed.  Atom
attributes are per-atom vector values.  See the page for
individual "compute" and "fix" commands to see what kinds of
quantities they generate.  See the optional *kind* keyword below for
how to force the fix ave/histo command to disambiguate if necessary.

Note that the output of this command is a single histogram for all
input values combined together, not one histogram per input value.
See below for details on the format of the output of this fix.

The input values must either be all scalars or all vectors (or
arrays), depending on the setting of the *mode* keyword.

If *mode* = scalar, then the input values must be scalars, or vectors
with a bracketed term appended, indicating the Ith value of the vector
is used.

If *mode* = vector, then the input values must be vectors, or arrays
with a bracketed term appended, indicating the Ith column of the array
is used.

If the fix ave/histo/weight command is used, exactly two values must
be specified.  If the values are vectors, they must be the same
length.  The first value (a scalar or vector) is what is histogrammed
into bins, in the same manner the fix ave/histo command operates.  The
second value (a scalar or vector) is used as a "weight".  This means
that instead of each value tallying a "1" to its bin, the
corresponding weight is tallied.  For example, the :math:`N^\text{th}` entry
(weight) in the second vector is tallied to the bin corresponding to the
:math:`N^\text{th}` entry in the first vector.

----------

For input values from a compute or fix or variable, the bracketed
index I can be specified using a wildcard asterisk with the index to
effectively specify multiple values.  This takes the form "\*" or
"\*n" or "m\*" or "m\*n".  If :math:`N` is the size of the vector
(for *mode* = scalar) or the number of columns in the array
(for *mode* = vector), then an asterisk with no numeric values means all
indices from 1 to :math:`N`\ .
A leading asterisk means all indices from 1 to n (inclusive).  A
trailing asterisk means all indices from m to :math:`N` (inclusive).  A middle
asterisk means all indices from m to n (inclusive).

Using a wildcard is the same as if the individual elements of the
vector or columns of the array had been listed one by one.  For example, the
following two fix ave/histo commands are equivalent, since the :doc:`compute
com/chunk <compute_com_chunk>` command creates a global array with three
columns:

.. code-block:: LAMMPS

   compute myCOM all com/chunk
   fix 1 all ave/histo 100 1 100 -10.0 10.0 100 c_myCOM[*] file tmp1.com mode vector
   fix 2 all ave/histo 100 1 100 -10.0 10.0 100 c_myCOM[1] c_myCOM[2] c_myCOM[3] file tmp2.com mode vector

.. note::

   For a vector-style variable, only the wildcard forms "\*n" or
   "m\*n" are allowed.  You must specify the upper bound, because
   vector-style variable lengths are not determined until the variable
   is evaluated.  If n is specified larger than the vector length
   turns out to be, zeroes are output for missing vector values.

----------

The :math:`N_\text{every}`, :math:`N_\text{repeat}`, and :math:`N_\text{freq}`
arguments specify on what time steps the input values will be used in order to
contribute to the histogram.  The final histogram is generated on time steps
that are multiple of :math:`N_\text{freq}`\ .  It is averaged over
:math:`N_\text{repeat}` histograms, computed in the preceding portion of the
simulation every :math:`N_\text{every}` time steps.
:math:`N_\text{freq}` must be a multiple of :math:`N_\text{every}` and
:math:`N_\text{every}` must be non-zero even if :math:`N_\text{repeat}` is 1.
Also, the time steps contributing to the histogram value cannot overlap
(i.e., :math:`N_\text{repeat}\times N_\text{every}` cannot exceed
:math:`N_\text{freq}`).

For example, if :math:`N_\text{every}=2`, :math:`N_\text{repeat}=6`, and
:math:`N_\text{freq}=100`, then input values on time steps 90, 92, 94, 96, 98,
and 100 will be used to compute the final histogram on timestep 100.
Similarly for timesteps 190, 192, 194, 196, 198, and 200 on timestep 200, etc.
If :math:`N_\text{repeat}=1` and :math:`N_\text{freq} = 100`, then no time
averaging of the histogram is done; a histogram is simply generated on
timesteps 100, 200, etc.

----------

The atom attribute values (*x*, *y*, *z*, *vx*, *vy*, *vz*, *fx*, *fy*, and
*fz*) are self-explanatory.  Note that other atom attributes can be used as
inputs to this fix by using the
:doc:`compute property/atom <compute_property_atom>` command and then
specifying an input value from that compute.

If a value begins with "c\_", a compute ID must follow which has been
previously defined in the input script.  If *mode* = scalar, then if
no bracketed term is appended, the global scalar calculated by the
compute is used.  If a bracketed term is appended, the Ith element of
the global vector calculated by the compute is used.  If *mode* =
vector, then if no bracketed term is appended, the global or per-atom
or local vector calculated by the compute is used.  If a bracketed
term is appended, the Ith column of the global or per-atom or local
array calculated by the compute is used.  See the discussion above for
how I can be specified with a wildcard asterisk to effectively specify
multiple values.

Note that there is a :doc:`compute reduce <compute_reduce>` command
that can sum per-atom quantities into a global scalar or vector, which
can then be accessed by fix ave/histo.  It can also be a compute defined
not in your input script, but by :doc:`thermodynamic output <thermo_style>`
or other fixes such as :doc:`fix nvt <fix_nh>`
or :doc:`fix temp/rescale <fix_temp_rescale>`.  See the doc pages for
these commands which give the IDs of these computes.  Users can also
write code for their own compute styles and
:doc:`add them to LAMMPS <Modify>`.

If a value begins with "f\_", a fix ID must follow which has been
previously defined in the input script.  If *mode* = scalar, then if
no bracketed term is appended, the global scalar calculated by the fix
is used.  If a bracketed term is appended, the Ith element of the
global vector calculated by the fix is used.  If *mode* = vector, then
if no bracketed term is appended, the global or per-atom or local
vector calculated by the fix is used.  If a bracketed term is
appended, the :math:`I^\text{th}` column of the global or per-atom or local
array calculated by the fix is used.  See the discussion above for how
:math:`I` can be specified with a wildcard asterisk to effectively specify
multiple values.

Note that some fixes only produce their values on certain timesteps,
which must be compatible with :math:`N_\text{every}`, else an error will
result.  Users can also write code for their own fix styles and
:doc:`add them to LAMMPS <Modify>`.

If a value begins with "v\_", a variable name must follow which has
been previously defined in the input script.  If *mode* = scalar, then
only equal-style or vector-style variables can be used, which both
produce global values.  In this mode, a vector-style variable requires
a bracketed term to specify the :math:`I^\text{th}` element of the vector
calculated by the variable.  If *mode* = vector, then only vector-style or
atom-style variables can be used, which produce a global or per-atom
vector respectively.  The vector-style variable must be used without a
bracketed term.  See the :doc:`variable <variable>` command for details.

Note that variables of style *equal*, *vector*, and *atom* define a
formula which can reference individual atom properties or
thermodynamic keywords, or they can invoke other computes, fixes, or
variables when they are evaluated, so this is a very general means of
specifying quantities to histogram.

----------

Additional optional keywords also affect the operation of this fix.

If the *mode* keyword is set to *scalar*, then all input values must
be global scalars, or elements of global vectors.  If the *mode*
keyword is set to *vector*, then all input values must be global or
per-atom or local vectors, or columns of global or per-atom or local
arrays.

The *kind* keyword only needs to be set if a compute or fix produces
more than one kind of output (global, per-atom, local).  If this is
not the case, then LAMMPS will determine what kind of input is
provided and whether all the input arguments are consistent.  If a
compute or fix produces more than one kind of output, the *kind*
keyword should be used to specify which output will be used.  The
remaining input arguments must still be consistent.

The *beyond* keyword determines how input values that fall outside the
*lo* to *hi* bounds are treated.  Values such that *lo* :math:`\le` value
:math:`\le` *hi* are assigned to one bin.  Values on a bin boundary are
assigned to the lower of the two bins.  If *beyond* is set to *ignore* then
values :math:`<` *lo* and values :math:`>` *hi* are ignored (i.e., they are not
binned). If *beyond* is set to *end*, then values :math:`<` *lo* are counted in
the first bin and values :math:`>` *hi* are counted in the last bin.
If *beyond* is set to *extend*, then two extra bins are created so that there
are :math:`N_\text{bins}+2` total bins.  Values :math:`<` *lo* are counted in
the first bin and values :math:`>` *hi* are counted in the last bin
:math:`(N_\text{bins}+2)`\ .  Values between
*lo* and *hi* (inclusive) are counted in bins 2 through
:math:`N_\text{bins}+1`\ .  The "coordinate" stored and printed for these two
extra bins is *lo* and *hi*\ .

The *ave* keyword determines how the histogram produced every
:math:`N_\text{freq}` steps are averaged with histograms produced on previous
steps that were multiples of :math:`N_\text{freq}`, before they are accessed by
another output command or written to a file.

If the *ave* setting is *one*, then the histograms produced on
timesteps that are multiples of :math:`N_\text{freq}` are independent of each
other; they are output as-is without further averaging.

If the *ave* setting is *running*, then the histograms produced on
timesteps that are multiples of :math:`N_\text{freq}` are summed and averaged
in a cumulative sense before being output.  Each bin value in the histogram
is thus the average of the bin value produced on that timestep with all
preceding values for the same bin.  This running average begins when the fix is
defined; it can only be restarted by deleting the fix via the
:doc:`unfix <unfix>` command, or by re-defining the fix by re-specifying it.

If the *ave* setting is *window*, then the histograms produced on
timesteps that are multiples of :math:`N_\text{freq}` are summed within a
moving "window" of time, so that the last :math:`M` histograms are used to
produce the output (e.g., if :math:`M = 3` and :math:`N_\text{freq} = 1000`,
then the output on step 10000 will be the combined histogram of the individual
histograms on steps 8000, 9000, and 10000.  Outputs on early steps will be sums
over less than :math:`M` histograms if they are not available.

The *start* keyword specifies what timestep histogramming will begin
on.  The default is step 0.  Often input values can be 0.0 at time 0,
so setting *start* to a larger value can avoid including a 0.0 in
a running or windowed histogram.

The *file* keyword allows a filename to be specified.  Every *Nfreq*
steps, one histogram is written to the file.  This includes a leading
line that contains the timestep, number of bins, the total count of
values contributing to the histogram, the count of values that were
not histogrammed (see the *beyond* keyword), the minimum value
encountered, and the maximum value encountered.  The min/max values
include values that were not histogrammed.  Following the leading
line, one line per bin is written into the file.  Each line contains
the bin #, the coordinate for the center of the bin (between *lo* and
*hi*\ ), the count of values in the bin, and the normalized count.  The
normalized count is the bin count divided by the total count (not
including values not histogrammed), so that the normalized values sum
to 1.0 across all bins.

The *overwrite* keyword will continuously overwrite the output file
with the latest output, so that it only contains one timestep worth of
output.  This option can only be used with the *ave running* setting.

The *title1*, *title2*, and *title3* keywords allow specification of
the strings that will be printed as the first three lines of the output
file, assuming the *file* keyword was used.  LAMMPS uses default
values for each of these, so they do not need to be specified.

By default, these header lines are as follows:

.. parsed-literal::

   # Histogram for fix ID
   # TimeStep Number-of-bins Total-counts Missing-counts Min-value Max-value
   # Bin Coord Count Count/Total

In the first line, ID is replaced with the fix-ID.  The second line
describes the six values that are printed at the first of each section
of output.  The third describes the four values printed for each bin in
the histogram.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to
:doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.

This fix produces a global vector and global array which can be
accessed by various :doc:`output commands <Howto_output>`.  The values
can only be accessed on timesteps that are multiples of :math:`N_\text{freq}`
since that is when a histogram is generated.  The global vector has four
values:

* 1 = total counts in the histogram
* 2 = values that were not histogrammed (see *beyond* keyword)
* 3 = min value of all input values, including ones not histogrammed
* 4 = max value of all input values, including ones not histogrammed

The global array has :math:`N_\text{bins}` rows and three columns.  The
first column has the bin coordinate, the second column has the count of
values in that histogram bin, and the third column has the bin count
divided by the total count (not including missing counts), so that the
values in the third column sum to 1.0.

The vector and array values calculated by this fix are all treated as
intensive.  If this is not the case (e.g., due to histogramming
per-atom input values), then you will need to account for that when
interpreting the values produced by this fix.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.
This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`compute <compute>`, :doc:`fix ave/atom <fix_ave_atom>`,
:doc:`fix ave/chunk <fix_ave_chunk>`, :doc:`fix ave/time <fix_ave_time>`,
:doc:`variable <variable>`, :doc:`fix ave/correlate <fix_ave_correlate>`,

Default
"""""""

none

The option defaults are mode = scalar, kind = figured out from input
arguments, ave = one, start = 0, no file output, beyond = ignore, and
title 1,2,3 = strings as described above.
